#!/usr/bin/env python3
"""Run the prospective Gate-A archive-root benchmark without reducing its verdict.

The command adapter at the top of this file is the only place that assumes the not-yet-landed
archive-root CLI.  A final/promotable run requires `--interface-status final`, structured
collection-build output, and at least one adversarial fixture hook.  The driver refuses to
overwrite a run directory and writes a complete SHA-256 artifact manifest only after every check
and command succeeds.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shutil
import socket
import subprocess
import tempfile
import tomllib
from pathlib import Path
from typing import Any, Mapping, Sequence

from archive_root_gate_common import (
    GateError,
    compare_v1_v2,
    load_json,
    parse_archive,
    require,
    run_timed,
    self_test_blake3,
    sha256_file,
    tree_digest,
    write_artifact_manifest,
    write_json,
)


DATE = "2026-09-01"
ROOT_DOMAIN = b"gravlax-aie-directory-root-v2\0"
SAMPLES = tuple("ABCDEFGH")
SOURCE = {
    "A": (110_183_760, "029de28423c59a8b9ca283e73687a9501cd7533cfd508788b61b15e79b751463"),
    "B": (83_674_178, "4feced152203593f22ddb5a13207dbcd7dd0be42e08d229c5209a9433105e474"),
    "C": (68_177_196, "9067a7d94acc472085a882b139e512a8e28dd1b39f761f9b9d2b4c6c0a697076"),
    "D": (397_844_161, "4f69add3177ebee40584aea2b4c1e8e40a81a1b64a32c455c27e502cd1e0b1f9"),
    "E": (64_156_624, "3afe043952a1d1c1b4d831365d3c12b4bdbd2c040d09e4c5a99defde964d76e5"),
    "F": (91_947_535, "8850c98d1453965fbce69fbb5cad17a71ee8714e5c1d2303dabff1c4e1e8ec55"),
    "G": (436_377_659, "19e9672274ef24bb194926ac0fc66e3c3c1b788ece74f1e478e9e04f7a7ddee3"),
    "H": (73_675_578, "d26731d109a6c0b74afe95e32f1471a290fa8e274f1e75c4ec7e55d70673375c"),
}
D0_BYTES = 111_087_381
D0_SHA256 = "d7c11d4258f7b78b5c3a20ebc5349e3a1f1f52a0279425701d14acbf79948842"
GTF_BYTES = 3_323_462_848
GTF_SHA256 = "ff32fd55c6799b3b94fe10aa17b2b5d4da952fa1de12fe44afadf32e949ec914"
BARCODES_BYTES = 115_512_960
BARCODES_SHA256 = "843a6f7038db8cb3c06f3dc21cc69d04139ffa10689780518d7a9e42dc2e819b"
LOCI = {
    "dense": "chr9:129925157-129927194",
    "sparse": "chr9:129861033-129861542",
    "absent": "chr9:129861034-129861542",
    "region": "chr9:84800000-85040000",
}

# CLI ASSUMPTIONS, isolated for review. No final PASS is possible unless the protocol records
# `interface_status=final`. An adapter JSON may replace these arrays without invoking a shell.
DEFAULT_ADAPTER: dict[str, Any] = {
    "schema": "gravlax.archive-root-cli-adapter.v1",
    "seal": ["{aie}", "seal-archive", "{input}", "--out", "{output}", "--json"],
    "inspect": ["{aie}", "inspect-archive", "{archive}", "--json"],
    "inspect_full": [
        "{aie}",
        "inspect-archive",
        "{archive}",
        "--json",
        "--verify-content",
    ],
    "collection_build_extra": ["--json"],
    "assumptions": [
        "seal-archive emits one strict JSON document",
        "inspect-archive emits one strict JSON document",
        "collection build --json reports identity_content_bytes_read, source_bytes_read, and source_sections_read",
        "v2 collection build uses the committed root when --source-digest is absent",
    ],
}


REQUIRED_FIXTURE_CASES = {
    "zero_sections": "accept",
    "one_section": "accept",
    "many_sections": "accept",
    "unknown_optional_section": "accept",
    "maximum_legal_name": "accept",
    "unselected_corruption_lazy_read": "accept",
    "bad_magic": "reject",
    "future_version": "reject",
    "root_mutation": "reject",
    "directory_byte_mutation": "reject",
    "entry_order_mutation": "reject",
    "name_mutation": "reject",
    "offset_mutation": "reject",
    "raw_length_mutation": "reject",
    "compressed_length_mutation": "reject",
    "payload_digest_mutation": "reject",
    "inline_header_mutation": "reject",
    "terminator_mutation": "reject",
    "footer_mutation": "reject",
    "payload_mutation": "reject",
    "swapped_payloads": "reject",
    "overlap": "reject",
    "gap": "reject",
    "truncation": "reject",
    "trailing_bytes": "reject",
    "decompressed_length_mismatch": "reject",
    "allocation_bomb": "reject",
    "compression_bomb": "reject",
    "selected_corruption_normal_read": "reject",
    "unselected_corruption_full_verify": "reject",
}


def parse_sample(value: str) -> tuple[str, Path]:
    try:
        sample, raw_path = value.split("=", 1)
    except ValueError as error:
        raise argparse.ArgumentTypeError("archive must be ID=PATH") from error
    if sample not in SAMPLES:
        raise argparse.ArgumentTypeError(f"archive ID must be one of {SAMPLES}")
    return sample, Path(raw_path)


def render(template: Sequence[str], values: Mapping[str, str]) -> list[str]:
    result = []
    for token in template:
        try:
            rendered = token.format_map(values)
        except KeyError as error:
            raise GateError(f"unknown CLI-adapter placeholder {error.args[0]!r}") from error
        require("{" not in rendered and "}" not in rendered, f"unresolved adapter token {rendered!r}")
        result.append(rendered)
    return result


def git_output(repo: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["git", "-C", str(repo), *arguments], capture_output=True, text=True, check=False
    )
    require(completed.returncode == 0, f"git {' '.join(arguments)} failed: {completed.stderr}")
    return completed.stdout.strip()


def command_version(command: Sequence[str]) -> str:
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    require(completed.returncode == 0, f"version command failed: {list(command)}")
    value = (completed.stdout or completed.stderr).strip()
    require(bool(value), f"version command was empty: {list(command)}")
    return value


def cargo_package_versions(lockfile: Path, names: set[str]) -> dict[str, list[str]]:
    with lockfile.open("rb") as handle:
        value = tomllib.load(handle)
    packages = value.get("package")
    require(type(packages) is list, "Cargo.lock package list is missing")
    result = {name: [] for name in names}
    for package in packages:
        if type(package) is dict and package.get("name") in names:
            version = package.get("version")
            require(type(version) is str, "Cargo.lock package version is malformed")
            result[package["name"]].append(version)
    require(all(result.values()), f"Cargo.lock lacks required packages: {[key for key, rows in result.items() if not rows]}")
    return {key: sorted(set(rows)) for key, rows in sorted(result.items())}


def host_memory_kib() -> int | None:
    path = Path("/proc/meminfo")
    if not path.is_file():
        return None
    match = re.search(r"^MemTotal:\s+(\d+)\s+kB$", path.read_text(), re.MULTILINE)
    return int(match.group(1)) if match else None


def is_promotable_run(
    interface_status: str,
    *,
    code_dirty: bool,
    scripts_dirty: bool,
    fixture_hooks: Sequence[Path],
) -> bool:
    """Return the single promotion predicate recorded in the immutable protocol.

    Keeping this predicate independently testable makes it impossible for the deliberately
    permissive dirty/development smoke interface to be mistaken for a Gate-A treatment run.
    The reducer adds a second, independent fail-closed check of the recorded state.
    """
    return (
        interface_status == "final"
        and not code_dirty
        and not scripts_dirty
        and bool(fixture_hooks)
    )


def planned_output_directories(run: Path) -> tuple[Path, ...]:
    """Enumerate every parent the driver, a child CLI, or a fixture hook writes beneath."""
    return (
        run,
        run / "archives",
        run / "archive-audits",
        run / "commands",
        run / "collections",
        run / "collections" / "build-benchmark",
        run / "replay",
        run / "replay" / "canonical",
        run / "replay" / "benchmark",
        run / "fixtures",
    )


def prepare_output_tree(run: Path) -> None:
    """Create the complete output-parent tree once, before the first measured command."""
    require(not run.exists(), f"refusing to overwrite {run}")
    for directory in planned_output_directories(run):
        directory.mkdir(parents=True, exist_ok=True)
    require(
        all(directory.is_dir() for directory in planned_output_directories(run)),
        "failed to create the complete output-parent tree",
    )


def query_command(aie: Path, collection: Path, kind: str) -> list[str]:
    if kind in ("dense", "sparse", "absent"):
        return [
            str(aie), "collection", "junction", str(collection), LOCI[kind],
            "--json", "--explain", "--top", "0",
        ]
    if kind == "region":
        return [
            str(aie), "collection", "region", str(collection), LOCI["region"],
            "--json", "--explain", "--top", "0",
        ]
    if kind == "jset":
        return [
            str(aie), "collection", "jset", str(collection),
            "--include", LOCI["dense"], "--exclude", LOCI["sparse"],
            "--json", "--explain", "--top", "0",
        ]
    raise GateError(f"unknown query kind {kind}")


def replay_command(
    aie: Path,
    archive: Path,
    gtf: Path,
    barcodes: Path,
    output: Path,
    mode: str,
) -> list[str]:
    require(mode in {"gene-stream", "gene-eager", "velocity-stream", "velocity-eager"}, "bad replay mode")
    command = [
        str(aie), "replay-rows", str(archive), "--gtf", str(gtf),
        "--barcodes", str(barcodes), "--out-dir", str(output),
    ]
    if mode.endswith("eager"):
        command.append("--eager")
    if mode.startswith("velocity"):
        command.append("--velocity")
    return command


def collection_build_command(
    aie: Path,
    adapter: Mapping[str, Any],
    samples: Sequence[str],
    archives: Mapping[str, Path],
    output: Path,
    arm: str,
    base: Path | None = None,
) -> list[str]:
    require(arm in {"v1", "v2"}, "invalid collection arm")
    command = [str(aie), "collection", "build"]
    if base is not None:
        command.extend(["--base", str(base)])
    for sample in samples:
        command.extend(["--sample", f"{sample}={archives[sample]}"])
        if arm == "v1":
            command.extend(["--source-digest", f"{sample}={SOURCE[sample][1]}"])
    command.extend(["--out", str(output)])
    extra = adapter.get("collection_build_extra")
    require(type(extra) is list and all(type(token) is str for token in extra), "bad collection adapter")
    command.extend(extra)
    return command


def validate_fixture_result(path: Path, expected_binary_sha256: str) -> dict[str, str]:
    value = load_json(path)
    require(type(value) is dict, f"{path}: fixture result is not an object")
    require(value.get("schema") == "gravlax.archive-root-adversarial.v1", f"{path}: fixture schema differs")
    require(value.get("binary_sha256") == expected_binary_sha256, f"{path}: fixture binary differs")
    cases = value.get("cases")
    require(type(cases) is list and cases, f"{path}: fixture cases are empty")
    seen: dict[str, str] = {}
    for case in cases:
        require(type(case) is dict, f"{path}: malformed fixture case")
        require(
            set(case)
            == {
                "id",
                "expect",
                "command",
                "exit_status",
                "stdout",
                "stderr",
                "fixture_bytes",
                "fixture_sha256",
                "mutation",
                "rejected_before_scientific_output",
            },
            f"{path}: fixture case fields differ",
        )
        case_id = case["id"]
        require(type(case_id) is str and case_id not in seen, f"{path}: duplicate fixture case")
        expect = case["expect"]
        require(expect in {"accept", "reject"}, f"{path}: invalid fixture expectation")
        seen[case_id] = expect
        require(type(case["command"]) is list and all(type(v) is str for v in case["command"]), "bad fixture command")
        require(type(case["exit_status"]) is int, f"{path}: fixture exit status is not an integer")
        require(type(case["fixture_bytes"]) is int and case["fixture_bytes"] >= 0, f"{path}: invalid fixture size")
        require(re.fullmatch(r"[0-9a-f]{64}", str(case["fixture_sha256"])) is not None, f"{path}: invalid fixture SHA-256")
        require(type(case["mutation"]) is str and case["mutation"], f"{path}: missing fixture mutation description")
        for stream in ("stdout", "stderr"):
            relative = case[stream]
            require(type(relative) is str and not Path(relative).is_absolute(), f"{path}: unsafe fixture path")
            stream_path = path.parent / relative
            require(stream_path.is_file() and stream_path.resolve().is_relative_to(path.parent.resolve()), "missing fixture stream")
        if expect == "reject":
            require(case["exit_status"] != 0, f"{path}: corruption was accepted")
            require(case["rejected_before_scientific_output"] is True, f"{path}: partial output was emitted")
            require((path.parent / case["stdout"]).stat().st_size == 0, f"{path}: rejected case wrote stdout")
        else:
            require(case["exit_status"] == 0, f"{path}: valid fixture was rejected")
            require(case["rejected_before_scientific_output"] is False, f"{path}: accepted case marked rejected")
            load_json(path.parent / case["stdout"])
    return seen


def load_collection_build_stdout(path: Path, interface_status: str) -> dict[str, Any] | None:
    """Accept legacy text only for an explicitly non-promotable development run."""
    try:
        value = load_json(path)
    except GateError:
        require(
            interface_status == "draft",
            f"{path}: final interface requires structured collection-build JSON",
        )
        return None
    require(type(value) is dict, f"{path}: collection-build JSON is not an object")
    return value


def validate_seal_stdout(
    value: Any,
    *,
    source: Path,
    sealed: Path,
    source_blake3: str | None,
    audit: Mapping[str, Any],
) -> None:
    require(type(value) is dict, "seal output is not an object")
    require(
        set(value)
        == {
            "schema",
            "input",
            "output",
            "source_format_version",
            "output_format_version",
            "sections",
            "compressed_payload_bytes_copied",
            "input_bytes",
            "output_bytes",
            "source_full_file_blake3",
            "source_identity_content_bytes_read",
            "archive_root",
            "encoded_sections_identity",
            "payload_verification_bytes_read",
            "elapsed_seconds",
        },
        "seal output fields changed",
    )
    require(value["schema"] == "gravlax.archive.seal.v1", "seal schema changed")
    require(Path(value["input"]).resolve() == source, "seal input path differs")
    require(Path(value["output"]).resolve() == sealed, "seal output path differs")
    require(value["source_format_version"] == 1 and value["output_format_version"] == 2, "seal versions differ")
    require(value["sections"] == audit["section_count"], "seal section count differs")
    copied = sum(row["compressed_bytes"] for row in audit["sections"])
    require(value["compressed_payload_bytes_copied"] == copied, "seal copied-byte count differs")
    require(value["payload_verification_bytes_read"] == copied, "seal verification-byte count differs")
    require(value["input_bytes"] == audit["source_bytes"], "seal input-byte count differs")
    require(value["output_bytes"] == audit["sealed_bytes"], "seal output-byte count differs")
    require(
        value["source_identity_content_bytes_read"] == audit["source_bytes"],
        "legacy identity scan did not account for the complete source file",
    )
    require(
        type(value["source_full_file_blake3"]) is str
        and re.fullmatch(r"[0-9a-f]{64}", value["source_full_file_blake3"]),
        "seal source BLAKE3 is malformed",
    )
    if source_blake3 is not None:
        require(value["source_full_file_blake3"] == source_blake3, "source BLAKE3 differs from frozen input")
    require(
        value["archive_root"]
        == {"scheme": "aie-directory-root-v2", "blake3": audit["archive_root"]},
        "seal root differs from independent audit",
    )
    encoded = value["encoded_sections_identity"]
    require(
        type(encoded) is dict
        and encoded.get("scheme") == "aie-encoded-sections-v1"
        and re.fullmatch(r"[0-9a-f]{64}", str(encoded.get("blake3"))) is not None,
        "seal encoded-section identity is malformed",
    )
    require(type(value["elapsed_seconds"]) in (int, float) and value["elapsed_seconds"] >= 0, "bad seal elapsed time")


def validate_inspect_stdout(
    value: Any,
    *,
    archive: Path,
    audit: Mapping[str, Any],
    full: bool,
) -> None:
    require(type(value) is dict, "inspect output is not an object")
    require(
        set(value)
        == {
            "schema",
            "archive",
            "format_version",
            "file_bytes",
            "sections",
            "native_identity",
            "encoded_sections_identity",
            "verification",
            "elapsed_seconds",
        },
        "inspect output fields changed",
    )
    require(value["schema"] == "gravlax.archive.identity.v1", "inspect schema changed")
    require(Path(value["archive"]).resolve() == archive, "inspect archive path differs")
    require(value["format_version"] == 2, "inspect did not report v2")
    require(value["file_bytes"] == audit["sealed_bytes"] and value["sections"] == audit["section_count"], "inspect size/count differs")
    require(
        value["native_identity"]
        == {"scheme": "aie-directory-root-v2", "blake3": audit["archive_root"]},
        "inspect native root differs",
    )
    encoded = value["encoded_sections_identity"]
    require(type(encoded) is dict and encoded.get("scheme") == "aie-encoded-sections-v1", "inspect encoded identity scheme differs")
    verification = value["verification"]
    require(
        type(verification) is dict
        and set(verification)
        == {
            "directory_and_root",
            "all_payloads",
            "identity_content_bytes_read",
            "ordinary_reads_verify_selected_payloads_only",
        },
        "inspect verification fields changed",
    )
    payload_bytes = sum(row["compressed_bytes"] for row in audit["sections"])
    require(verification["directory_and_root"] is True, "inspect did not validate directory/root")
    require(verification["all_payloads"] is full, "inspect payload-verification flag differs")
    require(
        verification["identity_content_bytes_read"] == (payload_bytes if full else 0),
        "inspect payload-read accounting differs",
    )
    require(verification["ordinary_reads_verify_selected_payloads_only"] is True, "lazy verification contract differs")


_PLAN_FIELDS = {
    "archive",
    "actual_archive_bytes_read",
    "planned_compressed_bytes",
    "chunks_decoded",
    "unique_chunks_decoded",
    "independent_chunk_decodes",
}


def scientific_query_projection(value: Any) -> Any:
    """Remove only source-location, timing, and execution-plan fields from a query result."""
    if isinstance(value, list):
        return [scientific_query_projection(item) for item in value]
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            if key in {"planning", "explain"} or key in _PLAN_FIELDS or key.endswith("_seconds"):
                continue
            result[key] = scientific_query_projection(item)
        return result
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--aie", type=Path, required=True)
    parser.add_argument("--code-repo", type=Path)
    parser.add_argument("--scripts-repo", type=Path)
    parser.add_argument("--archive", action="append", type=parse_sample, default=[])
    parser.add_argument("--d0-archive", type=Path)
    parser.add_argument("--gtf", type=Path)
    parser.add_argument("--barcodes", type=Path)
    parser.add_argument("--time-binary", type=Path, default=Path("/usr/bin/time"))
    parser.add_argument("--scratch-dir", type=Path)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--warm-blocks", type=int, default=8)
    parser.add_argument("--cli-adapter", type=Path)
    parser.add_argument("--interface-status", choices=("draft", "final"), default="draft")
    parser.add_argument(
        "--fixture-hook",
        action="append",
        type=Path,
        default=[],
        help=(
            "Executable accepting --aie, --out-dir, --root-domain-hex, --fixture-archive, "
            "and --source-archive, and writing result.json"
        ),
    )
    parser.add_argument("--allow-dirty-code", action="store_true")
    parser.add_argument("--allow-dirty-scripts", action="store_true")
    args = parser.parse_args()

    self_test_blake3()
    require(args.threads == 24, "the frozen gate requires exactly 24 threads")
    require(args.warm_blocks == 8, "the frozen gate requires exactly eight paired warm blocks")
    root = args.project_root.resolve()
    run = args.run_dir.resolve()
    aie = args.aie.resolve()
    code_repo = (args.code_repo or root / "work/gravlax").resolve()
    scripts_repo = (args.scripts_repo or root / "gravlax-paper-scripts").resolve()
    time_binary = args.time_binary.resolve()
    require(aie.is_file() and os.access(aie, os.X_OK), f"missing executable {aie}")
    require(time_binary.is_file() and os.access(time_binary, os.X_OK), f"missing time binary {time_binary}")
    require((code_repo / ".git").exists(), f"not a Git checkout: {code_repo}")
    require((scripts_repo / ".git").exists(), f"not a Git checkout: {scripts_repo}")
    require(not run.exists(), f"refusing to overwrite {run}")

    supplied = dict(args.archive)
    require(len(supplied) == len(args.archive), "duplicate --archive ID")
    source_archives = {
        sample: (
            supplied.get(sample)
            or root / "runs/post-v1/sez-transcript-end-atlas-r1" / f"donor-{sample}" / "sez.called.aie"
        ).resolve()
        for sample in SAMPLES
    }
    require(set(supplied).issubset(SAMPLES), "unknown archive IDs")
    for sample, path in source_archives.items():
        require(path.is_file(), f"missing source archive {sample}: {path}")
        require(path.stat().st_size == SOURCE[sample][0], f"source archive {sample} byte size differs")
        parsed = parse_archive(path, ROOT_DOMAIN, verify_payloads=False)
        require(parsed.version == 1, f"source archive {sample} is not v1")

    d0 = (args.d0_archive or root / "runs/archive/d0.aie").resolve()
    gtf = (args.gtf or root / "annotations/gencode.v49.annotation.gtf").resolve()
    barcodes = (
        args.barcodes or root / "runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv"
    ).resolve()
    for path in (d0, gtf, barcodes):
        require(path.is_file(), f"missing replay input {path}")
    require(d0.stat().st_size == D0_BYTES and sha256_file(d0) == D0_SHA256, "D0 archive identity differs")
    require(gtf.stat().st_size == GTF_BYTES and sha256_file(gtf) == GTF_SHA256, "v49 GTF identity differs")
    require(
        barcodes.stat().st_size == BARCODES_BYTES and sha256_file(barcodes) == BARCODES_SHA256,
        "D0 barcode dictionary identity differs",
    )
    require(parse_archive(d0, ROOT_DOMAIN, verify_payloads=False).version == 1, "D0 source is not v1")

    adapter = DEFAULT_ADAPTER if args.cli_adapter is None else load_json(args.cli_adapter.resolve())
    require(type(adapter) is dict and adapter.get("schema") == DEFAULT_ADAPTER["schema"], "CLI adapter schema differs")
    for key in ("seal", "inspect", "inspect_full", "collection_build_extra", "assumptions"):
        require(type(adapter.get(key)) is list and all(type(v) is str for v in adapter[key]), f"bad adapter {key}")

    commit = git_output(code_repo, "rev-parse", "HEAD")
    require(re.fullmatch(r"[0-9a-f]{40}", commit) is not None, "invalid Gravlax commit")
    dirty = git_output(code_repo, "status", "--porcelain")
    if dirty:
        require(args.allow_dirty_code, "Gravlax worktree is dirty")
    scripts_commit = git_output(scripts_repo, "rev-parse", "HEAD")
    scripts_dirty = git_output(scripts_repo, "status", "--porcelain")
    if scripts_dirty:
        require(args.allow_dirty_scripts, "reproducibility worktree is dirty")
    promotable = is_promotable_run(
        args.interface_status,
        code_dirty=bool(dirty),
        scripts_dirty=bool(scripts_dirty),
        fixture_hooks=args.fixture_hook,
    )

    prepare_output_tree(run)
    environment = {"RAYON_NUM_THREADS": str(args.threads)}
    binary_sha = sha256_file(aie)
    protocol = {
        "schema": "gravlax.archive-root-gate-protocol.v1",
        "date": DATE,
        "interface_status": args.interface_status,
        "promotable_run": promotable,
        "root_domain_ascii": "gravlax-aie-directory-root-v2\\0",
        "root_domain_hex": ROOT_DOMAIN.hex(),
        "root_preimage": "domain || AIE0 || u32_le(version=2) || u64_le(directory_offset) || exact_directory_bytes",
        "threads": args.threads,
        "warm_paired_blocks": args.warm_blocks,
        "schedule": "odd blocks v1-v2; even blocks v2-v1",
        "project_root": str(root),
        "run_dir": str(run),
        "aie": {"path": str(aie), "bytes": aie.stat().st_size, "sha256": binary_sha},
        "gravlax": {"repo": str(code_repo), "commit": commit, "dirty": bool(dirty)},
        "paper_scripts": {
            "repo": str(scripts_repo),
            "commit": scripts_commit,
            "dirty": bool(scripts_dirty),
        },
        "adapter": adapter,
        "inputs": {
            "archives": [
                {"id": sample, "path": str(source_archives[sample]), "bytes": SOURCE[sample][0], "expected_blake3": SOURCE[sample][1]}
                for sample in SAMPLES
            ],
            "d0": {"path": str(d0), "bytes": D0_BYTES, "sha256": D0_SHA256},
            "gtf": {"path": str(gtf), "bytes": GTF_BYTES, "sha256": GTF_SHA256},
            "barcodes": {"path": str(barcodes), "bytes": BARCODES_BYTES, "sha256": BARCODES_SHA256},
        },
        "fixture_hooks": [str(path.resolve()) for path in args.fixture_hook],
        "host": {
            "hostname": socket.gethostname(),
            "platform": platform.platform(),
            "python": platform.python_version(),
            "cpu_count": os.cpu_count(),
            "processor": platform.processor(),
            "memory_total_kib": host_memory_kib(),
        },
        "tools": {
            "aie": command_version([str(aie), "--version"]),
            "rustc": command_version(["rustc", "--version"]),
            "cargo": command_version(["cargo", "--version"]),
            "gnu_time": command_version([str(time_binary), "--version"]).splitlines()[0],
            "cargo_lock_packages": cargo_package_versions(
                code_repo / "Cargo.lock", {"blake3", "zstd", "zstd-safe", "zstd-sys"}
            ),
        },
    }
    write_json(run / "protocol.json", protocol)

    all_v1 = {**source_archives, "D0": d0}
    all_v2: dict[str, Path] = {}
    scratch_parent = args.scratch_dir.resolve() if args.scratch_dir else None
    if scratch_parent is not None:
        scratch_parent.mkdir(parents=True, exist_ok=True)
    for sample, source in all_v1.items():
        sealed = run / "archives" / f"{sample}.v2.aie"
        all_v2[sample] = sealed
        command = render(adapter["seal"], {"aie": str(aie), "input": str(source), "output": str(sealed)})
        run_timed(run / "commands/seal" / sample, command, time_binary=time_binary, cwd=root, environment=environment)
        seal_stdout = load_json(run / "commands/seal" / sample / "stdout.txt")
        audit = compare_v1_v2(source, sealed, ROOT_DOMAIN)
        validate_seal_stdout(
            seal_stdout,
            source=source,
            sealed=sealed,
            source_blake3=SOURCE[sample][1] if sample in SOURCE else None,
            audit=audit,
        )
        with tempfile.TemporaryDirectory(prefix=f"gravlax-root-{sample}-", dir=scratch_parent) as temporary:
            repeated = Path(temporary) / f"{sample}.repeat.v2.aie"
            repeat_command = render(adapter["seal"], {"aie": str(aie), "input": str(source), "output": str(repeated)})
            run_timed(
                run / "commands/seal-repeat" / sample,
                repeat_command,
                time_binary=time_binary,
                cwd=root,
                environment=environment,
            )
            audit["repeat_sha256"] = sha256_file(repeated)
        audit["sealed_sha256"] = sha256_file(sealed)
        require(
            audit["repeat_sha256"] == audit["sealed_sha256"],
            f"{sample}: repeated seal is not byte-identical by SHA-256",
        )
        write_json(run / "archive-audits" / f"{sample}.json", audit)

        for mode, template in (("lazy", adapter["inspect"]), ("full", adapter["inspect_full"])):
            inspect_command = render(template, {"aie": str(aie), "archive": str(sealed)})
            command_dir = run / "commands/inspect" / sample / mode
            run_timed(command_dir, inspect_command, time_binary=time_binary, cwd=root, environment=environment)
            inspect_stdout = load_json(command_dir / "stdout.txt")
            validate_inspect_stdout(
                inspect_stdout,
                archive=sealed,
                audit=audit,
                full=mode == "full",
            )
            require(
                inspect_stdout["encoded_sections_identity"]
                == seal_stdout["encoded_sections_identity"],
                f"{sample}: seal and inspect encoded-section identities differ",
            )

    canonical_collections: dict[tuple[str, str], Path] = {}
    for arm, archives in (("v1", source_archives), ("v2", {sample: all_v2[sample] for sample in SAMPLES})):
        for form, samples, base_form in (
            ("root", SAMPLES, None),
            ("base", SAMPLES[:4], None),
            ("extension", SAMPLES[4:], "base"),
        ):
            output = run / "collections" / f"{arm}-{form}.aicollection"
            base = canonical_collections.get((arm, base_form)) if base_form else None
            command = collection_build_command(aie, adapter, samples, archives, output, arm, base)
            run_timed(
                run / "commands/collection-canonical" / f"{arm}-{form}",
                command,
                time_binary=time_binary,
                cwd=root,
                environment=environment,
            )
            load_collection_build_stdout(
                run / "commands/collection-canonical" / f"{arm}-{form}" / "stdout.txt",
                args.interface_status,
            )
            canonical_collections[(arm, form)] = output

        reverse = run / "collections" / f"{arm}-root-reverse.aicollection"
        reverse_command = collection_build_command(
            aie, adapter, tuple(reversed(SAMPLES)), archives, reverse, arm
        )
        reverse_dir = run / "commands/collection-canonical" / f"{arm}-root-reverse"
        run_timed(
            reverse_dir,
            reverse_command,
            time_binary=time_binary,
            cwd=root,
            environment=environment,
        )
        load_collection_build_stdout(reverse_dir / "stdout.txt", args.interface_status)
        require(
            sha256_file(reverse) == sha256_file(canonical_collections[(arm, "root")]),
            f"{arm}: reversed-input root build is not byte-identical",
        )

    canonical_query_values: dict[tuple[str, str, str], Any] = {}
    for form in ("root", "extension"):
        for kind in ("dense", "sparse", "absent", "region", "jset"):
            for arm in ("v1", "v2"):
                command_dir = run / "commands/query-canonical" / f"{arm}-{form}-{kind}"
                run_timed(
                    command_dir,
                    query_command(aie, canonical_collections[(arm, form)], kind),
                    time_binary=time_binary,
                    cwd=root,
                    environment=environment,
                )
                query_value = load_json(command_dir / "stdout.txt")
                canonical_query_values[(arm, form, kind)] = query_value

    query_equivalence: dict[str, Any] = {}
    for kind in ("dense", "sparse", "absent", "region", "jset"):
        projections = {
            f"{arm}-{form}": scientific_query_projection(
                canonical_query_values[(arm, form, kind)]
            )
            for arm in ("v1", "v2")
            for form in ("root", "extension")
        }
        reference = projections["v1-root"]
        require(
            all(value == reference for value in projections.values()),
            f"{kind}: v1/v2 or fresh-root/extension scientific results differ",
        )
        query_equivalence[kind] = reference
    write_json(
        run / "archive-audits/query-equivalence.json",
        {
            "schema": "gravlax.archive-root-query-equivalence.v1",
            "arms": ["v1-root", "v1-extension", "v2-root", "v2-extension"],
            "queries": query_equivalence,
        },
    )

    build_schedule = ["block\tposition\tarm\tlabel"]
    for block in range(1, args.warm_blocks + 1):
        order = ("v1", "v2") if block % 2 else ("v2", "v1")
        for position, arm in enumerate(order, 1):
            label = f"b{block:02d}-p{position}-{arm}"
            output = run / "collections/build-benchmark" / f"{label}.aicollection"
            archives = source_archives if arm == "v1" else {sample: all_v2[sample] for sample in SAMPLES}
            command = collection_build_command(aie, adapter, SAMPLES, archives, output, arm)
            run_timed(
                run / "commands/collection-build-benchmark" / label,
                command,
                time_binary=time_binary,
                cwd=root,
                environment=environment,
            )
            load_collection_build_stdout(
                run / "commands/collection-build-benchmark" / label / "stdout.txt",
                args.interface_status,
            )
            build_schedule.append(f"{block}\t{position}\t{arm}\t{label}")
    (run / "collection-build-schedule.tsv").write_text("\n".join(build_schedule) + "\n")

    replay_trees: dict[tuple[str, str], dict[str, Any]] = {}
    for arm, archive in (("v1", d0), ("v2", all_v2["D0"])):
        for mode in ("gene-stream", "gene-eager", "velocity-stream", "velocity-eager"):
            output = run / "replay/canonical" / f"{arm}-{mode}"
            command_dir = run / "commands/replay-canonical" / f"{arm}-{mode}"
            run_timed(
                command_dir,
                replay_command(aie, archive, gtf, barcodes, output, mode),
                time_binary=time_binary,
                cwd=root,
                environment=environment,
            )
            replay_trees[(arm, mode)] = tree_digest(output)
            write_json(run / "archive-audits" / f"replay-{arm}-{mode}.json", replay_trees[(arm, mode)])
    for mode in ("gene-stream", "gene-eager", "velocity-stream", "velocity-eager"):
        require(replay_trees[("v1", mode)] == replay_trees[("v2", mode)], f"{mode}: v1/v2 replay differs")
    for arm in ("v1", "v2"):
        require(replay_trees[(arm, "gene-stream")] == replay_trees[(arm, "gene-eager")], f"{arm}: Gene eager differs")
        require(replay_trees[(arm, "velocity-stream")] == replay_trees[(arm, "velocity-eager")], f"{arm}: Velocity eager differs")

    replay_schedule = ["block\tposition\tarm\tlabel"]
    for block in range(1, args.warm_blocks + 1):
        order = ("v1", "v2") if block % 2 else ("v2", "v1")
        for position, arm in enumerate(order, 1):
            label = f"b{block:02d}-p{position}-{arm}"
            archive = d0 if arm == "v1" else all_v2["D0"]
            output = run / "replay/benchmark" / label
            command_dir = run / "commands/replay-benchmark" / label
            run_timed(
                command_dir,
                replay_command(aie, archive, gtf, barcodes, output, "gene-stream"),
                time_binary=time_binary,
                cwd=root,
                environment=environment,
            )
            require(tree_digest(output) == replay_trees[("v1", "gene-stream")], f"{label}: timed replay differs")
            replay_schedule.append(f"{block}\t{position}\t{arm}\t{label}")
    (run / "replay-schedule.tsv").write_text("\n".join(replay_schedule) + "\n")

    query_schedule = ["kind\tblock\tposition\tarm\tlabel"]
    for kind in ("dense", "sparse", "jset"):
        for block in range(1, args.warm_blocks + 1):
            order = ("v1", "v2") if block % 2 else ("v2", "v1")
            for position, arm in enumerate(order, 1):
                label = f"{kind}-b{block:02d}-p{position}-{arm}"
                command_dir = run / "commands/query-benchmark" / label
                run_timed(
                    command_dir,
                    query_command(aie, canonical_collections[(arm, "root")], kind),
                    time_binary=time_binary,
                    cwd=root,
                    environment=environment,
                )
                value = load_json(command_dir / "stdout.txt")
                require(
                    scientific_query_projection(value)
                    == scientific_query_projection(canonical_query_values[("v1", "root", kind)]),
                    f"{label}: timed query result differs from the canonical result",
                )
                query_schedule.append(f"{kind}\t{block}\t{position}\t{arm}\t{label}")
    (run / "query-schedule.tsv").write_text("\n".join(query_schedule) + "\n")

    observed_fixture_cases: dict[str, str] = {}
    hook_names: set[str] = set()
    for hook in args.fixture_hook:
        hook = hook.resolve()
        require(hook.is_file() and os.access(hook, os.X_OK), f"fixture hook is not executable: {hook}")
        name = re.sub(r"[^A-Za-z0-9_.-]+", "-", hook.stem)
        require(name and name not in hook_names, "duplicate fixture hook name")
        hook_names.add(name)
        output = run / "fixtures" / name
        command = [
            str(hook), "--aie", str(aie), "--out-dir", str(output),
            "--root-domain-hex", ROOT_DOMAIN.hex(),
            "--fixture-archive", str(all_v2["A"]),
            "--source-archive", str(source_archives["A"]),
        ]
        run_timed(
            run / "commands/fixture-hooks" / name,
            command,
            time_binary=time_binary,
            cwd=root,
            environment=environment,
        )
        hook_cases = validate_fixture_result(output / "result.json", binary_sha)
        overlap = set(observed_fixture_cases) & set(hook_cases)
        require(not overlap, f"fixture case IDs repeated across hooks: {sorted(overlap)}")
        observed_fixture_cases.update(hook_cases)
    if args.interface_status == "final":
        missing = set(REQUIRED_FIXTURE_CASES) - set(observed_fixture_cases)
        require(not missing, f"fixture coverage missing {sorted(missing)}")
        mismatched = {
            case: (REQUIRED_FIXTURE_CASES[case], observed_fixture_cases[case])
            for case in REQUIRED_FIXTURE_CASES
            if observed_fixture_cases.get(case) != REQUIRED_FIXTURE_CASES[case]
        }
        require(not mismatched, f"fixture expectations differ: {mismatched}")

    write_json(
        run / "driver-completion.json",
        {
            "schema": "gravlax.archive-root-driver-completion.v1",
            "status": "COMPLETE" if promotable else "DEVELOPMENT_ONLY",
            "promotable_run": promotable,
            "fixture_cases": [
                {"id": case, "expect": observed_fixture_cases[case]}
                for case in sorted(observed_fixture_cases)
            ],
            "root_domain_hex": ROOT_DOMAIN.hex(),
        },
    )
    manifest = write_artifact_manifest(run)
    print(json.dumps({"run_dir": str(run), "artifact_manifest": manifest}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except GateError as error:
        raise SystemExit(f"archive-root benchmark failed: {error}") from error
