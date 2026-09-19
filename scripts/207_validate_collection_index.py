#!/usr/bin/env python3
"""Fail-closed reduction of the locked eight-archive collection-index benchmark."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
import shlex
from pathlib import Path
from typing import Any


EXPECTED_RUN = "collection-index-v2-r5"
EXPECTED_COMMIT = "6b92503cffbd48da8312bdd75c1f1f70d5ff2a27"
EXPECTED_BINARY_SHA256 = "3c5faaf38d8b9ea14001a8561b78b026452804d454b16a5f2e11dbc089605a4d"
SAMPLES = tuple("ABCDEFGH")
SOURCE = {
    "A": (110_183_760, "029de28423c59a8b9ca283e73687a9501cd7533cfd508788b61b15e79b751463", 853),
    "B": (83_674_178, "4feced152203593f22ddb5a13207dbcd7dd0be42e08d229c5209a9433105e474", 855),
    "C": (68_177_196, "9067a7d94acc472085a882b139e512a8e28dd1b39f761f9b9d2b4c6c0a697076", 851),
    "D": (397_844_161, "4f69add3177ebee40584aea2b4c1e8e40a81a1b64a32c455c27e502cd1e0b1f9", 871),
    "E": (64_156_624, "3afe043952a1d1c1b4d831365d3c12b4bdbd2c040d09e4c5a99defde964d76e5", 841),
    "F": (91_947_535, "8850c98d1453965fbce69fbb5cad17a71ee8714e5c1d2303dabff1c4e1e8ec55", 849),
    "G": (436_377_659, "19e9672274ef24bb194926ac0fc66e3c3c1b788ece74f1e478e9e04f7a7ddee3", 874),
    "H": (73_675_578, "d26731d109a6c0b74afe95e32f1471a290fa8e274f1e75c4ec7e55d70673375c", 854),
}
REFERENCE = {
    "algo": "aie-genome-blake3-v1",
    "digest": "2817ea0bc4a2b919ef6007efdb5504757dd9cb76bf550228e510112a41cb7472",
    "stamped": True,
}
LOCI = {
    "dense": "chr9:129925157-129927194",
    "sparse": "chr9:129861033-129861542",
    "absent": "chr9:129861034-129861542",
    "region": "chr9:84800000-85040000",
}
EXPECTED_TOTALS = {
    "dense": {"cells": 180, "umis": 182},
    "sparse": {"cells": 4, "umis": 4},
    "absent": {"cells": 0, "umis": 0},
    "region": {"cells": 20_019, "molecules": 74_848, "umis": 74_352},
    "jset": {
        "both": 0,
        "exclude_only": 4,
        "include_only": 182,
        "informative_umis": 186,
        "usage_fraction": 182 / 186,
    },
}
EXPECTED_INDEX = {
    "root": {
        "archive_routes": 1_941_890,
        "chunk_postings": 1_945_703,
        "global_junctions": None,
        "interval_chunks": 6_848,
        "junction_bin_bp": 16_000_000,
        "junction_segments": 282,
        "segment_junction_rows": 1_344_270,
    },
    "chain": {
        "archive_routes": 1_941_890,
        "chunk_postings": 1_945_703,
        "global_junctions": None,
        "interval_chunks": 6_848,
        "junction_bin_bp": 16_000_000,
        "junction_segments": 552,
        "segment_junction_rows": 1_521_092,
    },
}
COLLECTION_FILES = (
    "sez4-base.aicollection",
    "sez8-extension.aicollection",
    "sez8-reverse.aicollection",
    "sez8-root.aicollection",
)
QUERY_KINDS = ("dense", "sparse", "absent", "region", "jset")


class ValidationError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def require_keys(value: dict[str, Any], expected: set[str], label: str) -> None:
    require(set(value) == expected, f"{label}: fields changed: {sorted(value)}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> Any:
    try:
        raw = path.read_bytes()
        require(raw.endswith(b"\n"), f"{path}: missing final newline")
        return json.loads(raw)
    except (OSError, json.JSONDecodeError, UnicodeDecodeError) as error:
        raise ValidationError(f"cannot read JSON {path}: {error}") from error


def parse_wall(value: str) -> float:
    fields = [float(field) for field in value.split(":")]
    require(1 <= len(fields) <= 3, f"invalid GNU-time wall value: {value!r}")
    seconds = fields[-1]
    if len(fields) >= 2:
        seconds += 60 * fields[-2]
    if len(fields) == 3:
        seconds += 3600 * fields[-3]
    return seconds


def read_time(path: Path) -> dict[str, Any]:
    text = path.read_text()
    patterns = {
        "command": r'^\s*Command being timed: "(.*)"$',
        "wall": r"^\s*Elapsed \(wall clock\) time.*: ([0-9:.]+)$",
        "rss": r"^\s*Maximum resident set size \(kbytes\): (\d+)$",
        "status": r"^\s*Exit status: (\d+)$",
    }
    found = {
        key: re.findall(pattern, text, flags=re.MULTILINE)
        for key, pattern in patterns.items()
    }
    for key, values in found.items():
        require(len(values) == 1, f"{path}: expected one {key} record, got {len(values)}")
    record = {
        "command": found["command"][0],
        "wall_seconds": parse_wall(found["wall"][0]),
        "max_rss_kib": int(found["rss"][0]),
        "exit_status": int(found["status"][0]),
        "record_sha256": sha256(path),
    }
    require(record["exit_status"] == 0, f"{path}: nonzero exit status")
    require(record["max_rss_kib"] > 0, f"{path}: maximum RSS is not positive")
    return record


def expected_files() -> set[str]:
    names = {
        "artifacts.sha256",
        "base-after.sha256",
        "base-before.sha256",
        *COLLECTION_FILES,
    }
    for stem in ("root", "reverse", "base", "extension"):
        names.add(f"{stem}-build.stdout.txt")
        names.add(f"{stem}-build.time.txt")
    for label in ("root-inspect", "chain-inspect", "root-inspect-verify"):
        names.add(f"{label}.json")
        names.add(f"{label}.time.txt")
    for arm in ("root", "chain"):
        for kind in QUERY_KINDS:
            names.add(f"{arm}-{kind}.json")
            names.add(f"{arm}-{kind}.time.txt")
    for label in ("root-support-pruned", "root-verify-dense"):
        names.add(f"{label}.json")
        names.add(f"{label}.time.txt")
    for kind in QUERY_KINDS:
        names.add(f"naive-{kind}.time.txt")
    for sample in SAMPLES:
        for kind in QUERY_KINDS:
            names.add(f"naive/{sample}-{kind}.json")
        names.add(f"naive/{sample}-sparse.stderr.txt")
        names.add(f"naive/{sample}-absent.stderr.txt")
        names.add(f"naive/{sample}-absent.stdout.txt")
    return names


def validate_file_set(run: Path) -> dict[str, Any]:
    actual = {
        str(path.relative_to(run))
        for path in run.rglob("*")
        if path.is_file()
    }
    expected = expected_files()
    require(actual == expected, f"run file set changed; missing={sorted(expected-actual)}, extra={sorted(actual-expected)}")
    entries = []
    for relative in sorted(actual):
        path = run / relative
        entries.append((relative, path.stat().st_size, sha256(path)))
    canonical = b"".join(
        relative.encode() + b"\0" + str(size).encode() + b"\0" + digest.encode() + b"\n"
        for relative, size, digest in entries
    )
    return {
        "files": len(entries),
        "bytes": sum(size for _, size, _ in entries),
        "manifest_sha256": hashlib.sha256(canonical).hexdigest(),
        "manifest_encoding": "relative_path NUL decimal_bytes NUL sha256 LF, sorted by relative_path",
    }


def parse_sha_manifest(path: Path) -> dict[str, tuple[str, Path]]:
    result: dict[str, tuple[str, Path]] = {}
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        require(match is not None, f"{path}:{line_number}: malformed SHA-256 line")
        digest, raw_path = match.groups()
        name = Path(raw_path).name
        require(name not in result, f"{path}: duplicate basename {name}")
        result[name] = (digest, Path(raw_path))
    return result


def validate_checksums(run: Path, binary: Path) -> dict[str, Any]:
    manifest = parse_sha_manifest(run / "artifacts.sha256")
    expected_names = {binary.name, *COLLECTION_FILES}
    require(set(manifest) == expected_names, "artifacts.sha256 entries changed")
    binary_digest, binary_recorded = manifest[binary.name]
    require(binary.resolve() == binary_recorded.resolve(), "binary path differs from artifacts.sha256")
    require(binary_digest == EXPECTED_BINARY_SHA256, "unexpected final binary digest")
    require(sha256(binary) == binary_digest, "binary SHA-256 does not match manifest")
    collections: dict[str, dict[str, Any]] = {}
    for name in COLLECTION_FILES:
        digest, recorded_path = manifest[name]
        path = run / name
        require(path.resolve() == recorded_path.resolve(), f"{name}: manifest path differs")
        require(sha256(path) == digest, f"{name}: SHA-256 does not match manifest")
        collections[name] = {"bytes": path.stat().st_size, "sha256": digest}

    before_text = (run / "base-before.sha256").read_text()
    after_text = (run / "base-after.sha256").read_text()
    require(before_text == after_text, "base index changed while extension was built")
    base_line = re.fullmatch(r"([0-9a-f]{64})  (.+)\n", before_text)
    require(base_line is not None, "malformed base before/after digest")
    require(Path(base_line.group(2)).resolve() == (run / COLLECTION_FILES[0]).resolve(), "base digest path differs")
    require(base_line.group(1) == collections[COLLECTION_FILES[0]]["sha256"], "base digest differs from artifact manifest")

    require(collections["sez8-root.aicollection"] == collections["sez8-reverse.aicollection"],
            "root/reverse indexes are not byte deterministic")
    return {
        "binary": {"path": str(binary), "bytes": binary.stat().st_size, "sha256": binary_digest},
        "collections": collections,
        "root_reverse_byte_identical": True,
        "base_unchanged_after_extension": True,
        "artifact_manifest_sha256": sha256(run / "artifacts.sha256"),
    }


def validate_build_stdout(run: Path) -> dict[str, Any]:
    expected = {
        "root": ("sez8-root.aicollection", 8, 1_344_270, 1_941_890, 1_945_703, 31_233_997, 8_931_229),
        "reverse": ("sez8-reverse.aicollection", 8, 1_344_270, 1_941_890, 1_945_703, 31_233_997, 8_931_229),
        "base": ("sez4-base.aicollection", 4, 766_890, 968_922, 971_016, 17_188_103, 4_751_118),
        "extension": ("sez8-extension.aicollection", 4, 754_202, 972_968, 974_687, 17_000_065, 4_718_653),
    }
    result = {}
    pattern = re.compile(
        r"built (.+): (\d+) new archives, (\d+) segment junctions, (\d+) archive routes, "
        r"(\d+) chunk postings; (\d+) raw / (\d+) file bytes in ([0-9.]+)s\n"
    )
    for label, values in expected.items():
        match = pattern.fullmatch((run / f"{label}-build.stdout.txt").read_text())
        require(match is not None, f"{label}: build summary changed")
        path, *numbers = match.groups()
        observed = tuple(int(value) for value in numbers[:-1])
        require(Path(path).resolve() == (run / values[0]).resolve(), f"{label}: build output path differs")
        require(observed == values[1:], f"{label}: build counts differ: {observed}")
        require(float(numbers[-1]) >= 0, f"{label}: invalid internal build time")
        result[label] = {
            "new_archives": observed[0],
            "junction_rows": observed[1],
            "archive_routes": observed[2],
            "chunk_postings": observed[3],
            "raw_bytes": observed[4],
            "file_bytes": observed[5],
            "internal_seconds": float(numbers[-1]),
        }
    return result


def validate_archives(archives: list[dict[str, Any]]) -> None:
    require(type(archives) is list and len(archives) == len(SAMPLES), "inspect: archive count differs")
    require([archive.get("id") for archive in archives] == list(SAMPLES), "inspect: archives are not sorted A--H")
    for archive in archives:
        require_keys(archive, {"blake3", "bytes", "chunks", "id", "path"}, f"archive {archive.get('id')}")
        sample = archive["id"]
        expected_bytes, expected_digest, expected_chunks = SOURCE[sample]
        require((archive["bytes"], archive["blake3"], archive["chunks"])
                == (expected_bytes, expected_digest, expected_chunks), f"archive {sample}: identity differs")
        path = Path(archive["path"])
        require(path.as_posix().endswith(f"/sez-transcript-end-atlas-r1/donor-{sample}/sez.called.aie"),
                f"archive {sample}: unexpected logical path")
        require(path.is_file() and path.stat().st_size == expected_bytes, f"archive {sample}: source file missing or wrong size")


def validate_inspect(run: Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    root = load_json(run / "root-inspect.json")
    chain = load_json(run / "chain-inspect.json")
    verified = load_json(run / "root-inspect-verify.json")
    top_keys = {"archives", "chromosomes", "elapsed_seconds", "file_bytes", "format_version", "guard", "index", "layers", "path", "reference", "schema"}
    for label, value in (("root", root), ("chain", chain), ("verified root", verified)):
        require(type(value) is dict, f"{label} inspect is not an object")
        require_keys(value, top_keys, f"{label} inspect")
        require(value["schema"] == "gravlax.collection.v2" and value["format_version"] == 2,
                f"{label}: collection schema/version differs")
        require(value["reference"] == REFERENCE, f"{label}: reference identity differs")
        require(type(value["elapsed_seconds"]) in (int, float) and value["elapsed_seconds"] >= 0,
                f"{label}: invalid elapsed time")
        validate_archives(value["archives"])
    require(root["archives"] == chain["archives"] == verified["archives"], "root/chain archive manifests differ")
    require(root["chromosomes"] == chain["chromosomes"] == verified["chromosomes"], "root/chain chromosome dictionaries differ")
    require(len(root["chromosomes"]) == len(set(root["chromosomes"])) > 0, "chromosome dictionary is empty or duplicated")
    require(root["index"] == EXPECTED_INDEX["root"], "root index totals differ")
    require(chain["index"] == EXPECTED_INDEX["chain"], "chain index totals differ")
    require(verified["index"] == root["index"], "verified inspect index differs")
    require(root["guard"] == {"content_digest_recorded": True, "content_digest_verified": False,
                              "filesystem_identity": "size + nanosecond mtime + nanosecond ctime + device + inode"},
            "root guard differs")
    require(chain["guard"] == root["guard"], "chain guard differs")
    require(verified["guard"] == {**root["guard"], "content_digest_verified": True}, "verify-content guard differs")
    require(Path(root["path"]).resolve() == (run / "sez8-root.aicollection").resolve(), "root inspect path differs")
    require(Path(chain["path"]).resolve() == (run / "sez8-extension.aicollection").resolve(), "chain inspect path differs")
    require(Path(verified["path"]).resolve() == (run / "sez8-root.aicollection").resolve(), "verified inspect path differs")
    require(root["file_bytes"] == (run / "sez8-root.aicollection").stat().st_size, "root file byte accounting differs")
    chain_bytes = sum((run / name).stat().st_size for name in ("sez4-base.aicollection", "sez8-extension.aicollection"))
    require(chain["file_bytes"] == chain_bytes, "chain aggregate byte accounting differs")
    require(verified["file_bytes"] == root["file_bytes"], "verified root byte accounting differs")

    for label, layers, expected_counts, expected_paths in (
        ("root", root["layers"], [(8, 1_344_270, 282)], ["sez8-root.aicollection"]),
        ("chain", chain["layers"], [(4, 766_890, 276), (4, 754_202, 276)],
         ["sez4-base.aicollection", "sez8-extension.aicollection"]),
    ):
        require(len(layers) == len(expected_counts), f"{label}: layer count differs")
        for index, (layer, counts, name) in enumerate(zip(layers, expected_counts, expected_paths)):
            require_keys(layer, {"archives", "junction_rows", "junction_segments", "path", "root_digest"},
                         f"{label} layer {index}")
            require((layer["archives"], layer["junction_rows"], layer["junction_segments"]) == counts,
                    f"{label} layer {index}: totals differ")
            require(Path(layer["path"]).resolve() == (run / name).resolve(), f"{label} layer {index}: path differs")
            require(re.fullmatch(r"[0-9a-f]{64}", layer["root_digest"]) is not None,
                    f"{label} layer {index}: malformed root digest")
    root_normalized = copy.deepcopy(root)
    verified_normalized = copy.deepcopy(verified)
    for value in (root_normalized, verified_normalized):
        value.pop("elapsed_seconds")
        value.pop("guard")
    require(root_normalized == verified_normalized, "verified inspect changed semantic content")
    return root, chain, verified


def check_count_row(row: dict[str, Any], label: str) -> None:
    keys = ("include_only", "exclude_only", "both")
    require(all(type(row.get(key)) is int and row[key] >= 0 for key in keys), f"{label}: malformed counts")
    informative = row["include_only"] + row["exclude_only"]
    require(row.get("informative_umis") == informative, f"{label}: informative total differs")
    expected_usage = None if informative == 0 else row["include_only"] / informative
    require(row.get("usage_fraction") == expected_usage, f"{label}: usage fraction differs")


def validate_plan(plan: dict[str, Any], samples: list[dict[str, Any]], layers: int, kind: str,
                  expected_identity_bytes: int = 0) -> None:
    common = {"actual_archive_bytes_read", "archives_opened", "archives_pruned", "archives_total",
              "collection_layers", "collection_load_seconds", "identity_content_bytes_read",
              "identity_guard_seconds", "planned_compressed_bytes", "route_planning_seconds",
              "source_execution_seconds", "total_seconds", "unique_chunks_decoded"}
    if kind in {"junction", "jset"}:
        common |= {"archive_catalogue_sections_read", "archive_posting_sections_read", "support_bound_pruned"}
    else:
        common |= {"archive_chunk_index_sections_read"}
    if kind == "jset":
        common |= {"chunk_decode_reduction_fraction", "independent_chunk_decodes"}
    require_keys(plan, common, f"{kind} planning")
    require(plan["archives_total"] == 8 and plan["archives_opened"] + plan["archives_pruned"] == 8,
            f"{kind}: archive planning counts differ")
    require(plan["collection_layers"] == layers, f"{kind}: collection layer count differs")
    require(plan["identity_content_bytes_read"] == expected_identity_bytes,
            f"{kind}: identity-content byte accounting differs")
    require(plan["actual_archive_bytes_read"] == sum(row["actual_archive_bytes_read"] for row in samples),
            f"{kind}: actual archive byte sum differs")
    require(plan["planned_compressed_bytes"] == sum(row["planned_compressed_bytes"] for row in samples),
            f"{kind}: planned compressed byte sum differs")
    if kind == "jset":
        unique = sum(row["unique_chunks_decoded"] for row in samples)
        independent = sum(row["independent_chunk_decodes"] for row in samples)
        require(plan["unique_chunks_decoded"] == unique and plan["independent_chunk_decodes"] == independent,
                "jset: chunk accounting differs")
        expected_reduction = 0.0 if independent == 0 else 1 - unique / independent
        require(math.isclose(plan["chunk_decode_reduction_fraction"], expected_reduction, abs_tol=1e-15),
                "jset: chunk reduction fraction differs")
    else:
        require(plan["unique_chunks_decoded"] == sum(row["chunks_decoded"] for row in samples),
                f"{kind}: unique chunk accounting differs")
    for key in ("collection_load_seconds", "identity_guard_seconds", "route_planning_seconds",
                "source_execution_seconds", "total_seconds"):
        require(type(plan[key]) in (int, float) and plan[key] >= 0, f"{kind}: invalid {key}")
    components = sum(plan[key] for key in ("collection_load_seconds", "identity_guard_seconds",
                                            "route_planning_seconds", "source_execution_seconds"))
    require(math.isclose(plan["total_seconds"], components, rel_tol=1e-4, abs_tol=5e-6),
            f"{kind}: phase times do not sum to total")


def validate_top_cells(rows: list[dict[str, Any]], cells: int, umis: int, label: str) -> None:
    require(type(rows) is list and len(rows) == cells, f"{label}: top-cell cardinality differs")
    barcodes = []
    for row in rows:
        require_keys(row, {"barcode", "umis"}, f"{label} cell")
        require(type(row["barcode"]) is str and row["barcode"], f"{label}: bad barcode")
        require(type(row["umis"]) is int and row["umis"] > 0, f"{label}: bad UMI count")
        barcodes.append(row["barcode"])
    require(len(barcodes) == len(set(barcodes)), f"{label}: duplicate cell barcode")
    require(sum(row["umis"] for row in rows) == umis, f"{label}: cell UMIs do not sum")


def validate_junction(value: dict[str, Any], archives: list[dict[str, Any]], layers: int,
                      label: str, expected_identity_bytes: int = 0,
                      support_pruned: bool = False) -> dict[str, Any]:
    require_keys(value, {"acceptor", "chrom", "collection_schema", "coordinates", "donor", "explain",
                         "min_support", "planning", "samples", "schema", "support_upper_bound", "totals"}, label)
    require(value["schema"] == "gravlax.collection.junction.v1" and
            value["collection_schema"] == "gravlax.collection.v2", f"{label}: schema differs")
    expected_kind = "dense" if value["donor"] == 129_925_157 else (
        "sparse" if value["donor"] == 129_861_033 else "absent")
    require(f"{value['chrom']}:{value['donor']}-{value['acceptor']}" == LOCI[expected_kind],
            f"{label}: locus differs")
    require(value["totals"] == ({"cells": 0, "umis": 0} if support_pruned else EXPECTED_TOTALS[expected_kind]),
            f"{label}: totals differ")
    expected_support = {"dense": 210, "sparse": 4, "absent": 0}[expected_kind]
    require(value["support_upper_bound"] == expected_support, f"{label}: support upper bound differs")
    require(value["min_support"] == (211 if support_pruned else 0), f"{label}: minimum support differs")
    samples = value["samples"]
    require([row.get("sample") for row in samples] == list(SAMPLES), f"{label}: sample order differs")
    for row, archive in zip(samples, archives):
        require_keys(row, {"actual_archive_bytes_read", "archive", "cells", "chunks_decoded",
                           "planned_compressed_bytes", "present", "sample", "supporting_children",
                           "top_cells", "umis"}, f"{label}/{row.get('sample')}")
        require(row["archive"] == archive["path"], f"{label}/{row['sample']}: archive path differs")
        require(all(type(row[key]) is int and row[key] >= 0 for key in
                    ("actual_archive_bytes_read", "cells", "chunks_decoded", "planned_compressed_bytes",
                     "supporting_children", "umis")), f"{label}/{row['sample']}: negative or malformed count")
        validate_top_cells(row["top_cells"], row["cells"], row["umis"], f"{label}/{row['sample']}")
    if not support_pruned:
        require(sum(row["cells"] for row in samples) == value["totals"]["cells"] and
                sum(row["umis"] for row in samples) == value["totals"]["umis"], f"{label}: sample totals do not sum")
        require(sum(row["supporting_children"] for row in samples) == expected_support,
                f"{label}: sample support does not sum")
    validate_plan(value["planning"], samples, layers, "junction", expected_identity_bytes)
    require(value["planning"]["support_bound_pruned"] is support_pruned, f"{label}: support-pruned flag differs")
    explain = value["explain"]
    require(len(explain) == 8, f"{label}: explanation count differs")
    for detail, row in zip(explain, samples):
        require_keys(detail, {"actual_archive_bytes_read", "decision", "planned_compressed_bytes",
                              "posting_chunks", "sample"}, f"{label} explain/{row['sample']}")
        require(detail["sample"] == row["sample"] and detail["actual_archive_bytes_read"] == row["actual_archive_bytes_read"]
                and detail["planned_compressed_bytes"] == row["planned_compressed_bytes"]
                and detail["posting_chunks"] == row["chunks_decoded"], f"{label}: explanation differs from sample row")
        expected_decision = "prune_support_bound" if support_pruned else ("open" if row["present"] else "prune_absent")
        require(detail["decision"] == expected_decision, f"{label}/{row['sample']}: route decision differs")
    return {
        "locus": LOCI[expected_kind],
        "support_upper_bound": expected_support,
        "totals": value["totals"],
        "samples": [
            {key: row[key] for key in ("sample", "present", "supporting_children", "umis", "cells",
                                       "chunks_decoded", "top_cells")}
            for row in samples
        ],
    }


def validate_region(value: dict[str, Any], archives: list[dict[str, Any]], layers: int, label: str) -> dict[str, Any]:
    require_keys(value, {"anchor_semantics", "chrom", "collection_schema", "coordinates", "end", "explain",
                         "planning", "samples", "schema", "start", "totals"}, label)
    require(value["schema"] == "gravlax.collection.region.v1" and
            value["collection_schema"] == "gravlax.collection.v2", f"{label}: schema differs")
    require(value["anchor_semantics"] is True and f"{value['chrom']}:{value['start']}-{value['end']}" == LOCI["region"],
            f"{label}: region differs")
    require(value["totals"] == EXPECTED_TOTALS["region"], f"{label}: totals differ")
    samples = value["samples"]
    require([row.get("sample") for row in samples] == list(SAMPLES), f"{label}: sample order differs")
    for row, archive in zip(samples, archives):
        require_keys(row, {"actual_archive_bytes_read", "archive", "cells", "chunks_decoded", "molecules",
                           "planned_compressed_bytes", "present", "sample", "supporting_children", "top_cells", "umis"},
                     f"{label}/{row.get('sample')}")
        require(row["archive"] == archive["path"] and row["present"] is True and row["supporting_children"] == 0,
                f"{label}/{row['sample']}: identity/presence differs")
        require(all(type(row[key]) is int and row[key] >= 0 for key in
                    ("actual_archive_bytes_read", "cells", "chunks_decoded", "molecules",
                     "planned_compressed_bytes", "umis")), f"{label}/{row['sample']}: malformed count")
        validate_top_cells(row["top_cells"], row["cells"], row["umis"], f"{label}/{row['sample']}")
    for key in ("cells", "molecules", "umis"):
        require(sum(row[key] for row in samples) == value["totals"][key], f"{label}: {key} total does not sum")
    validate_plan(value["planning"], samples, layers, "region")
    require(len(value["explain"]) == 8, f"{label}: explanation count differs")
    for detail, row in zip(value["explain"], samples):
        require_keys(detail, {"actual_archive_bytes_read", "chunks", "decision", "planned_compressed_bytes", "sample"},
                     f"{label} explain/{row['sample']}")
        require(detail == {"actual_archive_bytes_read": row["actual_archive_bytes_read"],
                           "chunks": row["chunks_decoded"], "decision": "open",
                           "planned_compressed_bytes": row["planned_compressed_bytes"], "sample": row["sample"]},
                f"{label}/{row['sample']}: explanation differs")
    return {
        "locus": LOCI["region"],
        "totals": value["totals"],
        "samples": [{key: row[key] for key in ("sample", "molecules", "umis", "cells", "chunks_decoded",
                                                "top_cells")} for row in samples],
    }


def validate_jset(value: dict[str, Any], archives: list[dict[str, Any]], layers: int, label: str) -> dict[str, Any]:
    require_keys(value, {"collection_schema", "coordinates", "explain", "junctions", "min_support", "planning",
                         "samples", "schema", "semantics", "totals"}, label)
    require(value["schema"] == "gravlax.collection.jset.v1" and
            value["collection_schema"] == "gravlax.collection.v2", f"{label}: schema differs")
    require(value["min_support"] == 0 and value["totals"] == EXPECTED_TOTALS["jset"], f"{label}: totals/minimum differ")
    require(value["semantics"] == {"both_in_usage_denominator": False,
                                    "class_categories": ["include_only", "exclude_only", "both"],
                                    "informative_umis": "include_only + exclude_only",
                                    "usage_fraction": "include_only / informative_umis"}, f"{label}: semantics differ")
    require(value["junctions"] == [
        {"locus": LOCI["dense"], "passes_min_support": True, "present_samples": 8,
         "side": "include", "support_upper_bound": 210},
        {"locus": LOCI["sparse"], "passes_min_support": True, "present_samples": 2,
         "side": "exclude", "support_upper_bound": 4},
    ], f"{label}: junction catalogue differs")
    samples = value["samples"]
    require([row.get("sample") for row in samples] == list(SAMPLES), f"{label}: sample order differs")
    for row, archive in zip(samples, archives):
        require_keys(row, {"actual_archive_bytes_read", "archive", "cells", "independent_chunk_decodes",
                           "planned_compressed_bytes", "present_components", "sample", "top_cells", "totals",
                           "unique_chunks_decoded"}, f"{label}/{row.get('sample')}")
        require(row["archive"] == archive["path"], f"{label}/{row['sample']}: archive path differs")
        check_count_row(row["totals"], f"{label}/{row['sample']}")
        require(type(row["top_cells"]) is list and len(row["top_cells"]) == row["cells"],
                f"{label}/{row['sample']}: cell count differs")
        seen = set()
        sums = {key: 0 for key in ("include_only", "exclude_only", "both")}
        for cell in row["top_cells"]:
            require_keys(cell, {"barcode", "both", "exclude_only", "include_only", "informative_umis", "usage_fraction"},
                         f"{label}/{row['sample']} cell")
            require(cell["barcode"] not in seen, f"{label}/{row['sample']}: duplicate barcode")
            seen.add(cell["barcode"])
            check_count_row(cell, f"{label}/{row['sample']}/{cell['barcode']}")
            require(cell["informative_umis"] + cell["both"] > 0, f"{label}: explicit zero cell")
            for key in sums:
                sums[key] += cell[key]
        require(all(sums[key] == row["totals"][key] for key in sums), f"{label}/{row['sample']}: cell totals differ")
    aggregate = {key: sum(row["totals"][key] for row in samples) for key in ("include_only", "exclude_only", "both")}
    require(all(aggregate[key] == value["totals"][key] for key in aggregate), f"{label}: sample totals do not sum")
    check_count_row(value["totals"], f"{label} total")
    validate_plan(value["planning"], samples, layers, "jset")
    require(value["planning"]["support_bound_pruned"] is False, f"{label}: unexpected support pruning")
    require(len(value["explain"]) == 8, f"{label}: explanation count differs")
    for detail, row in zip(value["explain"], samples):
        require_keys(detail, {"actual_archive_bytes_read", "decision", "planned_compressed_bytes",
                              "present_components", "sample", "unique_chunks"}, f"{label} explain/{row['sample']}")
        require(detail == {"actual_archive_bytes_read": row["actual_archive_bytes_read"], "decision": "open",
                           "planned_compressed_bytes": row["planned_compressed_bytes"],
                           "present_components": row["present_components"], "sample": row["sample"],
                           "unique_chunks": row["unique_chunks_decoded"]}, f"{label}/{row['sample']}: explanation differs")
    return {
        "junctions": value["junctions"],
        "totals": value["totals"],
        "samples": [{key: row[key] for key in ("sample", "present_components", "cells", "totals", "top_cells",
                                                "unique_chunks_decoded", "independent_chunk_decodes")}
                    for row in samples],
    }


def validate_naive_junction(path: Path, collection_row: dict[str, Any], locus: str, explicit_absence: bool) -> None:
    value = load_json(path)
    if not collection_row["present"]:
        require(value == {"present": False, "umis": 0, "cells": 0}, f"{path}: absence marker differs")
        return
    require_keys(value, {"acceptor", "cell_rows", "cell_rows_truncated", "cells", "chrom", "coordinates",
                         "donor", "posting_chunks", "schema", "supporting_children", "umis"}, str(path))
    chrom, bounds = locus.split(":")
    donor, acceptor = (int(part) for part in bounds.split("-"))
    require(value["schema"] == "gravlax.query.junction.v1" and value["chrom"] == chrom
            and value["donor"] == donor and value["acceptor"] == acceptor, f"{path}: locus/schema differs")
    require(value["cell_rows_truncated"] is False, f"{path}: naive cell rows were truncated")
    require((value["umis"], value["cells"], value["supporting_children"], value["posting_chunks"], value["cell_rows"])
            == (collection_row["umis"], collection_row["cells"], collection_row["supporting_children"],
                collection_row["chunks_decoded"], collection_row["top_cells"]), f"{path}: collection/naive result differs")
    require(not explicit_absence, f"{path}: expected explicit absence but query is present")


def validate_naive(run: Path, root_values: dict[str, dict[str, Any]]) -> int:
    by_kind = {kind: {row["sample"]: row for row in root_values[kind]["samples"]}
               for kind in ("dense", "sparse", "absent")}
    checks = 0
    absent_marker = {"present": False, "umis": 0, "cells": 0}
    for sample in SAMPLES:
        for kind in ("dense", "sparse", "absent"):
            path = run / "naive" / f"{sample}-{kind}.json"
            row = by_kind[kind][sample]
            validate_naive_junction(path, row, LOCI[kind], kind == "absent")
            checks += 1
            if kind == "sparse":
                stderr = (run / "naive" / f"{sample}-sparse.stderr.txt").read_text()
                if row["present"]:
                    require("not present in the archive" not in stderr and f"junction {LOCI[kind]}: open" in stderr,
                            f"sample {sample}: successful sparse marker differs")
                else:
                    require(load_json(path) == absent_marker and
                            f"junction {LOCI[kind]} not present in the archive" in stderr,
                            f"sample {sample}: sparse absence marker differs")
            if kind == "absent":
                require((run / "naive" / f"{sample}-absent.stdout.txt").read_bytes() == b"",
                        f"sample {sample}: absent query unexpectedly wrote stdout")
                stderr = (run / "naive" / f"{sample}-absent.stderr.txt").read_text()
                require(load_json(path) == absent_marker and f"junction {LOCI[kind]} not present in the archive" in stderr,
                        f"sample {sample}: absent marker differs")

        region = load_json(run / "naive" / f"{sample}-region.json")
        region_row = next(row for row in root_values["region"]["samples"] if row["sample"] == sample)
        require_keys(region, {"anchor_semantics", "cell_rows", "cell_rows_truncated", "cells", "chrom", "chunks_decoded",
                              "coordinates", "end", "molecules", "schema", "scope", "start", "umis"},
                     f"naive region/{sample}")
        require(region["schema"] == "gravlax.query.region.v1" and region["anchor_semantics"] is True and
                f"{region['chrom']}:{region['start']}-{region['end']}" == LOCI["region"], f"naive region/{sample}: locus differs")
        require(region["cell_rows"] == [] and region["cell_rows_truncated"] is True,
                f"naive region/{sample}: --top 0 behavior differs")
        require((region["molecules"], region["umis"], region["cells"], region["chunks_decoded"])
                == (region_row["molecules"], region_row["umis"], region_row["cells"], region_row["chunks_decoded"]),
                f"naive region/{sample}: collection totals differ")
        checks += 1

        jset = load_json(run / "naive" / f"{sample}-jset.json")
        jset_row = next(row for row in root_values["jset"]["samples"] if row["sample"] == sample)
        require_keys(jset, {"cell_rows", "cell_rows_truncated", "cells", "coordinates", "exclusion_junctions",
                            "inclusion_junctions", "planning", "schema", "scope", "semantics", "totals"},
                     f"naive jset/{sample}")
        require(jset["schema"] == "gravlax.query.jset.v1" and jset["cell_rows_truncated"] is False,
                f"naive jset/{sample}: schema/truncation differs")
        require((jset["cells"], jset["totals"], jset["cell_rows"])
                == (jset_row["cells"], jset_row["totals"], jset_row["top_cells"]),
                f"naive jset/{sample}: collection result differs")
        present_components = sum(row["present"] for row in jset["inclusion_junctions"] + jset["exclusion_junctions"])
        require(present_components == jset_row["present_components"], f"naive jset/{sample}: presence count differs")
        require(jset["planning"]["unique_chunk_decodes"] == jset_row["unique_chunks_decoded"] and
                jset["planning"]["independent_chunk_decodes"] == jset_row["independent_chunk_decodes"],
                f"naive jset/{sample}: chunk planning differs")
        checks += 1
    require(checks == 40, "naive comparison count differs")
    return checks


def validate_timings(run: Path, binary: Path) -> dict[str, Any]:
    names = sorted(path.stem.removesuffix(".time") for path in run.glob("*.time.txt"))
    expected = sorted({
        "root-build", "reverse-build", "base-build", "extension-build",
        "root-inspect", "chain-inspect", "root-inspect-verify",
        *(f"{arm}-{kind}" for arm in ("root", "chain") for kind in QUERY_KINDS),
        "root-support-pruned", "root-verify-dense",
        *(f"naive-{kind}" for kind in QUERY_KINDS),
    })
    require(names == expected, "timing record set differs")
    result = {name: read_time(run / f"{name}.time.txt") for name in names}
    for name, record in result.items():
        if name.startswith("naive-"):
            kind = name.removeprefix("naive-")
            tokens = shlex.split(record["command"])
            require(tokens[0].endswith("scripts/206_benchmark_collection_index.sh"), f"{name}: driver path differs")
            require(tokens[1:] == ["__naive-kind", kind, str(run / "naive")],
                    f"{name}: aggregate naive command differs")
        else:
            tokens = shlex.split(record["command"])
            require(Path(tokens[0]).resolve() == binary.resolve(), f"{name}: timed binary differs")

    def exact(label: str, tail: list[str]) -> None:
        observed = shlex.split(result[label]["command"])
        require(observed == [str(binary), *tail], f"{label}: timed command differs")

    root_path = str(run / "sez8-root.aicollection")
    chain_path = str(run / "sez8-extension.aicollection")
    query_tail = {
        "dense": ["collection", "junction", "{}", LOCI["dense"], "--json", "--explain", "--top", "0"],
        "sparse": ["collection", "junction", "{}", LOCI["sparse"], "--json", "--explain", "--top", "0"],
        "absent": ["collection", "junction", "{}", LOCI["absent"], "--json", "--explain", "--top", "0"],
        "region": ["collection", "region", "{}", LOCI["region"], "--json", "--explain", "--top", "0"],
        "jset": ["collection", "jset", "{}", "--include", LOCI["dense"], "--exclude", LOCI["sparse"],
                 "--json", "--explain", "--top", "0"],
    }
    for arm, collection_path in (("root", root_path), ("chain", chain_path)):
        for kind in QUERY_KINDS:
            exact(f"{arm}-{kind}", [collection_path if token == "{}" else token for token in query_tail[kind]])
    exact("root-support-pruned", ["collection", "junction", root_path, LOCI["dense"], "--min-support", "211",
                                  "--json", "--explain", "--top", "0"])
    exact("root-verify-dense", ["collection", "junction", root_path, LOCI["dense"], "--verify-content",
                               "--json", "--explain", "--top", "0"])
    exact("root-inspect", ["collection", "inspect", root_path])
    exact("chain-inspect", ["collection", "inspect", chain_path])
    exact("root-inspect-verify", ["collection", "inspect", root_path, "--verify-content"])

    def build_tail(samples: tuple[str, ...], output: str, base: str | None = None) -> list[str]:
        tail = ["collection", "build"]
        if base is not None:
            tail.extend(["--base", base])
        project_root = run.parents[2]
        for sample in samples:
            archive = project_root / "runs" / "post-v1" / "sez-transcript-end-atlas-r1" / f"donor-{sample}" / "sez.called.aie"
            tail.extend(["--sample", f"{sample}={archive}", "--source-digest", f"{sample}={SOURCE[sample][1]}"])
        tail.extend(["--out", str(run / output)])
        return tail

    exact("root-build", build_tail(SAMPLES, "sez8-root.aicollection"))
    exact("reverse-build", build_tail(tuple(reversed(SAMPLES)), "sez8-reverse.aicollection"))
    exact("base-build", build_tail(SAMPLES[:4], "sez4-base.aicollection"))
    exact("extension-build", build_tail(SAMPLES[4:], "sez8-extension.aicollection",
                                        str(run / "sez4-base.aicollection")))
    return result


def query_measurement(value: dict[str, Any], root_time: dict[str, Any], chain_time: dict[str, Any],
                      naive_time: dict[str, Any], source_bytes: int) -> dict[str, Any]:
    planning = value["planning"]
    speedup = None if root_time["wall_seconds"] == 0 else naive_time["wall_seconds"] / root_time["wall_seconds"]
    return {
        "totals": value["totals"],
        "root": {
            "external_wall_seconds": root_time["wall_seconds"],
            "max_rss_kib": root_time["max_rss_kib"],
            "internal_total_seconds": planning["total_seconds"],
            "archives_opened": planning["archives_opened"],
            "archives_pruned": planning["archives_pruned"],
            "unique_chunks_decoded": planning["unique_chunks_decoded"],
            "planned_compressed_bytes": planning["planned_compressed_bytes"],
            "actual_archive_bytes_read": planning["actual_archive_bytes_read"],
            "actual_fraction_of_source_archives": planning["actual_archive_bytes_read"] / source_bytes,
        },
        "chain": {"external_wall_seconds": chain_time["wall_seconds"], "max_rss_kib": chain_time["max_rss_kib"]},
        "naive_eight_archive_aggregate": {"external_wall_seconds": naive_time["wall_seconds"],
                                          "max_rss_kib": naive_time["max_rss_kib"]},
        "naive_over_root_wall_ratio": speedup,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--gravlax-commit", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run = args.run_dir.resolve()
    binary = args.binary.resolve()
    require(run.name == EXPECTED_RUN, f"only the locked {EXPECTED_RUN} run may be frozen")
    require(args.gravlax_commit == EXPECTED_COMMIT, "unexpected Gravlax commit")
    require(binary.is_file(), "binary does not exist")

    file_set = validate_file_set(run)
    checksums = validate_checksums(run, binary)
    build = validate_build_stdout(run)
    root_inspect, chain_inspect, verified_inspect = validate_inspect(run)
    archives = root_inspect["archives"]
    source_bytes = sum(archive["bytes"] for archive in archives)
    require(source_bytes == 1_326_036_691, "source archive byte total differs")
    project_root = run.parents[2]
    logical_archives = [
        {**archive, "path": str(Path(archive["path"]).resolve().relative_to(project_root))}
        for archive in archives
    ]

    root_raw = {kind: load_json(run / f"root-{kind}.json") for kind in QUERY_KINDS}
    chain_raw = {kind: load_json(run / f"chain-{kind}.json") for kind in QUERY_KINDS}
    root_logic = {
        "dense": validate_junction(root_raw["dense"], archives, 1, "root dense"),
        "sparse": validate_junction(root_raw["sparse"], archives, 1, "root sparse"),
        "absent": validate_junction(root_raw["absent"], archives, 1, "root absent"),
        "region": validate_region(root_raw["region"], archives, 1, "root region"),
        "jset": validate_jset(root_raw["jset"], archives, 1, "root jset"),
    }
    chain_logic = {
        "dense": validate_junction(chain_raw["dense"], archives, 2, "chain dense"),
        "sparse": validate_junction(chain_raw["sparse"], archives, 2, "chain sparse"),
        "absent": validate_junction(chain_raw["absent"], archives, 2, "chain absent"),
        "region": validate_region(chain_raw["region"], archives, 2, "chain region"),
        "jset": validate_jset(chain_raw["jset"], archives, 2, "chain jset"),
    }
    require(root_logic == chain_logic, "fresh root and immutable base+extension chain differ logically")
    naive_checks = validate_naive(run, root_logic)

    pruned = load_json(run / "root-support-pruned.json")
    pruned_logic = validate_junction(pruned, archives, 1, "support-pruned dense", support_pruned=True)
    require(pruned["planning"]["archives_opened"] == 0 and pruned["planning"]["archives_pruned"] == 8
            and pruned["planning"]["actual_archive_bytes_read"] == 0
            and pruned["planning"]["planned_compressed_bytes"] == 0
            and all(row["present"] for row in pruned["samples"]), "support-bound pruning did not avoid all source reads")

    verified_dense = load_json(run / "root-verify-dense.json")
    verified_logic = validate_junction(verified_dense, archives, 1, "verified dense",
                                       expected_identity_bytes=source_bytes)
    require(verified_logic == root_logic["dense"], "verify-content changed the dense result")
    require(verified_dense["planning"]["identity_content_bytes_read"] == source_bytes,
            "verify-content did not hash every source byte exactly once")

    timings = validate_timings(run, binary)
    root_bytes = root_inspect["file_bytes"]
    chain_bytes = chain_inspect["file_bytes"]
    root_fraction = root_bytes / source_bytes
    chain_fraction = chain_bytes / source_bytes
    require(root_fraction < 0.05 and chain_fraction < 0.05, "collection sidecar exceeds the 5% gate")

    query_result = {
        kind: query_measurement(root_raw[kind], timings[f"root-{kind}"], timings[f"chain-{kind}"],
                                timings[f"naive-{kind}"], source_bytes)
        for kind in QUERY_KINDS
    }
    result = {
        "schema": "gravlax.collection-index-benchmark.v1",
        "status": "PASS",
        "identity": {
            "gravlax_commit": args.gravlax_commit,
            "benchmark_run": f"runs/post-v1/{EXPECTED_RUN}",
            "binary": {**checksums["binary"], "path": str(binary.relative_to(project_root))},
            "validated_run_file_set": file_set,
            "artifact_manifest_sha256": checksums["artifact_manifest_sha256"],
        },
        "inputs": {
            "archives": logical_archives,
            "archive_count": len(archives),
            "archive_bytes": source_bytes,
            "archive_chunks": sum(archive["chunks"] for archive in archives),
            "reference": REFERENCE,
        },
        "collection": {
            "format_version": root_inspect["format_version"],
            "root": {**checksums["collections"]["sez8-root.aicollection"],
                     "root_digest": root_inspect["layers"][0]["root_digest"],
                     "fraction_of_source_archives": root_fraction,
                     "index": root_inspect["index"]},
            "reverse": checksums["collections"]["sez8-reverse.aicollection"],
            "base": {**checksums["collections"]["sez4-base.aicollection"],
                     "root_digest": chain_inspect["layers"][0]["root_digest"]},
            "extension": {**checksums["collections"]["sez8-extension.aicollection"],
                          "root_digest": chain_inspect["layers"][1]["root_digest"]},
            "chain": {"bytes": chain_bytes, "fraction_of_source_archives": chain_fraction,
                      "layers": 2, "index": chain_inspect["index"]},
            "build_summaries": build,
        },
        "queries": query_result,
        "support_pruning": {
            "locus": pruned_logic["locus"],
            "support_upper_bound": pruned_logic["support_upper_bound"],
            "min_support": 211,
            "archives_opened": 0,
            "archives_pruned": 8,
            "actual_archive_bytes_read": 0,
            "external_wall_seconds": timings["root-support-pruned"]["wall_seconds"],
            "max_rss_kib": timings["root-support-pruned"]["max_rss_kib"],
        },
        "content_verification": {
            "inspect_guard_verified": verified_inspect["guard"]["content_digest_verified"],
            "query_identity_content_bytes_read": source_bytes,
            "one_full_pass_over_source_archives": True,
            "result_unchanged": True,
            "external_wall_seconds": timings["root-verify-dense"]["wall_seconds"],
            "max_rss_kib": timings["root-verify-dense"]["max_rss_kib"],
        },
        "timings": {
            name: {key: value for key, value in record.items() if key != "command"}
            for name, record in timings.items()
        },
        "gates": {
            "root_reverse_byte_deterministic": checksums["root_reverse_byte_identical"],
            "base_unchanged_by_extension": checksums["base_unchanged_after_extension"],
            "root_chain_logically_identical_for_five_queries": True,
            "forty_per_sample_results_equal_naive": naive_checks == 40,
            "explicit_absence_markers_validated": True,
            "support_bound_prunes_all_archives": True,
            "verify_content_reads_each_source_byte_once": True,
            "root_sidecar_below_five_percent": root_fraction < 0.05,
            "chain_sidecars_below_five_percent": chain_fraction < 0.05,
            "all_twenty_four_timed_commands_exit_zero": len(timings) == 24,
            "expected_scientific_totals": True,
        },
        "execution_history": {
            "r1_r2": "failed benchmark-driver attempts; not scientific arms",
            "r3_r4": "successful pre-final dry runs; not scientific arms",
            "r5": "locked final benchmark with matched aggregate naive timings",
        },
    }
    require(all(result["gates"].values()), f"collection-index gate failed: {result['gates']}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    try:
        main()
    except ValidationError as error:
        raise SystemExit(f"validation failed: {error}") from error
