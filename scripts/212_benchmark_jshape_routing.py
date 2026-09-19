#!/usr/bin/env python3
"""Run the frozen Gate-B local-shape-routing campaign without reducing its verdict.

The frozen adapter requires the final Gravlax CLI and JSON contracts. A promotable run additionally
requires clean code and scripts commits, the frozen Gate-A result and 96-junction panel, the same
binary for both arms, exactly eight alternating warm blocks, and complete independent fixture-hook
coverage. The artifact manifest is written last.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import socket
import subprocess
from pathlib import Path
from typing import Any, Mapping, Sequence

import zstandard

from archive_root_gate_common import (
    GateError,
    load_json,
    parse_archive,
    require,
    run_timed,
    self_test_blake3,
    sha256_file,
    write_artifact_manifest,
    write_json,
)
from jshape_routing_gate_common import (
    DATE,
    DEFAULT_ADAPTER,
    EXTERNAL_FIXTURE_CASES,
    EXPECTED_TOTALS,
    GATE_A_RESULT_COMMIT,
    GATE_A_RESULT_SHA256,
    PANEL_METADATA_SHA256,
    PANEL_SHA256,
    REQUIRED_FIXTURE_CASES,
    RUST_FIXTURE_TESTS,
    ROOT_DOMAIN,
    ROOTED_ARCHIVES,
    SAMPLES,
    THREADS,
    WARM_BLOCKS,
    adapter_is_frozen_final,
    build_command,
    build_schedule,
    fixed_cases,
    inspect_command,
    is_promotable_run,
    load_frozen_panel,
    parse_sample,
    planned_output_directories,
    prepare_output_tree,
    query_command,
    scientific_projection,
    validate_adapter,
)


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


def host_memory_kib() -> int | None:
    path = Path("/proc/meminfo")
    if not path.is_file():
        return None
    match = re.search(r"^MemTotal:\s+(\d+)\s+kB$", path.read_text(), re.MULTILINE)
    return int(match.group(1)) if match else None


def write_build_schedule(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    lines = ["case_id\tblock\tposition\tarm\tlabel"]
    lines.extend(
        f"{row['case_id']}\t{row['block']}\t{row['position']}\t{row['arm']}\t{row['label']}"
        for row in rows
    )
    path.write_text("\n".join(lines) + "\n")


def write_query_schedule(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    cases: Mapping[str, Any],
) -> None:
    lines = ["case_id\tkind\tlocus\tstratum\tblock\tposition\tarm\tlabel"]
    for row in rows:
        case = cases[row["case_id"]]
        lines.append(
            "\t".join(
                map(
                    str,
                    (
                        case.case_id,
                        case.kind,
                        case.locus,
                        case.stratum or "",
                        row["block"],
                        row["position"],
                        row["arm"],
                        row["label"],
                    ),
                )
            )
        )
    path.write_text("\n".join(lines) + "\n")


def require_clean_success(command_dir: Path, *, json_stdout: bool = True) -> Any | None:
    require((command_dir / "stderr.txt").stat().st_size == 0, f"{command_dir}: unexpected stderr")
    if not json_stdout:
        require((command_dir / "stdout.txt").stat().st_size == 0, f"{command_dir}: unexpected stdout")
        return None
    value = load_json(command_dir / "stdout.txt")
    require(type(value) is dict, f"{command_dir}: stdout is not one JSON object")
    return value


def validate_fixture_result(path: Path, binary_sha256: str, route_sha256: str) -> dict[str, str]:
    value = load_json(path)
    require(
        type(value) is dict
        and set(value)
        == {
            "schema", "binary_sha256", "route_collection_sha256", "panel_sha256",
            "code_commit", "tools", "rust_tests", "cases",
        },
        f"{path}: fixture result fields differ",
    )
    require(value["schema"] == "gravlax.jshape-route-adversarial.v2", "fixture schema differs")
    require(value["binary_sha256"] == binary_sha256, "fixture binary differs")
    require(value["route_collection_sha256"] == route_sha256, "fixture route collection differs")
    require(value["panel_sha256"] == PANEL_SHA256, "fixture panel differs")
    protocol = load_json(path.parents[2] / "protocol.json")
    require(value["code_commit"] == protocol["gravlax"]["commit"], "fixture code commit differs")
    require(type(value["tools"]) is dict and set(value["tools"]) == {"python", "zstandard_python"}, "fixture tools differ")
    rust_tests = value["rust_tests"]
    require(type(rust_tests) is list and [row.get("test_name") for row in rust_tests if type(row) is dict] == list(RUST_FIXTURE_TESTS), "fixture Rust tests differ")
    code_repo = Path(protocol["gravlax"]["repo"]).resolve()
    for row in rust_tests:
        require(
            set(row) == {"test_name", "command", "exit_status", "stdout", "stderr", "source", "source_bytes", "source_sha256", "cases"},
            "fixture Rust-test fields differ",
        )
        require(row["exit_status"] == 0 and tuple(row["cases"]) == RUST_FIXTURE_TESTS[row["test_name"]], "fixture Rust-test result differs")
        expected_command = ["cargo", "test", "-p", "aie", "--bin", "aie", row["test_name"], "--", "--exact", "--test-threads=1"]
        require(row["command"] == expected_command, "fixture Rust-test command differs")
        source = code_repo / "crates/aie/src" / ("shaperoute.rs" if row["test_name"].startswith("shaperoute::") else "collectioncmd.rs")
        require(row["source"] == str(source) and row["source_bytes"] == source.stat().st_size and row["source_sha256"] == sha256_file(source), "fixture Rust-test source differs")
        for stream in (row["stdout"], row["stderr"]):
            require(type(stream) is str and not Path(stream).is_absolute() and (path.parent / stream).is_file(), "fixture Rust-test stream differs")
        require(re.search(r"test result: ok\. 1 passed;", (path.parent / row["stdout"]).read_text()) is not None, "fixture Rust test did not run exactly once")
    cases = value["cases"]
    require(type(cases) is list and cases, "fixture case list is empty")
    observed: dict[str, str] = {}
    expected_fields = {
        "id", "expect", "evidence_kind", "test_name", "command", "exit_status", "stdout", "stderr", "fixture_bytes",
        "fixture_sha256", "mutation", "rejected_before_scientific_output",
        "reconstruction_exact",
    }
    for case in cases:
        require(type(case) is dict and set(case) == expected_fields, "fixture case fields differ")
        case_id, expect = case["id"], case["expect"]
        require(type(case_id) is str and case_id not in observed, "duplicate fixture case")
        require(expect in {"accept", "reject"}, f"{case_id}: invalid expectation")
        observed[case_id] = expect
        require(type(case["command"]) is list and all(type(token) is str for token in case["command"]), f"{case_id}: command differs")
        require(type(case["exit_status"]) is int, f"{case_id}: exit status differs")
        require(type(case["fixture_bytes"]) is int and case["fixture_bytes"] >= 0, f"{case_id}: fixture size differs")
        require(re.fullmatch(r"[0-9a-f]{64}", str(case["fixture_sha256"])) is not None, f"{case_id}: fixture digest differs")
        require(type(case["mutation"]) is str and case["mutation"], f"{case_id}: mutation description missing")
        require(case["evidence_kind"] in {"rust_unit_test", "cli_mutation"}, f"{case_id}: evidence kind differs")
        for stream in ("stdout", "stderr"):
            relative = case[stream]
            require(type(relative) is str and not Path(relative).is_absolute(), f"{case_id}: unsafe stream path")
            stream_path = path.parent / relative
            require(stream_path.is_file() and stream_path.resolve().is_relative_to(path.parent.resolve()), f"{case_id}: stream missing")
        stdout = path.parent / case["stdout"]
        if case["evidence_kind"] == "rust_unit_test":
            require(case_id not in EXTERNAL_FIXTURE_CASES and case["test_name"] in RUST_FIXTURE_TESTS, f"{case_id}: Rust evidence mapping differs")
            require(case_id in RUST_FIXTURE_TESTS[case["test_name"]] and case["exit_status"] == 0, f"{case_id}: Rust predicate failed")
            require(case["rejected_before_scientific_output"] is (expect == "reject"), f"{case_id}: Rust rejection predicate differs")
            require(case["reconstruction_exact"] is (expect == "accept"), f"{case_id}: Rust exactness predicate differs")
            require(re.search(r"test result: ok\. 1 passed;", stdout.read_text()) is not None, f"{case_id}: Rust test evidence differs")
        elif expect == "reject":
            require(case_id in EXTERNAL_FIXTURE_CASES and case["test_name"] is None, f"{case_id}: external evidence mapping differs")
            require(case["exit_status"] != 0, f"{case_id}: invalid route accepted")
            require(case["rejected_before_scientific_output"] is True, f"{case_id}: partial output flag differs")
            require(case["reconstruction_exact"] is False, f"{case_id}: rejected route marked exact")
            require(stdout.stat().st_size == 0, f"{case_id}: rejected route emitted scientific stdout")
        else:
            require(case["exit_status"] == 0, f"{case_id}: valid fixture rejected")
            require(case["rejected_before_scientific_output"] is False, f"{case_id}: valid fixture marked rejected")
            require(case["reconstruction_exact"] is True, f"{case_id}: reconstruction did not match")
            require(type(load_json(stdout)) is dict, f"{case_id}: valid fixture stdout is not JSON")
    return observed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--aie", type=Path, required=True)
    parser.add_argument("--code-repo", type=Path)
    parser.add_argument("--scripts-repo", type=Path)
    parser.add_argument("--gate-a-result", type=Path)
    parser.add_argument("--panel", type=Path)
    parser.add_argument("--panel-metadata", type=Path)
    parser.add_argument("--archive", action="append", default=[])
    parser.add_argument("--time-binary", type=Path, default=Path("/usr/bin/time"))
    parser.add_argument("--threads", type=int, default=THREADS)
    parser.add_argument("--warm-blocks", type=int, default=WARM_BLOCKS)
    parser.add_argument("--cli-adapter", type=Path)
    parser.add_argument("--interface-status", choices=("draft", "final"), default="draft")
    parser.add_argument("--fixture-hook", action="append", type=Path, default=[])
    parser.add_argument("--allow-dirty-code", action="store_true")
    parser.add_argument("--allow-dirty-scripts", action="store_true")
    args = parser.parse_args()

    self_test_blake3()
    require(args.threads == THREADS, "the frozen Gate-B protocol requires exactly 24 threads")
    require(args.warm_blocks == WARM_BLOCKS, "the frozen Gate-B protocol requires eight paired blocks")
    root = args.project_root.resolve()
    run = args.run_dir.resolve()
    aie = args.aie.resolve()
    code_repo = (args.code_repo or root / "work/gravlax").resolve()
    scripts_repo = (args.scripts_repo or root / "gravlax-paper-scripts").resolve()
    gate_a_path = (args.gate_a_result or scripts_repo / "results/post-v1-archive-root-gate-a.json").resolve()
    panel_path = (args.panel or scripts_repo / "experiments/designs/archive-root-jshape-routing-panel.tsv").resolve()
    metadata_path = (args.panel_metadata or scripts_repo / "experiments/designs/archive-root-jshape-routing-panel.metadata.json").resolve()
    time_binary = args.time_binary.resolve()
    require(aie.is_file() and os.access(aie, os.X_OK), f"missing executable {aie}")
    require(time_binary.is_file() and os.access(time_binary, os.X_OK), f"missing time binary {time_binary}")
    require((code_repo / ".git").exists() and (scripts_repo / ".git").exists(), "missing Git checkout")
    require(not run.exists(), f"refusing to overwrite {run}")

    gate_a = load_json(gate_a_path)
    require(sha256_file(gate_a_path) == GATE_A_RESULT_SHA256, "Gate-A result SHA-256 differs")
    require(gate_a.get("schema") == "gravlax.archive-root-gate-a.v1" and gate_a.get("status") == "PASS", "Gate A did not pass")
    require(type(gate_a.get("gates")) is dict and all(gate_a["gates"].values()), "Gate-A result has a failed gate")
    require(gate_a_path.is_relative_to(scripts_repo), "Gate-A result is outside the scripts checkout")
    result_commit = git_output(
        scripts_repo,
        "log",
        "-n",
        "1",
        "--format=%H",
        "--",
        gate_a_path.relative_to(scripts_repo).as_posix(),
    )
    require(result_commit == GATE_A_RESULT_COMMIT, "Gate-A result-bearing commit differs")
    panel_cases, panel_metadata = load_frozen_panel(panel_path, metadata_path)
    require(sha256_file(metadata_path) == PANEL_METADATA_SHA256, "panel metadata identity differs")

    supplied: dict[str, Path] = {}
    for raw in args.archive:
        try:
            sample, path = parse_sample(raw)
        except (GateError, ValueError) as error:
            raise GateError(str(error)) from error
        require(sample not in supplied, f"duplicate archive {sample}")
        supplied[sample] = path
    recorded_run = Path(gate_a["identity"]["run_dir"]).resolve()
    archives = {
        sample: (supplied.get(sample) or recorded_run / "archives" / f"{sample}.v2.aie").resolve()
        for sample in SAMPLES
    }
    for sample, path in archives.items():
        expected = ROOTED_ARCHIVES[sample]
        require(path.is_file() and path.stat().st_size == expected["bytes"], f"{sample}: rooted archive size differs")
        require(sha256_file(path) == expected["sha256"], f"{sample}: rooted archive SHA-256 differs")
        parsed = parse_archive(path, ROOT_DOMAIN, verify_payloads=False)
        require(parsed.version == 2 and parsed.root == expected["archive_root"], f"{sample}: archive root differs")

    adapter = validate_adapter(DEFAULT_ADAPTER if args.cli_adapter is None else load_json(args.cli_adapter.resolve()))
    fixture_hooks = [path.resolve() for path in args.fixture_hook]
    required_fixture_hooks = [(scripts_repo / "scripts/214_jshape_routing_adversarial.py").resolve()]
    for hook in fixture_hooks:
        require(hook.is_relative_to(scripts_repo), f"fixture hook is outside the scripts checkout: {hook}")
        relative = hook.relative_to(scripts_repo).as_posix()
        git_output(scripts_repo, "ls-files", "--error-unmatch", "--", relative)
    code_commit = git_output(code_repo, "rev-parse", "HEAD")
    scripts_commit = git_output(scripts_repo, "rev-parse", "HEAD")
    code_dirty = bool(git_output(code_repo, "status", "--porcelain"))
    scripts_dirty = bool(git_output(scripts_repo, "status", "--porcelain"))
    if code_dirty:
        require(args.allow_dirty_code, "Gravlax worktree is dirty")
    if scripts_dirty:
        require(args.allow_dirty_scripts, "reproducibility worktree is dirty")
    promotable = is_promotable_run(
        args.interface_status,
        adapter=adapter,
        code_dirty=code_dirty,
        scripts_dirty=scripts_dirty,
        fixture_hooks=fixture_hooks,
        required_fixture_hooks=required_fixture_hooks,
    )
    if args.interface_status == "final":
        require(
            fixture_hooks == required_fixture_hooks,
            "the frozen final run requires exactly scripts/214_jshape_routing_adversarial.py",
        )
    require(
        adapter_is_frozen_final(adapter),
        "Gate-B CLI/JSON adapter is not the source-controlled FINAL contract",
    )

    prepare_output_tree(run)
    environment = {"RAYON_NUM_THREADS": str(args.threads)}
    binary_sha256 = sha256_file(aie)
    protocol = {
        "schema": "gravlax.jshape-route-gate-b-protocol.v1",
        "date": DATE,
        "interface_status": args.interface_status,
        "promotable_run": promotable,
        "threads": args.threads,
        "warm_paired_blocks": args.warm_blocks,
        "schedule": "odd blocks fallback-route; even blocks route-fallback",
        "panel_percentile": "nearest-rank p95 = sorted[ceil(0.95*n)-1]",
        "project_root": str(root),
        "run_dir": str(run),
        "aie": {"path": str(aie), "bytes": aie.stat().st_size, "sha256": binary_sha256},
        "gravlax": {"repo": str(code_repo), "commit": code_commit, "dirty": code_dirty},
        "paper_scripts": {"repo": str(scripts_repo), "commit": scripts_commit, "dirty": scripts_dirty},
        "adapter": adapter,
        "gate_a_result": {"path": str(gate_a_path), "commit": GATE_A_RESULT_COMMIT, "sha256": GATE_A_RESULT_SHA256},
        "panel": {"path": str(panel_path), "sha256": PANEL_SHA256, "metadata": str(metadata_path), "metadata_sha256": PANEL_METADATA_SHA256, "rows": len(panel_cases)},
        "archives": [
            {"id": sample, "path": str(archives[sample]), **ROOTED_ARCHIVES[sample]}
            for sample in SAMPLES
        ],
        "fixture_hooks": [str(path) for path in fixture_hooks],
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
            "python_zstandard": zstandard.__version__,
            "independent_blake3": {
                "implementation": "scripts/archive_root_gate_common.py",
                "implementation_sha256": sha256_file(scripts_repo / "scripts/archive_root_gate_common.py"),
                "known_vectors": "PASS",
            },
        },
    }
    write_json(run / "protocol.json", protocol)
    write_json(run / "audits/panel-metadata.json", panel_metadata)

    collections: dict[tuple[str, str], Path] = {}
    for arm in ("fallback", "route"):
        for form, sample_ids, base_form in (
            ("root", SAMPLES, None),
            ("repeat", SAMPLES, None),
            ("reverse", tuple(reversed(SAMPLES)), None),
            ("base", SAMPLES[:4], None),
            ("extension", SAMPLES[4:], "base"),
        ):
            output = run / "collections/canonical" / f"{arm}-{form}.aicollection"
            base = collections.get((arm, base_form)) if base_form else None
            before = sha256_file(base) if base is not None else None
            command = build_command(aie, archives, sample_ids, output, arm, adapter, base=base)
            command_dir = run / "commands/build-canonical" / f"{arm}-{form}"
            run_timed(command_dir, command, time_binary=time_binary, cwd=root, environment=environment)
            require_clean_success(command_dir)
            require(output.is_file(), f"{arm}-{form}: collection output missing")
            if base is not None:
                require(sha256_file(base) == before, f"{arm}: base changed during extension")
            collections[(arm, form)] = output
        root_hash = sha256_file(collections[(arm, "root")])
        require(sha256_file(collections[(arm, "repeat")]) == root_hash, f"{arm}: repeated build differs")
        require(sha256_file(collections[(arm, "reverse")]) == root_hash, f"{arm}: reversed-input build differs")
    write_json(
        run / "audits/collection-identities.json",
        {
            "schema": "gravlax.jshape-route-collection-identities.v1",
            "collections": {
                f"{arm}-{form}": {
                    "path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)
                }
                for (arm, form), path in sorted(collections.items())
            },
        },
    )

    for arm, form, mode in (
        ("fallback", "root", "fallback"),
        ("fallback", "extension", "fallback"),
        ("route", "root", "route"),
        ("route", "root", "route_full"),
        ("route", "extension", "route"),
        ("route", "extension", "route_full"),
    ):
        label = f"{arm}-{form}-{mode}"
        command_dir = run / "commands/inspect" / label
        run_timed(
            command_dir,
            inspect_command(aie, collections[(arm, form)], mode, adapter),
            time_binary=time_binary,
            cwd=root,
            environment=environment,
        )
        require_clean_success(command_dir)

    cases = fixed_cases() + panel_cases
    case_by_id = {case.case_id: case for case in cases}
    require(len(case_by_id) == len(cases) == 101, "Gate-B scientific case set differs")
    projections: dict[str, dict[str, Any]] = {}
    for case in cases:
        arms = {}
        for form in ("root", "extension"):
            for arm in ("fallback", "route"):
                label = f"{case.case_id}-{arm}-{form}"
                command_dir = run / "commands/query-canonical" / label
                run_timed(
                    command_dir,
                    query_command(aie, collections[(arm, form)], case, arm, adapter),
                    time_binary=time_binary,
                    cwd=root,
                    environment=environment,
                )
                raw = require_clean_success(command_dir)
                arms[f"{arm}-{form}"] = scientific_projection(raw)
        reference = arms["fallback-root"]
        require(all(value == reference for value in arms.values()), f"{case.case_id}: scientific result differs across arms/forms")
        if case.case_id in EXPECTED_TOTALS:
            require(reference.get("totals") == EXPECTED_TOTALS[case.case_id], f"{case.case_id}: frozen total differs")
        projections[case.case_id] = {
            "kind": case.kind,
            "locus": case.locus,
            "stratum": case.stratum,
            "scientific": reference,
        }
    write_json(
        run / "audits/scientific-equivalence.json",
        {
            "schema": "gravlax.jshape-route-scientific-equivalence.v1",
            "arms": ["fallback-root", "route-root", "fallback-extension", "route-extension"],
            "cases": projections,
        },
    )

    canonical_build_rows = build_schedule()
    for row in canonical_build_rows:
        arm, label = row["arm"], row["label"]
        output = run / "collections/benchmark" / f"{label}.aicollection"
        command_dir = run / "commands/build-benchmark" / label
        run_timed(
            command_dir,
            build_command(aie, archives, SAMPLES, output, arm, adapter),
            time_binary=time_binary,
            cwd=root,
            environment=environment,
        )
        require_clean_success(command_dir)
    write_build_schedule(run / "build-schedule.tsv", canonical_build_rows)

    fixed_timed_ids = ("dense", "sparse", "jset", "region")
    fixed_rows = build_schedule(fixed_timed_ids)
    for row in fixed_rows:
        case, arm, label = case_by_id[row["case_id"]], row["arm"], row["label"]
        command_dir = run / "commands/query-fixed-benchmark" / label
        run_timed(
            command_dir,
            query_command(aie, collections[(arm, "root")], case, arm, adapter),
            time_binary=time_binary,
            cwd=root,
            environment=environment,
        )
        value = require_clean_success(command_dir)
        require(scientific_projection(value) == projections[case.case_id]["scientific"], f"{label}: timed fixed result differs")
    write_query_schedule(run / "fixed-query-schedule.tsv", fixed_rows, case_by_id)

    panel_rows = build_schedule(case.case_id for case in panel_cases)
    for row in panel_rows:
        case, arm, label = case_by_id[row["case_id"]], row["arm"], row["label"]
        command_dir = run / "commands/query-panel-benchmark" / label
        run_timed(
            command_dir,
            query_command(aie, collections[(arm, "root")], case, arm, adapter),
            time_binary=time_binary,
            cwd=root,
            environment=environment,
        )
        value = require_clean_success(command_dir)
        require(scientific_projection(value) == projections[case.case_id]["scientific"], f"{label}: timed panel result differs")
    write_query_schedule(run / "panel-query-schedule.tsv", panel_rows, case_by_id)

    observed: dict[str, str] = {}
    hook_names = set()
    route_root = collections[("route", "root")]
    for hook in fixture_hooks:
        require(hook.is_file() and os.access(hook, os.X_OK), f"fixture hook unavailable: {hook}")
        name = re.sub(r"[^A-Za-z0-9_.-]+", "-", hook.stem)
        require(name and name not in hook_names, "duplicate fixture hook name")
        hook_names.add(name)
        output = run / "fixtures" / name
        command = [
            str(hook), "--aie", str(aie), "--out-dir", str(output),
            "--protocol", str(run / "protocol.json"),
            "--fallback-collection", str(collections[("fallback", "root")]),
            "--route-collection", str(route_root),
            "--route-chain", str(collections[("route", "extension")]),
            "--panel", str(panel_path),
        ]
        for sample in SAMPLES:
            command.extend(["--archive", f"{sample}={archives[sample]}"])
        command_dir = run / "commands/fixture-hooks" / name
        run_timed(command_dir, command, time_binary=time_binary, cwd=root, environment=environment)
        require_clean_success(command_dir, json_stdout=False)
        hook_cases = validate_fixture_result(output / "result.json", binary_sha256, sha256_file(route_root))
        require(not (set(observed) & set(hook_cases)), "fixture case ID repeated across hooks")
        observed.update(hook_cases)
    if args.interface_status == "final":
        require(observed == REQUIRED_FIXTURE_CASES, f"fixture coverage/expectations differ: {observed}")

    write_json(
        run / "driver-completion.json",
        {
            "schema": "gravlax.jshape-route-driver-completion.v1",
            "status": "COMPLETE" if promotable else "DEVELOPMENT_ONLY",
            "promotable_run": promotable,
            "panel_sha256": PANEL_SHA256,
            "fixture_cases": [
                {"id": case_id, "expect": observed[case_id]} for case_id in sorted(observed)
            ],
        },
    )
    manifest = write_artifact_manifest(run)
    print(json.dumps({"run_dir": str(run), "artifact_manifest": manifest}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except GateError as error:
        raise SystemExit(f"Gate-B benchmark failed: {error}") from error
