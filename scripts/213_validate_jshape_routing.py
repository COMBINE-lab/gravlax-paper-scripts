#!/usr/bin/env python3
"""Fail-closed reducer for the frozen Gate-B local-shape-routing campaign.

This reducer never runs Gravlax. It independently freezes the complete artifact set, checks every
recorded command and scientific answer, reconstructs the accounting and paired summaries, enforces
all frozen integrity/size/I/O/runtime/RSS gates, and writes compact PASS plus an excluded-artifact
inventory only after every check succeeds.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from archive_root_gate_common import (
    GateError,
    load_json,
    median,
    parse_archive,
    parse_gnu_time,
    require,
    self_test_blake3,
    sha256_file,
    validate_artifact_manifest,
)
from jshape_routing_gate_common import (
    AGGREGATE_SOURCE_BYTES,
    ALLOWED_ROUTE_BUILD_SECTIONS,
    COMPRESSED_SHAPES_AUDIT_BYTES,
    DEFAULT_ADAPTER,
    EXTERNAL_FIXTURE_CASES,
    EXPECTED_TOTALS,
    FORBIDDEN_BUILD_SECTION_PATTERNS,
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
    THRESHOLDS,
    WARM_BLOCKS,
    adapter_is_frozen_final,
    build_command,
    exact_keys,
    fixed_cases,
    inspect_command,
    load_frozen_panel,
    nearest_rank,
    nonnegative_integer,
    nonnegative_number,
    path_value,
    query_command,
    scientific_projection,
    validate_adapter,
)


COMMAND_FILES = {"command.json", "stdout.txt", "stderr.txt", "time.txt"}


def ratio(treatment: float, control: float, label: str) -> float:
    require(control > 0 and treatment >= 0, f"{label}: invalid ratio operands")
    return treatment / control


def logical_path(path: Path, project_root: Path) -> str:
    resolved = path.resolve()
    require(resolved.is_relative_to(project_root), f"artifact outside project root: {resolved}")
    relative = resolved.relative_to(project_root).as_posix()
    require(relative and "\t" not in relative and "\n" not in relative, "unsafe artifact path")
    return relative


def write_pass_outputs(
    *,
    result_path: Path,
    result: Mapping[str, Any],
    inventory_path: Path,
    inventory: str,
) -> None:
    """Atomically expose neither output, or both outputs with PASS written last."""
    result_path.parent.mkdir(parents=True, exist_ok=True)
    inventory_path.parent.mkdir(parents=True, exist_ok=True)
    result_temporary = result_path.with_name(result_path.name + ".tmp")
    inventory_temporary = inventory_path.with_name(inventory_path.name + ".tmp")
    require(not result_temporary.exists() and not inventory_temporary.exists(), "stale reducer temporary output")
    try:
        result_temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        inventory_temporary.write_text(inventory)
        os.replace(inventory_temporary, inventory_path)
        try:
            os.replace(result_temporary, result_path)
        except OSError:
            inventory_path.unlink(missing_ok=True)
            raise
    finally:
        result_temporary.unlink(missing_ok=True)
        inventory_temporary.unlink(missing_ok=True)


def command_files(
    directory: Path,
    *,
    expected_argv: Sequence[str] | None = None,
    expected_cwd: Path | None = None,
    json_stdout: bool = True,
) -> dict[str, Any]:
    require(directory.is_dir(), f"missing command directory {directory}")
    require({path.name for path in directory.iterdir()} == COMMAND_FILES, f"{directory}: command artifacts differ")
    command = load_json(directory / "command.json")
    exact_keys(command, {"schema", "argv", "cwd", "environment"}, str(directory / "command.json"))
    require(command["schema"] == "gravlax.archive-root-command.v1", f"{directory}: command schema differs")
    require(type(command["argv"]) is list and all(type(token) is str for token in command["argv"]), f"{directory}: argv malformed")
    require(command["environment"] == {"RAYON_NUM_THREADS": str(THREADS)}, f"{directory}: environment differs")
    if expected_argv is not None:
        require(command["argv"] == list(map(str, expected_argv)), f"{directory}: argv differs")
    if expected_cwd is not None:
        require(Path(command["cwd"]).resolve() == expected_cwd.resolve(), f"{directory}: cwd differs")
    resources = parse_gnu_time(directory / "time.txt")
    require((directory / "stderr.txt").stat().st_size == 0, f"{directory}: successful command wrote stderr")
    if json_stdout:
        stdout = load_json(directory / "stdout.txt")
        require(type(stdout) is dict, f"{directory}: stdout is not one JSON object")
    else:
        require((directory / "stdout.txt").stat().st_size == 0, f"{directory}: hook wrote stdout")
        stdout = None
    return {"command": command, "resources": resources, "stdout": stdout}


def parse_schedule(
    path: Path,
    cases: Mapping[str, Any],
    *,
    query: bool,
) -> list[dict[str, Any]]:
    lines = path.read_text().splitlines()
    expected_header = (
        "case_id\tkind\tlocus\tstratum\tblock\tposition\tarm\tlabel"
        if query
        else "case_id\tblock\tposition\tarm\tlabel"
    )
    require(lines and lines[0] == expected_header, f"{path}: schedule header differs")
    rows = []
    for raw in lines[1:]:
        fields = raw.split("\t")
        require(len(fields) == (8 if query else 5), f"{path}: malformed schedule row")
        if query:
            case_id, kind, locus, stratum, raw_block, raw_position, arm, label = fields
            require(case_id in cases, f"{path}: unknown query case")
            case = cases[case_id]
            require((kind, locus, stratum) == (case.kind, case.locus, case.stratum or ""), f"{path}: frozen case differs")
        else:
            case_id, raw_block, raw_position, arm, label = fields
            require(case_id == "build", f"{path}: build case differs")
        require(re.fullmatch(r"[1-8]", raw_block) is not None and raw_position in {"1", "2"}, f"{path}: block/position differs")
        block, position = int(raw_block), int(raw_position)
        expected_order = ("fallback", "route") if block % 2 else ("route", "fallback")
        require(arm == expected_order[position - 1], f"{path}: alternating order differs")
        require(label == f"{case_id}-b{block:02d}-p{position}-{arm}", f"{path}: label differs")
        rows.append({"case_id": case_id, "block": block, "position": position, "arm": arm, "label": label})
    expected_rows = len(cases) * 16 if query else 16
    require(len(rows) == expected_rows and len({row["label"] for row in rows}) == expected_rows, f"{path}: schedule cardinality differs")
    for case_id in cases:
        require(sum(row["case_id"] == case_id for row in rows) == 16, f"{path}: {case_id} row count differs")
    return rows


BUILD_FIELDS = {
    "schema", "output", "collection_format_version", "shape_routes_requested",
    "new_archives", "segment_junctions", "archive_routes", "chunk_postings",
    "raw_section_bytes", "file_bytes", "shape_routes", "elapsed_seconds", "source_io",
}
BUILD_ROUTE_FIELDS = {"archives", "sections", "exact_spans", "pairs", "raw_bytes", "compressed_bytes"}
BUILD_SOURCE_FIELDS = {
    "identity_content_bytes_read", "total_bytes_read", "shape_route_payload_bytes_read",
    "shape_route_source_bytes_read", "sections_read", "archives",
}
BUILD_ARCHIVE_SOURCE_FIELDS = {
    "id", "format_version", "identity_scheme", "identity_content_bytes_read",
    "total_bytes_read", "shape_route_payload_bytes_read", "shape_route_source_bytes_read",
    "sections_read",
}
QUERY_PHASE_KEYS = (
    "collection_load_seconds",
    "identity_guard_seconds",
    "route_planning_seconds",
    "source_execution_seconds",
    "total_seconds",
)


def validate_query_phases(
    value: Mapping[str, Any],
    paths: Mapping[str, str],
    label: str,
) -> dict[str, float]:
    result = {
        key: nonnegative_number(path_value(value, paths[key], f"{label}: {key}"), f"{label}.{key}")
        for key in QUERY_PHASE_KEYS
    }
    components = sum(result[key] for key in QUERY_PHASE_KEYS if key != "total_seconds")
    # The end-to-end timer also includes final aggregation and JSON construction. It must cover all
    # four disjoint measured phases, but is not required to equal their sum.
    require(result["total_seconds"] + 5e-6 >= components, f"{label}: phase timings exceed total")
    result["unattributed_seconds"] = max(0.0, result["total_seconds"] - components)
    return result


def validate_build_json(
    value: Mapping[str, Any],
    *,
    arm: str,
    output: Path,
    adapter: Mapping[str, Any],
    expected_samples: Sequence[str],
    expected_new_samples: Sequence[str] | None = None,
    expected_shape_payload_bytes: int | None = None,
) -> dict[str, Any]:
    require(type(value) is dict and value.get("schema") == adapter["schemas"]["build"], "build JSON schema differs")
    exact_keys(value, BUILD_FIELDS, "build JSON")
    require(Path(value.get("output", "")).resolve() == output.resolve(), "build JSON output path differs")
    require(value["collection_format_version"] == 4, "collection format version differs")
    expected_ids = sorted(expected_samples)
    expected_new_ids = sorted(expected_new_samples if expected_new_samples is not None else expected_samples)
    require(set(expected_new_ids).issubset(expected_ids), "new build archive IDs are not in the reported collection")
    require(value["new_archives"] == len(expected_new_ids), "build archive count differs")
    for key in ("segment_junctions", "archive_routes", "chunk_postings", "raw_section_bytes", "file_bytes"):
        nonnegative_integer(value[key], f"build {key}")
    elapsed_seconds = nonnegative_number(value["elapsed_seconds"], "build elapsed seconds")
    require(value["file_bytes"] == output.stat().st_size, "build file byte accounting differs")

    paths = adapter["build_paths"]
    enabled = path_value(value, paths["route_enabled"], "build route enabled")
    require(type(enabled) is bool and enabled is (arm == "route"), "build route arm differs")
    route = value["shape_routes"]
    exact_keys(route, BUILD_ROUTE_FIELDS, "build shape routes")
    route_archives = nonnegative_integer(path_value(value, paths["route_archives"], "build route archives"), "route archives")
    route_sections = nonnegative_integer(path_value(value, paths["route_sections"], "build route sections"), "route sections")
    exact_spans = nonnegative_integer(path_value(value, paths["route_exact_spans"], "build route spans"), "route exact spans")
    pairs = nonnegative_integer(path_value(value, paths["route_pairs"], "build route pairs"), "route pairs")
    raw_bytes = nonnegative_integer(path_value(value, paths["route_raw_bytes"], "build route raw bytes"), "route raw bytes")
    compressed_bytes = nonnegative_integer(path_value(value, paths["route_compressed_bytes"], "build route compressed bytes"), "route compressed bytes")

    source = value["source_io"]
    exact_keys(source, BUILD_SOURCE_FIELDS, "build source I/O")
    source_bytes = nonnegative_integer(path_value(value, paths["source_total_bytes"], "build source bytes"), "build source bytes")
    source_route_payload = nonnegative_integer(
        path_value(value, paths["source_route_payload_bytes"], "build route payload source bytes"),
        "build route payload source bytes",
    )
    source_route_source = nonnegative_integer(
        path_value(value, paths["source_route_source_bytes"], "build attributable route source bytes"),
        "build attributable route source bytes",
    )
    identity_content = nonnegative_integer(source["identity_content_bytes_read"], "build identity-content bytes")
    source_sections = path_value(value, paths["source_sections"], "build source sections")
    require(type(source_sections) is list and all(type(name) is str for name in source_sections), "build source sections malformed")
    require(source_sections == sorted(set(source_sections)), "build source sections unsorted or duplicated")
    for name in source_sections:
        require(name in ALLOWED_ROUTE_BUILD_SECTIONS, f"build read undeclared source section {name!r}")
        require(not any(pattern.search(name) for pattern in FORBIDDEN_BUILD_SECTION_PATTERNS), f"build read forbidden source section {name!r}")
    rows = source["archives"]
    require(type(rows) is list and [row.get("id") for row in rows if type(row) is dict] == expected_ids, "build source archive order/IDs differ")
    expected_sections = sorted(ALLOWED_ROUTE_BUILD_SECTIONS if arm == "route" else ALLOWED_ROUTE_BUILD_SECTIONS - {"shapes"})
    for expected_id, row in zip(expected_ids, rows):
        is_new = expected_id in expected_new_ids
        require(type(row) is dict, f"build source row {expected_id} malformed")
        exact_keys(row, BUILD_ARCHIVE_SOURCE_FIELDS, f"build source row {expected_id}")
        require(row["id"] == expected_id and row["format_version"] == 2, f"build source row {expected_id}: identity differs")
        require(row["identity_scheme"] == "aie-directory-root-v2", f"build source row {expected_id}: identity scheme differs")
        require(row["identity_content_bytes_read"] == 0, f"build source row {expected_id}: rooted identity scanned payload")
        require(
            row["sections_read"] == (expected_sections if is_new else []),
            f"build source row {expected_id}: sections differ",
        )
        row_total = nonnegative_integer(row["total_bytes_read"], f"build source row {expected_id}: total bytes")
        row_payload = nonnegative_integer(row["shape_route_payload_bytes_read"], f"build source row {expected_id}: route payload bytes")
        row_source = nonnegative_integer(row["shape_route_source_bytes_read"], f"build source row {expected_id}: route source bytes")
        if arm == "route" and is_new:
            require(row_payload > 0 and row_source >= row_payload and row_total >= row_source, f"build source row {expected_id}: route attribution differs")
        else:
            require(row_payload == 0 and row_source == 0, f"build source row {expected_id}: unexpected route attribution")
    require(identity_content == sum(row["identity_content_bytes_read"] for row in rows) == 0, "build aggregate identity bytes differ")
    require(source_bytes == sum(row["total_bytes_read"] for row in rows), "build aggregate source bytes differ")
    require(source_route_payload == sum(row["shape_route_payload_bytes_read"] for row in rows), "build aggregate route payload bytes differ")
    require(source_route_source == sum(row["shape_route_source_bytes_read"] for row in rows), "build aggregate route source bytes differ")
    require(source_sections == sorted({name for row in rows for name in row["sections_read"]}), "build aggregate source sections differ")

    if arm == "route":
        require(route_archives == len(expected_new_ids), "route archive count differs")
        require(route_sections > 0 and exact_spans > 0 and pairs > 0, "route build emitted no route rows")
        require(raw_bytes > 0 and compressed_bytes > 0, "route build emitted no route bytes")
        require(raw_bytes <= value["raw_section_bytes"] and compressed_bytes <= value["file_bytes"], "route bytes exceed collection totals")
        require("shapes" in source_sections, "route build did not read shapes")
        require(source_route_payload > 0 and source_route_source >= source_route_payload, "route build source attribution differs")
        require(source_bytes >= source_route_source, "route-attributable source bytes exceed total source bytes")
        if expected_shape_payload_bytes is not None:
            require(source_route_payload == expected_shape_payload_bytes, "route build shapes payload differs from frozen audit")
    else:
        require(
            (route_archives, route_sections, exact_spans, pairs, raw_bytes, compressed_bytes) == (0, 0, 0, 0, 0, 0),
            "fallback build emitted route state",
        )
        require("shapes" not in source_sections, "fallback build read shapes")
        require(source_route_payload == 0 and source_route_source == 0, "fallback build attributed route source bytes")
    return {
        "file_bytes": value["file_bytes"],
        "route_enabled": enabled,
        "route_archives": route_archives,
        "route_sections": route_sections,
        "route_exact_spans": exact_spans,
        "route_pairs": pairs,
        "route_raw_bytes": raw_bytes,
        "route_compressed_bytes": compressed_bytes,
        "source_total_bytes": source_bytes,
        "source_route_payload_bytes": source_route_payload,
        "source_route_source_bytes": source_route_source,
        "source_sections": source_sections,
        "elapsed_seconds": elapsed_seconds,
    }


def validate_query_json(
    value: Mapping[str, Any],
    *,
    arm: str,
    case: Any,
    adapter: Mapping[str, Any],
) -> dict[str, Any]:
    require(type(value) is dict and value.get("schema") == adapter["schemas"][case.kind], f"{case.case_id}: query schema differs")
    paths = adapter["query_paths"]
    expected_route_use = arm == "route" and case.kind in {"junction", "jset"}
    identity_bytes = nonnegative_integer(path_value(value, paths["source_identity_bytes"], f"{case.case_id}: source identity bytes"), "query source identity bytes")
    execution_bytes = nonnegative_integer(path_value(value, paths["source_execution_bytes"], f"{case.case_id}: source execution bytes"), "query source execution bytes")
    source_bytes = nonnegative_integer(path_value(value, paths["source_bytes"], f"{case.case_id}: source bytes"), "query source bytes")
    sidecar_bytes = nonnegative_integer(path_value(value, paths["sidecar_bytes"], f"{case.case_id}: sidecar bytes"), "query sidecar bytes")
    route_sidecar_bytes = nonnegative_integer(path_value(value, paths["route_sidecar_bytes"], f"{case.case_id}: route sidecar bytes"), "query route sidecar bytes")
    total_bytes = nonnegative_integer(path_value(value, paths["total_bytes"], f"{case.case_id}: total bytes"), "query total bytes")
    route_blocks = nonnegative_integer(path_value(value, paths["route_blocks_loaded"], f"{case.case_id}: route blocks"), "query route blocks")
    routed_archives = nonnegative_integer(path_value(value, paths["routed_archives"], f"{case.case_id}: routed archives"), "query routed archives")
    fallback_archives = nonnegative_integer(path_value(value, paths["fallback_archives"], f"{case.case_id}: fallback archives"), "query fallback archives")
    archives_opened = nonnegative_integer(path_value(value, paths["archives_opened"], f"{case.case_id}: archives opened"), "query archives opened")
    require(source_bytes == identity_bytes + execution_bytes, f"{case.case_id}: source byte components do not sum")
    require(total_bytes == source_bytes + sidecar_bytes, f"{case.case_id}: total logical bytes do not sum")
    require(route_sidecar_bytes <= sidecar_bytes, f"{case.case_id}: route payload exceeds collection-sidecar bytes")
    if case.kind == "region":
        require((route_blocks, routed_archives, fallback_archives, route_sidecar_bytes) == (0, 0, 0, 0), f"{case.case_id}: region used junction routes")
    elif expected_route_use and archives_opened > 0:
        require(routed_archives + fallback_archives == archives_opened, f"{case.case_id}: routed/fallback archive accounting differs")
        if case.kind == "junction":
            require(routed_archives == archives_opened and fallback_archives == 0, f"{case.case_id}: point junction unexpectedly fell back")
        else:
            require(routed_archives > 0, f"{case.case_id}: routed junction-set used no route")
        require(route_blocks >= routed_archives and route_sidecar_bytes > 0, f"{case.case_id}: routed block accounting differs")
    else:
        require(routed_archives == 0 and route_blocks == 0 and route_sidecar_bytes == 0, f"{case.case_id}: non-routed query loaded route state")
        expected_fallback = archives_opened if case.kind in {"junction", "jset"} else 0
        require(fallback_archives == expected_fallback, f"{case.case_id}: fallback archive accounting differs")
    phases = validate_query_phases(value, paths, f"{case.case_id}: phases")
    return {
        "source_identity_bytes": identity_bytes,
        "source_execution_bytes": execution_bytes,
        "source_bytes": source_bytes,
        "sidecar_bytes": sidecar_bytes,
        "route_sidecar_bytes": route_sidecar_bytes,
        "total_bytes": total_bytes,
        "route_blocks_loaded": route_blocks,
        "routed_archives": routed_archives,
        "fallback_archives": fallback_archives,
        "archives_opened": archives_opened,
        "phase_seconds": phases,
    }


def stable_query_accounting(value: Mapping[str, Any]) -> dict[str, Any]:
    return {key: item for key, item in value.items() if key != "phase_seconds"}


def require_absent_source_payload_zero(
    accounting: Mapping[tuple[str, str, str], Mapping[str, Any]],
) -> None:
    """Check source execution I/O, while allowing rooted identity-directory reads."""
    for arm in ("fallback", "route"):
        for form in ("root", "extension"):
            key = ("absent", arm, form)
            require(key in accounting, f"absent {arm} {form}: missing canonical accounting")
            require(
                accounting[key]["source_execution_bytes"] == 0,
                f"absent {arm} {form}: read source execution payload bytes",
            )


def canonical_build_reported_samples(form: str) -> tuple[str, ...]:
    """Return the canonical IDs reported after sorting inputs and flattening a base chain."""
    if form == "base":
        return SAMPLES[:4]
    if form in {"root", "repeat", "reverse", "extension"}:
        return SAMPLES
    raise GateError(f"unknown canonical build form {form}")


INSPECT_FIELDS = {
    "schema", "format_version", "path", "file_bytes", "layers", "reference",
    "chromosomes", "archives", "index", "guard", "io", "elapsed_seconds",
}
INSPECT_LAYER_FIELDS = {
    "path", "format_version", "root_digest", "archives", "junction_rows",
    "junction_segments", "shape_routes",
}
INSPECT_ROUTE_SUMMARY_FIELDS = {
    "archives", "sections", "exact_spans", "pairs", "compressed_bytes",
    "reconstruction_verified",
}
INSPECT_ARCHIVE_FIELDS = {
    "id", "path", "bytes", "archive_format_version", "native_identity",
    "encoded_sections_identity", "chunks", "shape_routes",
}
INSPECT_BINDING_FIELDS = {
    "codec_version", "archive_ordinal", "source_root", "shapes_digest", "n_shapes", "sections",
    "exact_spans", "pairs", "compressed_bytes", "reconstruction_verified", "descriptors",
}
INSPECT_DESCRIPTOR_FIELDS = {
    "first_span", "last_span", "section_name", "exact_spans", "pairs", "compressed_bytes",
}
INSPECT_INDEX_FIELDS = {
    "segment_junction_rows", "global_junctions", "archive_routes", "chunk_postings",
    "interval_chunks", "junction_bin_bp", "junction_segments", "shape_route_archives",
    "shape_route_sections", "shape_route_exact_spans", "shape_route_pairs",
    "shape_route_compressed_bytes",
}
INSPECT_GUARD_FIELDS = {
    "filesystem_identity", "content_digest_recorded", "content_digest_verified",
    "shape_route_payloads_verified", "shape_route_reconstruction_verified",
}
INSPECT_IO_FIELDS = {
    "collection_sidecar_bytes_read", "source_identity_content_bytes_read",
    "source_identity_total_bytes_read", "route_verification_source_bytes_read",
}


def validate_collection_inspect(
    value: Mapping[str, Any],
    *,
    arm: str,
    full: bool,
    adapter: Mapping[str, Any],
    collection: Path,
    layer_paths: Sequence[Path],
    source_directories: Mapping[str, Any],
) -> dict[str, Any]:
    require(type(value) is dict and value.get("schema") == adapter["schemas"]["inspect"], "inspect schema differs")
    exact_keys(value, INSPECT_FIELDS, "collection inspect")
    require(value["format_version"] == 4, "inspect collection format differs")
    require(Path(value["path"]).resolve() == collection.resolve(), "inspect collection path differs")
    require(value["file_bytes"] == sum(path.stat().st_size for path in layer_paths), "inspect collection byte count differs")
    nonnegative_number(value["elapsed_seconds"], "inspect elapsed seconds")
    paths = adapter["inspect_paths"]
    route_archives = nonnegative_integer(path_value(value, paths["route_archives"], "inspect route archives"), "inspect route archives")
    route_sections = nonnegative_integer(path_value(value, paths["route_sections"], "inspect route sections"), "inspect route sections")
    exact_spans = nonnegative_integer(path_value(value, paths["route_exact_spans"], "inspect route exact spans"), "inspect route exact spans")
    pairs = nonnegative_integer(path_value(value, paths["route_pairs"], "inspect route pairs"), "inspect route pairs")
    compressed = nonnegative_integer(path_value(value, paths["route_compressed_bytes"], "inspect route bytes"), "inspect route bytes")
    expected_route = arm == "route"
    require((route_archives > 0) is expected_route, "inspect route arm differs")

    layers = value["layers"]
    require(type(layers) is list and len(layers) == len(layer_paths), "inspect layer count differs")
    layer_route_summaries = []
    expected_local_ordinals: list[int] = []
    for layer_index, (layer, expected_path) in enumerate(zip(layers, layer_paths)):
        require(type(layer) is dict, f"inspect layer {layer_index} malformed")
        exact_keys(layer, INSPECT_LAYER_FIELDS, f"inspect layer {layer_index}")
        require(Path(layer["path"]).resolve() == expected_path.resolve(), f"inspect layer {layer_index}: path differs")
        require(layer["format_version"] == 4, f"inspect layer {layer_index}: format differs")
        require(re.fullmatch(r"[0-9a-f]{64}", str(layer["root_digest"])) is not None, f"inspect layer {layer_index}: root malformed")
        layer_archives = nonnegative_integer(layer["archives"], f"inspect layer {layer_index}: archives")
        expected_local_ordinals.extend(range(layer_archives))
        nonnegative_integer(layer["junction_rows"], f"inspect layer {layer_index}: junction rows")
        nonnegative_integer(layer["junction_segments"], f"inspect layer {layer_index}: junction segments")
        route_summary = layer["shape_routes"]
        require(type(route_summary) is dict, f"inspect layer {layer_index}: route summary malformed")
        exact_keys(route_summary, INSPECT_ROUTE_SUMMARY_FIELDS, f"inspect layer {layer_index}: route summary")
        numeric_summary = {
            key: nonnegative_integer(route_summary[key], f"inspect layer {layer_index}: {key}")
            for key in ("archives", "sections", "exact_spans", "pairs", "compressed_bytes")
        }
        require(route_summary["reconstruction_verified"] is full, f"inspect layer {layer_index}: reconstruction flag differs")
        if expected_route:
            require(numeric_summary["archives"] == layer_archives and all(numeric_summary.values()), f"inspect layer {layer_index}: route cardinalities differ")
        else:
            require(all(value == 0 for value in numeric_summary.values()), f"inspect layer {layer_index}: fallback has route state")
        layer_route_summaries.append({**numeric_summary, "reconstruction_verified": full})

    archive_rows = value["archives"]
    require(type(archive_rows) is list and len(archive_rows) == len(SAMPLES), "inspect archive count differs")
    require(len(expected_local_ordinals) == len(SAMPLES), "inspect layer archive totals differ")
    binding_rows = []
    for global_ordinal, (sample, local_ordinal, archive) in enumerate(zip(SAMPLES, expected_local_ordinals, archive_rows)):
        require(type(archive) is dict, f"inspect archive {sample} malformed")
        exact_keys(archive, INSPECT_ARCHIVE_FIELDS, f"inspect archive {sample}")
        source = source_directories[sample]
        require(archive["id"] == sample and Path(archive["path"]).resolve() == source.path.resolve(), f"inspect archive {sample}: path/ID differs")
        require(archive["bytes"] == source.file_bytes and archive["archive_format_version"] == 2, f"inspect archive {sample}: source size/format differs")
        exact_keys(archive["native_identity"], {"scheme", "blake3"}, f"inspect archive {sample}: native identity")
        require(
            archive["native_identity"] == {"scheme": "aie-directory-root-v2", "blake3": ROOTED_ARCHIVES[sample]["archive_root"]},
            f"inspect archive {sample}: native identity differs",
        )
        exact_keys(archive["encoded_sections_identity"], {"scheme", "blake3"}, f"inspect archive {sample}: encoded identity")
        require(archive["encoded_sections_identity"]["scheme"] == "aie-encoded-sections-v1", f"inspect archive {sample}: encoded scheme differs")
        require(re.fullmatch(r"[0-9a-f]{64}", str(archive["encoded_sections_identity"]["blake3"])) is not None, f"inspect archive {sample}: encoded digest malformed")
        nonnegative_integer(archive["chunks"], f"inspect archive {sample}: chunks")
        binding = archive["shape_routes"]
        if not expected_route:
            require(binding is None, f"inspect archive {sample}: fallback has route binding")
            continue
        require(type(binding) is dict, f"inspect archive {sample}: route binding missing")
        exact_keys(binding, INSPECT_BINDING_FIELDS, f"inspect archive {sample}: route binding")
        require(binding["codec_version"] == 1, f"inspect archive {sample}: route codec differs")
        require(binding["archive_ordinal"] == local_ordinal, f"inspect archive {sample}: bound archive ordinal differs")
        require(binding["source_root"] == ROOTED_ARCHIVES[sample]["archive_root"], f"inspect archive {sample}: bound root differs")
        shapes_entry = next((section for section in source.sections if section.name == "shapes"), None)
        require(shapes_entry is not None and binding["shapes_digest"] == shapes_entry.committed_blake3, f"inspect archive {sample}: bound shapes digest differs")
        require(type(binding["n_shapes"]) is int and binding["n_shapes"] > 0, f"inspect archive {sample}: n_shapes differs")
        require(binding["reconstruction_verified"] is full, f"inspect archive {sample}: reconstruction flag differs")
        descriptors = binding["descriptors"]
        require(type(descriptors) is list and descriptors, f"inspect archive {sample}: descriptors missing")
        previous_last = None
        descriptor_totals = {"sections": 0, "exact_spans": 0, "pairs": 0, "compressed_bytes": 0}
        for descriptor in descriptors:
            exact_keys(descriptor, INSPECT_DESCRIPTOR_FIELDS, f"inspect archive {sample}: descriptor")
            first = nonnegative_integer(descriptor["first_span"], f"inspect archive {sample}: first span")
            last = nonnegative_integer(descriptor["last_span"], f"inspect archive {sample}: last span")
            spans = nonnegative_integer(descriptor["exact_spans"], f"inspect archive {sample}: exact spans")
            descriptor_pairs = nonnegative_integer(descriptor["pairs"], f"inspect archive {sample}: pairs")
            descriptor_bytes = nonnegative_integer(descriptor["compressed_bytes"], f"inspect archive {sample}: compressed bytes")
            require(first <= last <= 0xFFFFFFFF and 1 <= spans <= 256, f"inspect archive {sample}: descriptor span bounds differ")
            require(descriptor_pairs >= spans and descriptor_bytes > 0, f"inspect archive {sample}: descriptor cardinalities differ")
            require(descriptor["section_name"] == f"s.{local_ordinal}.{first}", f"inspect archive {sample}: section name differs")
            require(previous_last is None or first > previous_last, f"inspect archive {sample}: descriptors overlap/unsorted")
            previous_last = last
            descriptor_totals["sections"] += 1
            descriptor_totals["exact_spans"] += spans
            descriptor_totals["pairs"] += descriptor_pairs
            descriptor_totals["compressed_bytes"] += descriptor_bytes
        for key, total in descriptor_totals.items():
            require(binding[key] == total, f"inspect archive {sample}: binding {key} differs")
        binding_rows.append(binding)

    aggregate = {
        "archives": route_archives,
        "sections": route_sections,
        "exact_spans": exact_spans,
        "pairs": pairs,
        "compressed_bytes": compressed,
    }
    for key in aggregate:
        require(aggregate[key] == sum(row[key] for row in layer_route_summaries), f"inspect aggregate {key} differs from layers")
    if expected_route:
        require(route_archives == len(SAMPLES) and all(value > 0 for value in aggregate.values()), "inspect routed aggregate differs")
        require(route_archives == len(binding_rows), "inspect aggregate archives differ from bindings")
        for key in aggregate.keys() - {"archives"}:
            require(aggregate[key] == sum(row[key] for row in binding_rows), f"inspect aggregate {key} differs from bindings")
    else:
        require(all(value == 0 for value in aggregate.values()) and not binding_rows, "inspect fallback aggregate differs")

    index = value["index"]
    require(type(index) is dict, "inspect index malformed")
    exact_keys(index, INSPECT_INDEX_FIELDS, "inspect index")
    require(index["global_junctions"] is None, "inspect global-junction sentinel differs")
    for key in INSPECT_INDEX_FIELDS - {"global_junctions"}:
        nonnegative_integer(index[key], f"inspect index {key}")
    guard = value["guard"]
    require(type(guard) is dict, "inspect guard malformed")
    exact_keys(guard, INSPECT_GUARD_FIELDS, "inspect guard")
    require(guard["content_digest_recorded"] is True and guard["content_digest_verified"] is False, "inspect source content guard differs")
    require(path_value(value, paths["payloads_verified"], "inspect route payload verification") is True, "inspect did not verify route payloads")
    require(path_value(value, paths["reconstruction_verified"], "inspect reconstruction") is full, "inspect reconstruction flag differs")
    io = value["io"]
    require(type(io) is dict, "inspect I/O malformed")
    exact_keys(io, INSPECT_IO_FIELDS, "inspect I/O")
    sidecar_bytes = nonnegative_integer(path_value(value, paths["sidecar_bytes"], "inspect sidecar bytes"), "inspect sidecar bytes")
    identity_content = nonnegative_integer(path_value(value, paths["identity_content_bytes"], "inspect identity content bytes"), "inspect identity content bytes")
    identity_total = nonnegative_integer(path_value(value, paths["identity_total_bytes"], "inspect identity total bytes"), "inspect identity total bytes")
    verification_source = nonnegative_integer(path_value(value, paths["verification_source_bytes"], "inspect verification source bytes"), "inspect verification source bytes")
    require(sidecar_bytes > 0 and identity_content == 0 and identity_total > 0, "inspect base I/O differs")
    require((verification_source > 0) is (full and expected_route), "inspect route-verification source I/O differs")
    binding_digest = hashlib.sha256(
        json.dumps(
            [{key: item for key, item in binding.items() if key != "reconstruction_verified"} for binding in binding_rows],
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()
    return {
        **aggregate,
        "binding_sha256": binding_digest,
        "layers": layer_route_summaries,
        "sidecar_bytes": sidecar_bytes,
        "identity_content_bytes": identity_content,
        "identity_total_bytes": identity_total,
        "verification_source_bytes": verification_source,
        "reconstruction_verified": full,
    }


def paired_summary(rows: Sequence[Mapping[str, Any]], records: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    by_arm = {"fallback": [], "route": []}
    for row in rows:
        by_arm[row["arm"]].append(records[row["label"]]["resources"])
    result = {}
    for field in ("wall_seconds", "user_seconds", "system_seconds", "max_rss_kib"):
        fallback = [record[field] for record in by_arm["fallback"]]
        route = [record[field] for record in by_arm["route"]]
        fallback_median, route_median = median(fallback), median(route)
        result[field] = {
            "fallback": fallback,
            "route": route,
            "fallback_median": fallback_median,
            "route_median": route_median,
            "route_over_fallback": ratio(route_median, fallback_median, field),
        }
    return result


def validate_fixtures(
    run: Path,
    *,
    protocol: Mapping[str, Any],
    binary_sha256: str,
    route_sha256: str,
) -> dict[str, Any]:
    observed: dict[str, str] = {}
    hooks = []
    fixture_root = run / "fixtures"
    hook_dirs = sorted(path for path in fixture_root.iterdir() if path.is_dir())
    require(hook_dirs, "no Gate-B fixture output")
    expected_hook_paths = [Path(path).resolve() for path in protocol["fixture_hooks"]]
    expected_by_name = {re.sub(r"[^A-Za-z0-9_.-]+", "-", path.stem): path for path in expected_hook_paths}
    require(set(expected_by_name) == {path.name for path in hook_dirs}, "fixture hook directories differ")
    archives = {row["id"]: Path(row["path"]).resolve() for row in protocol["archives"]}
    for hook_dir in hook_dirs:
        hook = expected_by_name[hook_dir.name]
        expected_command = [
            str(hook), "--aie", protocol["aie"]["path"], "--out-dir", str(hook_dir),
            "--protocol", str(run / "protocol.json"),
            "--fallback-collection", str(run / "collections/canonical/fallback-root.aicollection"),
            "--route-collection", str(run / "collections/canonical/route-root.aicollection"),
            "--route-chain", str(run / "collections/canonical/route-extension.aicollection"),
            "--panel", protocol["panel"]["path"],
        ]
        for sample in SAMPLES:
            expected_command.extend(["--archive", f"{sample}={archives[sample]}"])
        command_files(
            run / "commands/fixture-hooks" / hook_dir.name,
            expected_argv=expected_command,
            expected_cwd=Path(protocol["project_root"]),
            json_stdout=False,
        )
        result_path = hook_dir / "result.json"
        value = load_json(result_path)
        exact_keys(
            value,
            {
                "schema", "binary_sha256", "route_collection_sha256", "panel_sha256",
                "code_commit", "tools", "rust_tests", "cases",
            },
            "fixture result",
        )
        require(value["schema"] == "gravlax.jshape-route-adversarial.v2", "fixture schema differs")
        require(value["binary_sha256"] == binary_sha256 and value["route_collection_sha256"] == route_sha256, "fixture identities differ")
        require(value["panel_sha256"] == PANEL_SHA256, "fixture panel differs")
        require(value["code_commit"] == protocol["gravlax"]["commit"], "fixture code commit differs")
        exact_keys(value["tools"], {"python", "zstandard_python"}, "fixture tools")
        require(all(type(item) is str and item for item in value["tools"].values()), "fixture tool versions differ")
        rust_tests = value["rust_tests"]
        require(type(rust_tests) is list and [row.get("test_name") for row in rust_tests if type(row) is dict] == list(RUST_FIXTURE_TESTS), "fixture Rust-test set/order differs")
        code_repo = Path(protocol["gravlax"]["repo"]).resolve()
        rust_by_name = {}
        for row in rust_tests:
            exact_keys(
                row,
                {"test_name", "command", "exit_status", "stdout", "stderr", "source", "source_bytes", "source_sha256", "cases"},
                "fixture Rust test",
            )
            test_name = row["test_name"]
            expected_command = ["cargo", "test", "-p", "aie", "--bin", "aie", test_name, "--", "--exact", "--test-threads=1"]
            require(row["command"] == expected_command and row["exit_status"] == 0, f"fixture Rust test {test_name}: command/status differs")
            require(tuple(row["cases"]) == RUST_FIXTURE_TESTS[test_name], f"fixture Rust test {test_name}: case mapping differs")
            source = code_repo / "crates/aie/src" / ("shaperoute.rs" if test_name.startswith("shaperoute::") else "collectioncmd.rs")
            require(row["source"] == str(source) and row["source_bytes"] == source.stat().st_size and row["source_sha256"] == sha256_file(source), f"fixture Rust test {test_name}: source differs")
            for stream in ("stdout", "stderr"):
                stream_path = hook_dir / row[stream]
                require(type(row[stream]) is str and not Path(row[stream]).is_absolute() and stream_path.is_file(), f"fixture Rust test {test_name}: stream differs")
            require(re.search(r"test result: ok\. 1 passed;", (hook_dir / row["stdout"]).read_text()) is not None, f"fixture Rust test {test_name}: exact test did not pass")
            rust_by_name[test_name] = row
        cases = value["cases"]
        require(type(cases) is list and cases, "fixture cases empty")
        for case in cases:
            exact_keys(
                case,
                {
                    "id", "expect", "evidence_kind", "test_name", "command", "exit_status", "stdout", "stderr", "fixture_bytes",
                    "fixture_sha256", "mutation", "rejected_before_scientific_output",
                    "reconstruction_exact",
                },
                "fixture case",
            )
            case_id, expect = case["id"], case["expect"]
            require(case_id not in observed and expect in {"accept", "reject"}, "duplicate fixture case")
            observed[case_id] = expect
            stdout, stderr = hook_dir / case["stdout"], hook_dir / case["stderr"]
            require(stdout.is_file() and stderr.is_file(), f"{case_id}: fixture streams missing")
            require(re.fullmatch(r"[0-9a-f]{64}", str(case["fixture_sha256"])) is not None, f"{case_id}: fixture digest malformed")
            require(case["expect"] == REQUIRED_FIXTURE_CASES.get(case_id), f"{case_id}: fixture expectation differs")
            if case["evidence_kind"] == "rust_unit_test":
                test_name = case["test_name"]
                require(case_id not in EXTERNAL_FIXTURE_CASES and test_name in rust_by_name and case_id in RUST_FIXTURE_TESTS[test_name], f"{case_id}: Rust evidence mapping differs")
                require(case["command"] == rust_by_name[test_name]["command"] and case["exit_status"] == 0, f"{case_id}: Rust evidence command differs")
                require(case["fixture_bytes"] == rust_by_name[test_name]["source_bytes"] and case["fixture_sha256"] == rust_by_name[test_name]["source_sha256"], f"{case_id}: Rust evidence identity differs")
                require(case["rejected_before_scientific_output"] is (expect == "reject"), f"{case_id}: Rust rejection predicate differs")
                require(case["reconstruction_exact"] is (expect == "accept"), f"{case_id}: Rust exactness predicate differs")
            elif expect == "reject":
                require(case["evidence_kind"] == "cli_mutation" and case_id in EXTERNAL_FIXTURE_CASES and case["test_name"] is None, f"{case_id}: external evidence mapping differs")
                require(case["exit_status"] != 0 and case["rejected_before_scientific_output"] is True, f"{case_id}: invalid route accepted/partial")
                require(stdout.stat().st_size == 0 and case["reconstruction_exact"] is False, f"{case_id}: rejection contract differs")
            else:
                raise GateError(f"{case_id}: external acceptance evidence is not allowed")
        hooks.append({"name": hook_dir.name, "cases": len(cases), "result_sha256": sha256_file(result_path)})
    require(observed == REQUIRED_FIXTURE_CASES, f"fixture coverage/expectations differ: {observed}")
    return {"cases": len(observed), "hooks": hooks, "all_required_cases_behaved_as_expected": True}


def inventory_text(
    *,
    project_root: Path,
    run: Path,
    artifact_manifest: Mapping[str, Any],
    binary: Path,
    panel: Path,
) -> str:
    rows: list[tuple[str, str, int, str, str, str]] = []
    role = "generated by scripts/212_benchmark_jshape_routing.py and validated by scripts/213_validate_jshape_routing.py; do not commit the artifact"

    def add(path: Path, kind: str, size: int, scheme: str, digest: str, description: str = role) -> None:
        require(type(size) is int and size >= 0 and re.fullmatch(r"[0-9a-f]{64}", digest) is not None, "inventory identity malformed")
        rows.append((logical_path(path, project_root), kind, size, scheme, digest, description))

    add(run, "gate_b_run_tree", artifact_manifest["bytes"], "complete-artifact-manifest-sha256", artifact_manifest["manifest_sha256"])
    add(binary, "executable", binary.stat().st_size, "sha256", sha256_file(binary))
    add(panel, "frozen_query_panel", panel.stat().st_size, "sha256", sha256_file(panel), "committed prospective Gate-B panel; do not regenerate after treatment")
    source_by_id = {row["id"]: Path(row["path"]) for row in load_json(run / "protocol.json")["archives"]}
    for sample in SAMPLES:
        source = source_by_id[sample]
        add(source, "rooted_aie", source.stat().st_size, "sha256", sha256_file(source))
    for path in sorted((run / "collections").rglob("*.aicollection"), key=lambda item: item.relative_to(run).as_posix()):
        add(path, "collection", path.stat().st_size, "sha256", sha256_file(path))
    require(len({(row[0], row[1]) for row in rows}) == len(rows), "duplicate inventory row")
    rows.sort(key=lambda row: (row[0], row[1]))
    lines = ["logical_path\tkind\tbytes\tdigest_scheme\tdigest\trole_and_regeneration"]
    lines.extend("\t".join(map(str, row)) for row in rows)
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--artifact-inventory-out", type=Path, required=True)
    parser.add_argument("--gravlax-commit", required=True)
    parser.add_argument("--paper-scripts-commit", required=True)
    args = parser.parse_args()

    self_test_blake3()
    run, out, inventory_out = args.run_dir.resolve(), args.out.resolve(), args.artifact_inventory_out.resolve()
    require(run.is_dir(), f"missing run directory {run}")
    require(not out.exists() and not inventory_out.exists() and out != inventory_out, "refusing to overwrite reducer output")
    require(re.fullmatch(r"[0-9a-f]{40}", args.gravlax_commit) is not None, "invalid Gravlax commit")
    require(re.fullmatch(r"[0-9a-f]{40}", args.paper_scripts_commit) is not None, "invalid scripts commit")

    artifact_manifest = validate_artifact_manifest(run)
    protocol = load_json(run / "protocol.json")
    exact_keys(
        protocol,
        {
            "schema", "date", "interface_status", "promotable_run", "threads",
            "warm_paired_blocks", "schedule", "panel_percentile", "project_root", "run_dir",
            "aie", "gravlax", "paper_scripts", "adapter", "gate_a_result", "panel",
            "archives", "fixture_hooks", "host", "tools",
        },
        "Gate-B protocol",
    )
    require(protocol["schema"] == "gravlax.jshape-route-gate-b-protocol.v1", "protocol schema differs")
    require(protocol["interface_status"] == "final" and protocol["promotable_run"] is True, "run is not promotable")
    require(protocol["threads"] == THREADS and protocol["warm_paired_blocks"] == WARM_BLOCKS, "thread/block protocol differs")
    require(protocol["schedule"] == "odd blocks fallback-route; even blocks route-fallback", "schedule contract differs")
    require(protocol["panel_percentile"] == "nearest-rank p95 = sorted[ceil(0.95*n)-1]", "percentile contract differs")
    adapter = validate_adapter(protocol["adapter"])
    require(adapter_is_frozen_final(adapter), "source-controlled Gate-B adapter is not FINAL")
    require(protocol["gravlax"] == {"repo": protocol["gravlax"]["repo"], "commit": args.gravlax_commit, "dirty": False}, "Gravlax identity differs")
    require(protocol["paper_scripts"] == {"repo": protocol["paper_scripts"]["repo"], "commit": args.paper_scripts_commit, "dirty": False}, "scripts identity differs")
    expected_fixture_hook = (
        Path(protocol["paper_scripts"]["repo"]).resolve()
        / "scripts/214_jshape_routing_adversarial.py"
    ).resolve()
    require(
        protocol["fixture_hooks"] == [str(expected_fixture_hook)],
        "fixture hook is not the frozen scripts/214_jshape_routing_adversarial.py",
    )
    require(protocol["gate_a_result"] == {"path": protocol["gate_a_result"]["path"], "commit": GATE_A_RESULT_COMMIT, "sha256": GATE_A_RESULT_SHA256}, "Gate-A identity differs")
    gate_a_path = Path(protocol["gate_a_result"]["path"])
    require(gate_a_path.is_file() and sha256_file(gate_a_path) == GATE_A_RESULT_SHA256, "Gate-A result file differs")
    binary = Path(protocol["aie"]["path"]).resolve()
    binary_sha256 = sha256_file(binary)
    require(binary.is_file() and binary.stat().st_size == protocol["aie"]["bytes"] and binary_sha256 == protocol["aie"]["sha256"], "binary identity differs")
    project_root = Path(protocol["project_root"]).resolve()
    require(Path(protocol["run_dir"]).resolve() == run and run.is_relative_to(project_root), "run path differs")

    panel = Path(protocol["panel"]["path"]).resolve()
    panel_metadata = Path(protocol["panel"]["metadata"]).resolve()
    panel_cases, metadata = load_frozen_panel(panel, panel_metadata)
    require(protocol["panel"] == {"path": str(panel), "sha256": PANEL_SHA256, "metadata": str(panel_metadata), "metadata_sha256": PANEL_METADATA_SHA256, "rows": 96}, "protocol panel differs")
    require(load_json(run / "audits/panel-metadata.json") == metadata, "copied panel metadata differs")
    require(
        [(row["id"], row["bytes"], row["sha256"], row["archive_root"]) for row in protocol["archives"]]
        == [(sample, ROOTED_ARCHIVES[sample]["bytes"], ROOTED_ARCHIVES[sample]["sha256"], ROOTED_ARCHIVES[sample]["archive_root"]) for sample in SAMPLES],
        "rooted archive protocol identities differ",
    )
    archives = {row["id"]: Path(row["path"]).resolve() for row in protocol["archives"]}
    source_directories = {}
    for sample, path in archives.items():
        expected = ROOTED_ARCHIVES[sample]
        require(path.is_file() and path.stat().st_size == expected["bytes"] and sha256_file(path) == expected["sha256"], f"{sample}: rooted archive differs")
        parsed = parse_archive(path, ROOT_DOMAIN, verify_payloads=False)
        require(parsed.version == 2 and parsed.root == expected["archive_root"], f"{sample}: root differs")
        source_directories[sample] = parsed

    completion = load_json(run / "driver-completion.json")
    exact_keys(completion, {"schema", "status", "promotable_run", "panel_sha256", "fixture_cases"}, "driver completion")
    require(completion["schema"] == "gravlax.jshape-route-driver-completion.v1", "completion schema differs")
    require(completion["status"] == "COMPLETE" and completion["promotable_run"] is True and completion["panel_sha256"] == PANEL_SHA256, "driver did not complete promotably")
    require(completion["fixture_cases"] == [{"id": case, "expect": REQUIRED_FIXTURE_CASES[case]} for case in sorted(REQUIRED_FIXTURE_CASES)], "completion fixture coverage differs")

    case_list = fixed_cases() + panel_cases
    cases = {case.case_id: case for case in case_list}
    require(len(cases) == 101, "scientific case set differs")
    identities = load_json(run / "audits/collection-identities.json")
    exact_keys(identities, {"schema", "collections"}, "collection identity audit")
    require(identities["schema"] == "gravlax.jshape-route-collection-identities.v1", "collection identity schema differs")
    collections: dict[tuple[str, str], Path] = {}
    build_summaries = {}
    for arm in ("fallback", "route"):
        for form, sample_ids, base_form in (
            ("root", SAMPLES, None),
            ("repeat", SAMPLES, None),
            ("reverse", tuple(reversed(SAMPLES)), None),
            ("base", SAMPLES[:4], None),
            ("extension", SAMPLES[4:], "base"),
        ):
            label = f"{arm}-{form}"
            output = run / "collections/canonical" / f"{label}.aicollection"
            base = collections.get((arm, base_form)) if base_form else None
            expected = build_command(binary, archives, sample_ids, output, arm, adapter, base=base)
            record = command_files(run / "commands/build-canonical" / label, expected_argv=expected, expected_cwd=project_root)
            expected_payload = COMPRESSED_SHAPES_AUDIT_BYTES if arm == "route" and form in {"root", "repeat", "reverse"} else None
            summary = validate_build_json(
                record["stdout"],
                arm=arm,
                output=output,
                adapter=adapter,
                expected_samples=canonical_build_reported_samples(form),
                expected_new_samples=sample_ids,
                expected_shape_payload_bytes=expected_payload,
            )
            audit_row = identities["collections"].get(label)
            require(audit_row == {"path": str(output), "bytes": output.stat().st_size, "sha256": sha256_file(output)}, f"{label}: identity audit differs")
            collections[(arm, form)] = output
            build_summaries[label] = summary
        root_digest = sha256_file(collections[(arm, "root")])
        require(sha256_file(collections[(arm, "repeat")]) == root_digest, f"{arm}: repeat differs")
        require(sha256_file(collections[(arm, "reverse")]) == root_digest, f"{arm}: reverse differs")
    require(
        build_summaries["route-base"]["source_route_payload_bytes"]
        + build_summaries["route-extension"]["source_route_payload_bytes"]
        == COMPRESSED_SHAPES_AUDIT_BYTES,
        "route chain shapes payload differs from frozen audit",
    )

    inspect_summaries = {}
    for arm, form, mode in (
        ("fallback", "root", "fallback"),
        ("fallback", "extension", "fallback"),
        ("route", "root", "route"),
        ("route", "root", "route_full"),
        ("route", "extension", "route"),
        ("route", "extension", "route_full"),
    ):
        label = f"{arm}-{form}-{mode}"
        expected = inspect_command(binary, collections[(arm, form)], mode, adapter)
        record = command_files(run / "commands/inspect" / label, expected_argv=expected, expected_cwd=project_root)
        layer_paths = (
            [collections[(arm, "root")]]
            if form == "root"
            else [collections[(arm, "base")], collections[(arm, "extension")]]
        )
        inspect_summaries[label] = validate_collection_inspect(
            record["stdout"],
            arm=arm,
            full=mode == "route_full",
            adapter=adapter,
            collection=collections[(arm, form)],
            layer_paths=layer_paths,
            source_directories=source_directories,
        )
    for form in ("root", "extension"):
        lazy, full = inspect_summaries[f"route-{form}-route"], inspect_summaries[f"route-{form}-route_full"]
        stable_keys = {
            "archives", "sections", "exact_spans", "pairs", "compressed_bytes",
            "binding_sha256", "sidecar_bytes", "identity_content_bytes", "identity_total_bytes",
        }
        require(
            {key: lazy[key] for key in stable_keys} == {key: full[key] for key in stable_keys},
            f"{form}: full inspection changed route metadata/accounting",
        )
        require(lazy["verification_source_bytes"] == 0 and full["verification_source_bytes"] > 0, f"{form}: reconstruction I/O differs")
    require(
        inspect_summaries["route-root-route_full"]["verification_source_bytes"]
        == build_summaries["route-root"]["source_route_source_bytes"],
        "root route reconstruction/source-build attribution differs",
    )
    require(
        inspect_summaries["route-extension-route_full"]["verification_source_bytes"]
        == build_summaries["route-base"]["source_route_source_bytes"]
        + build_summaries["route-extension"]["source_route_source_bytes"],
        "chain route reconstruction/source-build attribution differs",
    )

    scientific = load_json(run / "audits/scientific-equivalence.json")
    exact_keys(scientific, {"schema", "arms", "cases"}, "scientific equivalence")
    require(scientific["schema"] == "gravlax.jshape-route-scientific-equivalence.v1", "scientific audit schema differs")
    require(scientific["arms"] == ["fallback-root", "route-root", "fallback-extension", "route-extension"], "scientific arms differ")
    require(set(scientific["cases"]) == set(cases), "scientific case IDs differ")
    canonical_accounting = {}
    for case_id, case in cases.items():
        audit = scientific["cases"][case_id]
        require(audit == {"kind": case.kind, "locus": case.locus, "stratum": case.stratum, "scientific": audit["scientific"]}, f"{case_id}: audit fields differ")
        if case_id in EXPECTED_TOTALS:
            require(audit["scientific"].get("totals") == EXPECTED_TOTALS[case_id], f"{case_id}: frozen totals differ")
        for form in ("root", "extension"):
            for arm in ("fallback", "route"):
                label = f"{case_id}-{arm}-{form}"
                expected = query_command(binary, collections[(arm, form)], case, arm, adapter)
                record = command_files(run / "commands/query-canonical" / label, expected_argv=expected, expected_cwd=project_root)
                require(scientific_projection(record["stdout"]) == audit["scientific"], f"{label}: scientific JSON differs")
                accounting = validate_query_json(record["stdout"], arm=arm, case=case, adapter=adapter)
                canonical_accounting[(case_id, arm, form)] = accounting

    build_rows = parse_schedule(run / "build-schedule.tsv", {"build": object()}, query=False)
    build_records = {}
    for row in build_rows:
        output = run / "collections/benchmark" / f"{row['label']}.aicollection"
        expected = build_command(binary, archives, SAMPLES, output, row["arm"], adapter)
        record = command_files(run / "commands/build-benchmark" / row["label"], expected_argv=expected, expected_cwd=project_root)
        validate_build_json(
            record["stdout"],
            arm=row["arm"],
            output=output,
            adapter=adapter,
            expected_samples=SAMPLES,
            expected_shape_payload_bytes=COMPRESSED_SHAPES_AUDIT_BYTES if row["arm"] == "route" else None,
        )
        build_records[row["label"]] = record
    build_performance = paired_summary(build_rows, build_records)
    require(build_performance["wall_seconds"]["route_over_fallback"] <= THRESHOLDS["construction_wall_ratio_max"], "route construction wall gate failed")
    require(build_performance["max_rss_kib"]["route_median"] <= THRESHOLDS["construction_peak_rss_kib_max"], "route construction RSS gate failed")

    fixed_timed = {case.case_id: case for case in case_list if case.case_id in {"dense", "sparse", "jset", "region"}}
    fixed_rows = parse_schedule(run / "fixed-query-schedule.tsv", fixed_timed, query=True)
    fixed_records = {}
    fixed_accounting = {}
    for row in fixed_rows:
        case, arm, label = cases[row["case_id"]], row["arm"], row["label"]
        expected = query_command(binary, collections[(arm, "root")], case, arm, adapter)
        record = command_files(run / "commands/query-fixed-benchmark" / label, expected_argv=expected, expected_cwd=project_root)
        require(scientific_projection(record["stdout"]) == scientific["cases"][case.case_id]["scientific"], f"{label}: timed scientific result differs")
        accounting = validate_query_json(record["stdout"], arm=arm, case=case, adapter=adapter)
        require(
            stable_query_accounting(accounting)
            == stable_query_accounting(canonical_accounting[(case.case_id, arm, "root")]),
            f"{label}: timed accounting differs",
        )
        fixed_records[label], fixed_accounting[label] = record, accounting

    fixed_performance = {}
    for case_id in fixed_timed:
        rows = [row for row in fixed_rows if row["case_id"] == case_id]
        summary = paired_summary(rows, fixed_records)
        fallback_accounting = fixed_accounting[next(row["label"] for row in rows if row["arm"] == "fallback")]
        route_accounting = fixed_accounting[next(row["label"] for row in rows if row["arm"] == "route")]
        source_ratio = ratio(route_accounting["source_bytes"], fallback_accounting["source_bytes"], f"{case_id} source bytes")
        total_ratio = ratio(route_accounting["total_bytes"], fallback_accounting["total_bytes"], f"{case_id} total bytes")
        if case_id in {"dense", "sparse", "jset"}:
            require(source_ratio <= THRESHOLDS["fixed_query_source_ratio_max"], f"{case_id}: source-I/O gate failed")
        else:
            require(total_ratio <= THRESHOLDS["region_regression_ratio_max"], "region byte-regression gate failed")
            require(summary["wall_seconds"]["route_over_fallback"] <= THRESHOLDS["region_regression_ratio_max"], "region wall-regression gate failed")
        if case_id == "dense":
            require(summary["wall_seconds"]["route_over_fallback"] <= THRESHOLDS["dense_wall_ratio_max"], "dense wall gate failed")
        if case_id == "sparse":
            require(summary["wall_seconds"]["route_over_fallback"] <= THRESHOLDS["sparse_wall_ratio_max"], "sparse wall gate failed")
        if case_id == "jset":
            require(summary["wall_seconds"]["route_over_fallback"] <= THRESHOLDS["jset_wall_ratio_max"], "jset wall gate failed")
        if case_id in {"dense", "jset"}:
            require(summary["max_rss_kib"]["route_over_fallback"] <= THRESHOLDS["dense_jset_rss_ratio_max"], f"{case_id}: RSS ratio gate failed")
            require(summary["max_rss_kib"]["route_median"] <= THRESHOLDS["dense_jset_rss_kib_max"], f"{case_id}: absolute RSS gate failed")
        if case_id == "sparse":
            require(summary["max_rss_kib"]["route_over_fallback"] <= THRESHOLDS["sparse_rss_ratio_max"], "sparse RSS ratio gate failed")
            require(summary["max_rss_kib"]["route_median"] <= THRESHOLDS["sparse_rss_kib_max"], "sparse absolute RSS gate failed")
        fixed_performance[case_id] = {"resources": summary, "source_ratio": source_ratio, "total_ratio": total_ratio}
    require_absent_source_payload_zero(canonical_accounting)

    panel_case_map = {case.case_id: case for case in panel_cases}
    panel_rows = parse_schedule(run / "panel-query-schedule.tsv", panel_case_map, query=True)
    panel_records = {}
    panel_accounting = {}
    for row in panel_rows:
        case, arm, label = cases[row["case_id"]], row["arm"], row["label"]
        expected = query_command(binary, collections[(arm, "root")], case, arm, adapter)
        record = command_files(run / "commands/query-panel-benchmark" / label, expected_argv=expected, expected_cwd=project_root)
        require(scientific_projection(record["stdout"]) == scientific["cases"][case.case_id]["scientific"], f"{label}: panel scientific result differs")
        accounting = validate_query_json(record["stdout"], arm=arm, case=case, adapter=adapter)
        require(
            stable_query_accounting(accounting)
            == stable_query_accounting(canonical_accounting[(case.case_id, arm, "root")]),
            f"{label}: panel accounting differs",
        )
        panel_records[label], panel_accounting[label] = record, accounting

    panel_by_case = {}
    source_ratios, total_ratios = [], []
    for case_id in panel_case_map:
        rows = [row for row in panel_rows if row["case_id"] == case_id]
        fallback = panel_accounting[next(row["label"] for row in rows if row["arm"] == "fallback")]
        route = panel_accounting[next(row["label"] for row in rows if row["arm"] == "route")]
        source_ratio = ratio(route["source_bytes"], fallback["source_bytes"], f"{case_id}: source bytes")
        total_ratio = ratio(route["total_bytes"], fallback["total_bytes"], f"{case_id}: total bytes")
        require(source_ratio <= THRESHOLDS["panel_source_per_query_ratio_max"], f"{case_id}: per-query source-I/O gate failed")
        source_ratios.append(source_ratio)
        total_ratios.append(total_ratio)
        panel_by_case[case_id] = {
            "stratum": cases[case_id].stratum,
            "locus": cases[case_id].locus,
            "source_ratio": source_ratio,
            "total_ratio": total_ratio,
            "resources": paired_summary(rows, panel_records),
        }
    source_median, source_p95 = median(source_ratios), nearest_rank(source_ratios, 0.95)
    total_median, total_p95 = median(total_ratios), nearest_rank(total_ratios, 0.95)
    require(source_median <= THRESHOLDS["panel_source_median_ratio_max"] and source_p95 <= THRESHOLDS["panel_source_p95_ratio_max"], "panel source-I/O distribution gate failed")
    require(total_median <= THRESHOLDS["panel_total_median_ratio_max"] and total_p95 <= THRESHOLDS["panel_total_p95_ratio_max"], "panel total-I/O distribution gate failed")
    aggregate_wall = {"fallback": [], "route": []}
    for block in range(1, WARM_BLOCKS + 1):
        for arm in ("fallback", "route"):
            labels = [row["label"] for row in panel_rows if row["block"] == block and row["arm"] == arm]
            require(len(labels) == 96, "panel block cardinality differs")
            aggregate_wall[arm].append(sum(panel_records[label]["resources"]["wall_seconds"] for label in labels))
    aggregate_wall_ratio = ratio(median(aggregate_wall["route"]), median(aggregate_wall["fallback"]), "panel aggregate wall")
    require(aggregate_wall_ratio <= THRESHOLDS["panel_aggregate_wall_ratio_max"], "panel aggregate wall gate failed")

    fallback_root_bytes = collections[("fallback", "root")].stat().st_size
    route_root_bytes = collections[("route", "root")].stat().st_size
    fallback_chain_bytes = sum(collections[("fallback", form)].stat().st_size for form in ("base", "extension"))
    route_chain_bytes = sum(collections[("route", form)].stat().st_size for form in ("base", "extension"))
    require(route_root_bytes >= fallback_root_bytes and route_chain_bytes >= fallback_chain_bytes, "route collection is smaller than fallback unexpectedly")
    root_premium = route_root_bytes - fallback_root_bytes
    chain_premium = route_chain_bytes - fallback_chain_bytes
    require(root_premium / AGGREGATE_SOURCE_BYTES <= THRESHOLDS["route_premium_fraction_max"], "root route premium gate failed")
    require(chain_premium / AGGREGATE_SOURCE_BYTES <= THRESHOLDS["route_premium_fraction_max"], "chain route premium gate failed")
    require(route_root_bytes / AGGREGATE_SOURCE_BYTES <= THRESHOLDS["complete_collection_fraction_max"], "root complete collection size gate failed")
    require(route_chain_bytes / AGGREGATE_SOURCE_BYTES <= THRESHOLDS["complete_collection_fraction_max"], "chain complete collection size gate failed")
    route_compressed = build_summaries["route-root"]["route_compressed_bytes"]
    require(route_compressed / COMPRESSED_SHAPES_AUDIT_BYTES <= THRESHOLDS["route_over_compressed_shapes_ratio_max"], "compressed route/shapes gate failed")
    require(inspect_summaries["route-root-route"]["compressed_bytes"] == route_compressed, "route build/inspect byte accounting differs")

    route_root_sha = sha256_file(collections[("route", "root")])
    fixtures = validate_fixtures(run, protocol=protocol, binary_sha256=binary_sha256, route_sha256=route_root_sha)
    inventory = inventory_text(project_root=project_root, run=run, artifact_manifest=artifact_manifest, binary=binary, panel=panel)
    inventory_sha256 = hashlib.sha256(inventory.encode()).hexdigest()
    result = {
        "schema": "gravlax.jshape-route-gate-b.v1",
        "status": "PASS",
        "identity": {
            "gravlax_commit": args.gravlax_commit,
            "paper_scripts_commit": args.paper_scripts_commit,
            "run_dir": str(run),
            "artifact_manifest": artifact_manifest,
            "binary": {"path": str(binary), "bytes": binary.stat().st_size, "sha256": binary_sha256},
            "gate_a_result": {"commit": GATE_A_RESULT_COMMIT, "sha256": GATE_A_RESULT_SHA256},
            "panel": {"path": logical_path(panel, project_root), "sha256": PANEL_SHA256, "rows": 96},
            "excluded_artifact_inventory": {"path": logical_path(inventory_out, project_root), "bytes": len(inventory.encode()), "sha256": inventory_sha256},
        },
        "protocol": {
            "date": protocol["date"], "threads": protocol["threads"],
            "warm_paired_blocks": protocol["warm_paired_blocks"], "host": protocol["host"],
            "tools": protocol["tools"], "percentile": protocol["panel_percentile"],
        },
        "representation": {
            "format_version": 1,
            "bucket": "exact intron span",
            "record": ["shape_id", "donor_offset"],
            "acceptor_offset": "checked donor_offset + exact bucket span",
            "authoritative_counts": "source molecule chunks",
        },
        "collections": {
            "builds": build_summaries,
            "inspections": inspect_summaries,
            "size": {
                "aggregate_source_bytes": AGGREGATE_SOURCE_BYTES,
                "fallback_root_bytes": fallback_root_bytes,
                "route_root_bytes": route_root_bytes,
                "root_route_premium_bytes": root_premium,
                "root_route_premium_fraction": root_premium / AGGREGATE_SOURCE_BYTES,
                "fallback_chain_bytes": fallback_chain_bytes,
                "route_chain_bytes": route_chain_bytes,
                "chain_route_premium_bytes": chain_premium,
                "chain_route_premium_fraction": chain_premium / AGGREGATE_SOURCE_BYTES,
                "compressed_route_bytes": route_compressed,
                "compressed_shapes_audit_bytes": COMPRESSED_SHAPES_AUDIT_BYTES,
                "route_over_compressed_shapes": route_compressed / COMPRESSED_SHAPES_AUDIT_BYTES,
            },
            "construction_performance": build_performance,
        },
        "scientific_exactness": {
            "cases": len(cases), "panel_cases": len(panel_cases),
            "root_chain_fallback_route_exact": True, "fixed_totals": EXPECTED_TOTALS,
        },
        "fixed_queries": fixed_performance,
        "panel": {
            "cases": panel_by_case,
            "source_ratio_median": source_median,
            "source_ratio_p95_nearest_rank": source_p95,
            "total_ratio_median": total_median,
            "total_ratio_p95_nearest_rank": total_p95,
            "aggregate_wall_seconds": aggregate_wall,
            "aggregate_wall_route_over_fallback": aggregate_wall_ratio,
        },
        "adversarial": fixtures,
        "gates": {
            "all_101_scientific_cases_exact_across_root_chain_and_arms": True,
            "route_reconstruction_exact": True,
            "route_binding_and_negative_fixtures_pass": True,
            "repeated_and_reversed_builds_byte_identical": True,
            "size_gates_pass": True,
            "construction_gates_pass": True,
            "source_and_total_io_gates_pass": True,
            "runtime_and_memory_gates_pass": True,
            "absent_query_reads_zero_source_payload_bytes": True,
        },
        "claim_guards": [
            "Derived routes never supply molecular counts.",
            "Point-junction and junction-set results do not imply replay or region acceleration.",
            "The 33,897,246-byte compressed-shapes value was frozen preimplementation as an input denominator; earlier route-size estimates were projections.",
            "Performance has no biological-validity implication.",
        ],
    }
    write_pass_outputs(result_path=out, result=result, inventory_path=inventory_out, inventory=inventory)
    print(json.dumps({"status": "PASS", "out": str(out), "artifact_inventory": str(inventory_out), "artifact_inventory_sha256": inventory_sha256}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except GateError as error:
        raise SystemExit(f"Gate-B validation failed: {error}") from error
