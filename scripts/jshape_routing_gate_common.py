#!/usr/bin/env python3
"""Frozen constants and strict primitives shared by the Gate-B driver and reducer.

The command adapter freezes the final Gravlax route CLI and structured JSON contracts. Gate
thresholds, panel identity, scientific cases, schedule, and fixture coverage are not
adapter-configurable.
"""

from __future__ import annotations

import csv
import math
import re
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from archive_root_gate_common import GateError, blake3_bytes, load_json, require, sha256_file


DATE = "2026-09-01"
ROOT_DOMAIN = b"gravlax-aie-directory-root-v2\0"
SELECTION_DOMAIN = b"gravlax-jshape-panel-v1\0"
SAMPLES = tuple("ABCDEFGH")
GATE_A_RESULT_COMMIT = "9cba7af125dfce2d0b0c0b14ea56e92e00215a5d"
GATE_A_RESULT_SHA256 = "3f75fb1f7d5cc8b726405017a320c3698b7384c5da1f833f6d1f0a29453429db"
PANEL_SHA256 = "55980706cd6fbf0c07a1e652b330e5aca6946b1365010c8ef4e8cb6d6c4f1ca1"
PANEL_BYTES = 11_578
PANEL_ROWS = 96
PANEL_METADATA_SHA256 = "c3a0cd6d7cd53dd5e734f425f3dbc5c106911a0bda3ab4627dcee1dad2031d37"
AGGREGATE_SOURCE_BYTES = 1_326_036_691
COMPRESSED_SHAPES_AUDIT_BYTES = 33_897_246
THREADS = 24
WARM_BLOCKS = 8


ROOTED_ARCHIVES: dict[str, dict[str, Any]] = {
    "A": {"bytes": 110_224_112, "sha256": "66f9cec66c0cf08a02ecf53641aa709465b6d306ac6eb76dd7605802b326f463", "archive_root": "a2b361ba9f436c9c37d108e5b8b09f478854d129ecbb3802d4bad6aef38f14d7"},
    "B": {"bytes": 83_708_098, "sha256": "991e906cadcbab61e1223c01065803bae1265eb628b74b11fba953340f4b40d1", "archive_root": "1682bbdba7b5ca241aa3aa4c9eda8022944a88a7cad9d2e8b972c6f151133ed6"},
    "C": {"bytes": 68_210_860, "sha256": "8ab72a252cad7ac7311a7aab09cb84aa9e5c338b79ee482df4faf73659b26738", "archive_root": "acbda19bd8ff2f472c9cc2f87c2dac05b260a46690d49c1e0b08a5bd667e7a17"},
    "D": {"bytes": 397_913_953, "sha256": "a478b0168d22e797525058af238c706d552d02f45a3d63a7a0becadf8409fa95", "archive_root": "1b07a8383a7b8d583155a7fd1799369a4ae963aa9d890160329f86365eaa8095"},
    "E": {"bytes": 64_188_688, "sha256": "11e2a42817da24eacc18c1f39d7b1f5996f89fcd700c6969b3fff8a8fbc9965a", "archive_root": "e6729534e9c10e7d0cff51f83980433187b9cbe62d6bbcbdcda4d1a9624d168a"},
    "F": {"bytes": 91_983_279, "sha256": "bfa29a8dcd599e192b13413024b75fdfb276b655834adbd7c1a3cd9457c4208f", "archive_root": "7ac8737019a626fe712dcdb72f4f7c16d48aaefb5c935759163437652f739971"},
    "G": {"bytes": 436_454_139, "sha256": "827ae5ae24e2ccfaf4b8352cf73db18fb964930ebf7b83b7517326f32ce673b8", "archive_root": "05ee7c9f298a4da7dc089dd77bbc171fa30cb2ffe0e52bf536dca51ca30f8941"},
    "H": {"bytes": 73_709_402, "sha256": "8b6894b69fdb6a9a5d6a17151a7edff2c1bcb7d3f124b577ba0c4d6604cc86c9", "archive_root": "21c7b1088a1c07bcd7900a537526c5d1356c48d2ddab9e61fe2b98954bdc6521"},
}


STRATA = {
    "one_archive": (1, 1, 24, 1_139_256),
    "two_to_three_archives": (2, 3, 24, 114_543),
    "four_to_seven_archives": (4, 7, 24, 68_436),
    "eight_archives": (8, 8, 24, 22_035),
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


THRESHOLDS: dict[str, float | int] = {
    "route_premium_fraction_max": 0.03,
    "complete_collection_fraction_max": 0.04,
    "route_over_compressed_shapes_ratio_max": 1.25,
    "construction_wall_ratio_max": 2.0,
    "construction_peak_rss_kib_max": 768 * 1024,
    "fixed_query_source_ratio_max": 0.20,
    "panel_source_median_ratio_max": 0.25,
    "panel_source_p95_ratio_max": 0.50,
    "panel_source_per_query_ratio_max": 1.05,
    "panel_total_median_ratio_max": 0.50,
    "panel_total_p95_ratio_max": 0.75,
    "panel_aggregate_wall_ratio_max": 0.75,
    "dense_wall_ratio_max": 0.80,
    "sparse_wall_ratio_max": 0.85,
    "jset_wall_ratio_max": 0.80,
    "region_regression_ratio_max": 1.05,
    "dense_jset_rss_ratio_max": 0.65,
    "dense_jset_rss_kib_max": 256 * 1024,
    "sparse_rss_ratio_max": 0.75,
    "sparse_rss_kib_max": 128 * 1024,
}


ALLOWED_ROUTE_BUILD_SECTIONS = {
    "meta",
    "chroms",
    "rans.tables",
    "shapes",
    "index.chunks",
    "index.junctions",
    "index.jpost",
}
FORBIDDEN_BUILD_SECTION_PATTERNS = (
    re.compile(r"^c[0-9]+$"),
    re.compile(r"^coc(?:\.|$)"),
    re.compile(r"^patterns(?:\.|$)"),
    re.compile(r"^edges(?:\.|$)"),
)


REQUIRED_FIXTURE_CASES: dict[str, str] = {
    "single_block_shape": "accept",
    "multiple_introns": "accept",
    "repeated_equal_span_distinct_offsets": "accept",
    "distinct_shapes_same_span": "accept",
    "forward_record": "accept",
    "reverse_record": "accept",
    "both_representatives": "accept",
    "multimapper": "accept",
    "umi_class_crossing_chunks": "accept",
    "absent_coordinate": "accept",
    "fallback_route_absent": "accept",
    "unknown_route_version": "reject",
    "incorrect_source_root": "reject",
    "incorrect_shapes_digest": "reject",
    "out_of_range_shape_id": "reject",
    "reversed_offsets": "reject",
    "inconsistent_source_intron": "reject",
    "unsorted_pairs": "reject",
    "duplicate_pairs": "reject",
    "wrong_archive_ordinal": "reject",
    "wrong_span_bucket": "reject",
    "count_overflow": "reject",
    "checked_addition_overflow": "reject",
    "truncation": "reject",
    "checksum_mismatch": "reject",
    "trailing_bytes": "reject",
}


RUST_FIXTURE_TESTS: dict[str, tuple[str, ...]] = {
    "shaperoute::tests::streaming_derivation_preserves_repeated_spans_offsets_and_shape_ids": (
        "multiple_introns",
        "repeated_equal_span_distinct_offsets",
        "distinct_shapes_same_span",
    ),
    "shaperoute::tests::single_block_shape_has_no_route_and_distinct_spans_partition_at_256": (
        "single_block_shape",
    ),
    "shaperoute::tests::block_codec_is_canonical_bounded_and_context_bound": (
        "duplicate_pairs",
    ),
    "shaperoute::tests::malformed_pairs_shapes_and_checked_additions_are_rejected": (
        "reversed_offsets",
        "unsorted_pairs",
        "count_overflow",
        "checked_addition_overflow",
    ),
    "shaperoute::tests::binding_names_directory_and_full_reconstruction_are_exact": (),
    "collectioncmd::tests::routed_v4_manifest_sections_are_deterministic_and_fully_reconstructable": (),
    "collectioncmd::tests::route_binding_codec_rejects_wrong_context_names_intervals_and_legacy_sources": (),
    "collectioncmd::tests::route_free_v4_is_deterministic_and_selects_the_exact_fallback": (
        "fallback_route_absent",
    ),
    "collectioncmd::tests::routed_build_is_identical_for_repeated_and_reversed_sample_inputs": (),
    "collectioncmd::tests::routed_and_fallback_reducers_match_on_geometry_multimappers_and_cross_chunk_umis": (
        "forward_record",
        "reverse_record",
        "both_representatives",
        "multimapper",
        "umi_class_crossing_chunks",
        "absent_coordinate",
    ),
}


EXTERNAL_FIXTURE_CASES: dict[str, str] = {
    "unknown_route_version": "manifest_codec_version",
    "incorrect_source_root": "manifest_source_root",
    "incorrect_shapes_digest": "manifest_shapes_digest",
    "out_of_range_shape_id": "block_shape_id",
    "inconsistent_source_intron": "block_reconstruction_mismatch",
    "wrong_archive_ordinal": "block_archive_ordinal",
    "wrong_span_bucket": "block_span_bucket",
    "truncation": "container_truncation",
    "checksum_mismatch": "route_payload_corruption",
    "trailing_bytes": "container_trailing_bytes",
}


require(
    set(REQUIRED_FIXTURE_CASES)
    == set(EXTERNAL_FIXTURE_CASES)
    | {case for cases in RUST_FIXTURE_TESTS.values() for case in cases},
    "Gate-B fixture evidence map does not partition the required cases",
)


# Frozen CLI/JSON adapter. Arrays are executed directly; no token is interpreted by a shell.
DEFAULT_ADAPTER: dict[str, Any] = {
    "schema": "gravlax.jshape-route-cli-adapter.v1",
    "status": "FINAL",
    "build_fallback_extra": ["--json"],
    "build_route_extra": ["--shape-routes", "--json"],
    "query_fallback_extra": [],
    "query_route_extra": [],
    "inspect_fallback": ["{aie}", "collection", "inspect", "{collection}"],
    "inspect_route": ["{aie}", "collection", "inspect", "{collection}"],
    "inspect_route_full": ["{aie}", "collection", "inspect", "{collection}", "--verify-routes"],
    "schemas": {
        "build": "gravlax.collection.build.v2",
        "inspect": "gravlax.collection.v4",
        "junction": "gravlax.collection.junction.v2",
        "region": "gravlax.collection.region.v2",
        "jset": "gravlax.collection.jset.v2",
    },
    "build_paths": {
        "route_enabled": "shape_routes_requested",
        "route_archives": "shape_routes.archives",
        "route_sections": "shape_routes.sections",
        "route_exact_spans": "shape_routes.exact_spans",
        "route_pairs": "shape_routes.pairs",
        "route_raw_bytes": "shape_routes.raw_bytes",
        "route_compressed_bytes": "shape_routes.compressed_bytes",
        "source_total_bytes": "source_io.total_bytes_read",
        "source_route_payload_bytes": "source_io.shape_route_payload_bytes_read",
        "source_route_source_bytes": "source_io.shape_route_source_bytes_read",
        "source_sections": "source_io.sections_read",
        "elapsed_seconds": "elapsed_seconds",
    },
    "inspect_paths": {
        "route_archives": "index.shape_route_archives",
        "route_sections": "index.shape_route_sections",
        "route_exact_spans": "index.shape_route_exact_spans",
        "route_pairs": "index.shape_route_pairs",
        "route_compressed_bytes": "index.shape_route_compressed_bytes",
        "payloads_verified": "guard.shape_route_payloads_verified",
        "reconstruction_verified": "guard.shape_route_reconstruction_verified",
        "sidecar_bytes": "io.collection_sidecar_bytes_read",
        "identity_content_bytes": "io.source_identity_content_bytes_read",
        "identity_total_bytes": "io.source_identity_total_bytes_read",
        "verification_source_bytes": "io.route_verification_source_bytes_read",
    },
    "query_paths": {
        "source_identity_bytes": "planning.source_archive_identity_bytes_read",
        "source_execution_bytes": "planning.source_archive_execution_bytes_read",
        "source_bytes": "planning.source_archive_bytes_read",
        "sidecar_bytes": "planning.collection_sidecar_bytes_read",
        "route_sidecar_bytes": "planning.shape_route_sidecar_payload_bytes_read",
        "total_bytes": "planning.total_logical_bytes_read",
        "route_blocks_loaded": "planning.route_blocks_loaded",
        "routed_archives": "planning.routed_archives",
        "fallback_archives": "planning.fallback_archives",
        "archives_opened": "planning.archives_opened",
        "collection_load_seconds": "planning.collection_load_seconds",
        "identity_guard_seconds": "planning.identity_guard_seconds",
        "route_planning_seconds": "planning.route_planning_seconds",
        "source_execution_seconds": "planning.source_execution_seconds",
        "total_seconds": "planning.total_seconds",
    },
    "assumptions": [
        "collection build --shape-routes opts into route construction; omission is exact fallback",
        "route sections are canonical s.<local_archive>.<first_span> sections with at most 256 spans",
        "each exact-span bucket stores sorted unique (shape_id, donor_offset) pairs",
        "bindings record archive ordinal, source root, shapes digest, n_shapes, and sorted descriptors",
        "structured query JSON reports source, sidecar, and total logical bytes plus phase timings",
        "collection inspect --verify-routes independently reconstructs every route from source shapes",
    ],
}


@dataclass(frozen=True)
class QueryCase:
    case_id: str
    kind: str
    locus: str
    stratum: str | None = None


def exact_keys(value: Mapping[str, Any], keys: set[str], label: str) -> None:
    require(type(value) is dict and set(value) == keys, f"{label}: fields changed: {sorted(value)}")


def parse_sample(value: str) -> tuple[str, Path]:
    try:
        sample, raw_path = value.split("=", 1)
    except ValueError as error:
        raise ValueError("archive must be ID=PATH") from error
    require(sample in SAMPLES, f"archive ID must be one of {SAMPLES}")
    return sample, Path(raw_path)


def render(template: Sequence[str], values: Mapping[str, str]) -> list[str]:
    result = []
    for token in template:
        require(type(token) is str and token, "adapter tokens must be nonempty strings")
        try:
            rendered = token.format_map(values)
        except KeyError as error:
            raise GateError(f"unknown adapter placeholder {error.args[0]!r}") from error
        require("{" not in rendered and "}" not in rendered, f"unresolved adapter token {rendered!r}")
        result.append(rendered)
    return result


def validate_adapter(value: Any) -> dict[str, Any]:
    exact_keys(
        value,
        {
            "schema", "status", "build_fallback_extra", "build_route_extra",
            "query_fallback_extra", "query_route_extra", "inspect_fallback", "inspect_route",
            "inspect_route_full", "schemas", "build_paths", "inspect_paths", "query_paths",
            "assumptions",
        },
        "Gate-B CLI adapter",
    )
    require(value["schema"] == DEFAULT_ADAPTER["schema"], "Gate-B adapter schema differs")
    require(value["status"] in {"DRAFT", "FINAL"}, "Gate-B adapter status differs")
    for key in (
        "build_fallback_extra", "build_route_extra", "query_fallback_extra", "query_route_extra",
        "inspect_fallback", "inspect_route", "inspect_route_full", "assumptions",
    ):
        require(type(value[key]) is list and all(type(item) is str for item in value[key]), f"bad adapter {key}")
    for key in ("schemas", "build_paths", "inspect_paths", "query_paths"):
        require(type(value[key]) is dict and all(type(k) is str and type(v) is str for k, v in value[key].items()), f"bad adapter {key}")
        require(set(value[key]) == set(DEFAULT_ADAPTER[key]), f"Gate-B adapter {key} fields differ")
    return value


def adapter_is_frozen_final(adapter: Mapping[str, Any]) -> bool:
    """Only the exact source-controlled final adapter can open treatment."""
    return DEFAULT_ADAPTER["status"] == "FINAL" and adapter == DEFAULT_ADAPTER


def is_promotable_run(
    interface_status: str,
    *,
    adapter: Mapping[str, Any],
    code_dirty: bool,
    scripts_dirty: bool,
    fixture_hooks: Sequence[Path],
    required_fixture_hooks: Sequence[Path],
) -> bool:
    return (
        interface_status == "final"
        and adapter_is_frozen_final(adapter)
        and not code_dirty
        and not scripts_dirty
        and list(fixture_hooks) == list(required_fixture_hooks)
        and bool(required_fixture_hooks)
    )


def selection_digest(chrom: str, donor: int, acceptor: int) -> str:
    return blake3_bytes(SELECTION_DOMAIN + chrom.encode("utf-8") + struct.pack("<II", donor, acceptor))


def load_frozen_panel(panel_path: Path, metadata_path: Path) -> tuple[list[QueryCase], dict[str, Any]]:
    require(panel_path.is_file() and panel_path.stat().st_size == PANEL_BYTES, "frozen Gate-B panel size differs")
    require(sha256_file(panel_path) == PANEL_SHA256, "frozen Gate-B panel SHA-256 differs")
    require(metadata_path.is_file() and sha256_file(metadata_path) == PANEL_METADATA_SHA256, "panel metadata SHA-256 differs")
    metadata = load_json(metadata_path)
    require(metadata.get("schema") == "gravlax.jshape-prospective-panel.v1", "panel metadata schema differs")
    require(metadata.get("status") == "FROZEN_AFTER_GATE_A_PASS_BEFORE_GATE_B_IMPLEMENTATION_OR_MEASUREMENT", "panel freeze status differs")
    require(metadata.get("panel", {}).get("sha256") == PANEL_SHA256, "metadata panel digest differs")
    require(metadata.get("panel", {}).get("rows") == PANEL_ROWS, "metadata panel cardinality differs")
    require(metadata.get("eligible_rows") == {name: row[3] for name, row in STRATA.items()}, "eligible strata differ")

    with panel_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(
            reader.fieldnames
            == ["stratum", "rank", "selection_blake3", "chrom", "donor", "acceptor", "archives_present", "archive_ids"],
            "panel columns differ",
        )
        raw_rows = list(reader)
    require(len(raw_rows) == PANEL_ROWS, "panel row count differs")
    seen_coordinates = set()
    cases = []
    for index, row in enumerate(raw_rows, 1):
        stratum = row["stratum"]
        require(stratum in STRATA, f"panel row {index}: unknown stratum")
        minimum, maximum, expected_rows, _ = STRATA[stratum]
        donor, acceptor = int(row["donor"]), int(row["acceptor"])
        require(0 <= donor < acceptor <= 0xFFFFFFFF, f"panel row {index}: invalid coordinates")
        coordinate = (row["chrom"], donor, acceptor)
        require(coordinate not in seen_coordinates, f"panel row {index}: repeated coordinate")
        seen_coordinates.add(coordinate)
        ids = row["archive_ids"].split(",")
        count = int(row["archives_present"])
        require(ids == sorted(set(ids)) and all(sample in SAMPLES for sample in ids), f"panel row {index}: archive IDs differ")
        require(count == len(ids) and minimum <= count <= maximum, f"panel row {index}: stratum membership differs")
        require(row["selection_blake3"] == selection_digest(*coordinate), f"panel row {index}: selection hash differs")
        cases.append(QueryCase(f"panel-{index:03d}", "junction", f"{coordinate[0]}:{donor}-{acceptor}", stratum))
    for stratum, (_, _, expected_rows, _) in STRATA.items():
        selected = [row for row in raw_rows if row["stratum"] == stratum]
        require(len(selected) == expected_rows, f"{stratum}: panel cardinality differs")
        require([int(row["rank"]) for row in selected] == list(range(1, expected_rows + 1)), f"{stratum}: ranks differ")
        ordering = [(row["selection_blake3"], row["chrom"], int(row["donor"]), int(row["acceptor"])) for row in selected]
        require(ordering == sorted(ordering), f"{stratum}: selection order differs")
    return cases, metadata


def fixed_cases() -> list[QueryCase]:
    return [
        QueryCase("dense", "junction", LOCI["dense"]),
        QueryCase("sparse", "junction", LOCI["sparse"]),
        QueryCase("absent", "junction", LOCI["absent"]),
        QueryCase("region", "region", LOCI["region"]),
        QueryCase("jset", "jset", f"include={LOCI['dense']};exclude={LOCI['sparse']}"),
    ]


def query_command(
    aie: Path,
    collection: Path,
    case: QueryCase,
    arm: str,
    adapter: Mapping[str, Any],
) -> list[str]:
    require(arm in {"fallback", "route"}, "invalid Gate-B arm")
    if case.kind == "junction":
        command = [str(aie), "collection", "junction", str(collection), case.locus]
    elif case.kind == "region":
        command = [str(aie), "collection", "region", str(collection), case.locus]
    elif case.kind == "jset":
        command = [
            str(aie), "collection", "jset", str(collection),
            "--include", LOCI["dense"], "--exclude", LOCI["sparse"],
        ]
    else:
        raise GateError(f"unknown Gate-B query kind {case.kind}")
    command.extend(["--json", "--explain", "--top", "0"])
    command.extend(adapter[f"query_{arm}_extra"])
    return command


def build_command(
    aie: Path,
    archives: Mapping[str, Path],
    samples: Sequence[str],
    output: Path,
    arm: str,
    adapter: Mapping[str, Any],
    *,
    base: Path | None = None,
) -> list[str]:
    require(arm in {"fallback", "route"}, "invalid Gate-B build arm")
    command = [str(aie), "collection", "build"]
    if base is not None:
        command.extend(["--base", str(base)])
    for sample in samples:
        require(sample in archives, f"missing archive {sample}")
        command.extend(["--sample", f"{sample}={archives[sample]}"])
    command.extend(["--out", str(output)])
    command.extend(adapter[f"build_{arm}_extra"])
    return command


def inspect_command(aie: Path, collection: Path, mode: str, adapter: Mapping[str, Any]) -> list[str]:
    require(mode in {"fallback", "route", "route_full"}, "invalid inspect mode")
    return render(adapter[f"inspect_{mode}"], {"aie": str(aie), "collection": str(collection)})


_EXECUTION_FIELDS = {
    "planning",
    "explain",
    "actual_archive_bytes_read",
    "planned_compressed_bytes",
    "actual_collection_bytes_read",
    "total_logical_bytes_read",
    "source_sections_read",
    "collection_sections_read",
    "chunks_decoded",
    "unique_chunks_decoded",
    "independent_chunk_decodes",
    "junction_shape_route_used",
    "shape_route_used",
    "shape_route_section",
    "shape_route_sections",
    "shape_route_sidecar_bytes_read",
}


def scientific_projection(value: Any) -> Any:
    """Remove only paths, plans, byte accounting, and timings; retain every scientific value."""
    if isinstance(value, list):
        return [scientific_projection(item) for item in value]
    if isinstance(value, dict):
        return {
            key: scientific_projection(item)
            for key, item in value.items()
            if key not in _EXECUTION_FIELDS
            and not key.endswith("_seconds")
            and key not in {"archive", "collection_path"}
        }
    return value


def path_value(value: Mapping[str, Any], dotted: str, label: str) -> Any:
    current: Any = value
    for token in dotted.split("."):
        require(type(current) is dict and token in current, f"{label}: missing {dotted}")
        current = current[token]
    return current


def nonnegative_integer(value: Any, label: str) -> int:
    require(type(value) is int and value >= 0, f"{label}: expected nonnegative integer")
    return value


def nonnegative_number(value: Any, label: str) -> float:
    require(type(value) in {int, float} and math.isfinite(value) and value >= 0, f"{label}: expected nonnegative number")
    return float(value)


def canonical_route_section_name(value: str) -> bool:
    return re.fullmatch(r"s\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)", value) is not None


def build_schedule(case_ids: Iterable[str] | None = None) -> list[dict[str, Any]]:
    ids = tuple(case_ids) if case_ids is not None else ("build",)
    rows = []
    for case_id in ids:
        for block in range(1, WARM_BLOCKS + 1):
            order = ("fallback", "route") if block % 2 else ("route", "fallback")
            for position, arm in enumerate(order, 1):
                rows.append(
                    {
                        "case_id": case_id,
                        "block": block,
                        "position": position,
                        "arm": arm,
                        "label": f"{case_id}-b{block:02d}-p{position}-{arm}",
                    }
                )
    return rows


def nearest_rank(values: Sequence[float], quantile: float) -> float:
    require(bool(values) and 0 < quantile <= 1, "invalid nearest-rank inputs")
    ordered = sorted(values)
    index = math.ceil(quantile * len(ordered)) - 1
    return ordered[index]


def planned_output_directories(run: Path) -> tuple[Path, ...]:
    return (
        run,
        run / "audits",
        run / "collections",
        run / "collections" / "canonical",
        run / "collections" / "benchmark",
        run / "commands",
        run / "commands" / "build-canonical",
        run / "commands" / "build-benchmark",
        run / "commands" / "inspect",
        run / "commands" / "query-canonical",
        run / "commands" / "query-fixed-benchmark",
        run / "commands" / "query-panel-benchmark",
        run / "commands" / "fixture-hooks",
        run / "fixtures",
    )


def prepare_output_tree(run: Path) -> None:
    require(not run.exists(), f"refusing to overwrite {run}")
    for directory in planned_output_directories(run):
        directory.mkdir(parents=True, exist_ok=True)
    require(all(path.is_dir() for path in planned_output_directories(run)), "failed to create Gate-B output tree")
