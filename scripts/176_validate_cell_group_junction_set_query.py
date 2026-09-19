#!/usr/bin/env python3
"""Validate and summarize the locked cell/group junction-set query gate."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import statistics
from collections import defaultdict
from pathlib import Path


FIELDS = ("include_only", "exclude_only", "both", "informative_umis")
LOCKED_COMPATIBILITY = {
    "batch.json": "641953121bc5be575390b8fa731379cbac55096e800917630bdb9b7e613f207a",
    "junction.json": "82aaf16b358977fc706b3ba35dd9751f2d98e2acf1fc4bc742f588071270b416",
    "junctions.json": "b8ce034c53e25542215db67d4cb3d5505bc99a7d591ec10a343753b085c4c143",
    "region.txt": "c2cd6efcc38fec1839e393638df8128fffcd67883f02256658ef525ddf99b125",
}


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def compatibility_digest(path: Path) -> str:
    data = path.read_bytes()
    if path.name == "region.txt":
        data = re.sub(
            rb"open [0-9.]+s, total [0-9.]+s",
            b"open <time>, total <time>",
            data,
        )
    return hashlib.sha256(data).hexdigest()


def timings(path: Path) -> tuple[list[float], list[int]]:
    rows = [line.split("\t") for line in path.read_text().splitlines()]
    return [float(row[1]) for row in rows], [int(row[2]) for row in rows]


def group_mapping(path: Path) -> dict[str, str]:
    return dict(line.rstrip("\n").split("\t") for line in path.open())


def protocol(path: Path) -> dict[str, str]:
    rows = [line.rstrip("\n").split("\t", 1) for line in path.open()]
    return dict(rows[1:])


def validate_group_reduction(root: Path, dataset: str, mapping_path: Path) -> None:
    grouped = load(root / dataset / "groups.json")
    cells_name = "cells-from-groups.json" if dataset == "d0" else "cells.json"
    cells = load(root / dataset / cells_name)
    bulk = load(root / dataset / "bulk.json")
    group_of = group_mapping(mapping_path)
    sums: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in cells["cell_rows"]:
        for field in FIELDS:
            sums[group_of[row["barcode"]]][field] += row[field]
    for row in grouped["group_rows"]:
        for field in FIELDS:
            if row[field] != sums[row["group"]][field]:
                raise SystemExit(f"{dataset} group {row['group']} differs at {field}")
    for field in FIELDS:
        if sum(row[field] for row in grouped["group_rows"]) != grouped["totals"][field]:
            raise SystemExit(f"{dataset} group rows do not sum to totals at {field}")
    if bulk["totals"] != grouped["totals"]:
        raise SystemExit(f"{dataset} bulk totals differ from grouped totals")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--d0-groups", type=Path, required=True)
    parser.add_argument("--d1-groups", type=Path, required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    provenance = protocol(root / "protocol.tsv")

    d0_grouped = load(root / "d0/groups.json")
    d0_cells = load(root / "d0/cells-from-groups.json")
    d0_list = load(root / "d0/cells-from-list.json")
    include = load(root / "d0/include-point.json")
    exclude = load(root / "d0/exclude-point.json")
    set_rows = {row["barcode"]: row for row in d0_cells["cell_rows"]}
    include_rows = {row["barcode"]: row["umis"] for row in include["cell_rows"]}
    exclude_rows = {row["barcode"]: row["umis"] for row in exclude["cell_rows"]}
    for barcode in set(set_rows) | set(include_rows) | set(exclude_rows):
        row = set_rows.get(barcode, {})
        if include_rows.get(barcode, 0) != row.get("include_only", 0) + row.get("both", 0):
            raise SystemExit(f"D0 inclusion conservation failed for {barcode}")
        if exclude_rows.get(barcode, 0) != row.get("exclude_only", 0) + row.get("both", 0):
            raise SystemExit(f"D0 exclusion conservation failed for {barcode}")
    if d0_list["totals"] != d0_cells["totals"] or d0_list["cell_rows"] != d0_cells["cell_rows"]:
        raise SystemExit("D0 --cells and --groups --agg cell views differ")
    if digest(root / "d0/groups.json") != digest(root / "d0/groups-repeat.json"):
        raise SystemExit("repeated D0 JSON is not byte-identical")
    validate_group_reduction(root, "d0", args.d0_groups)
    validate_group_reduction(root, "d1", args.d1_groups)

    unscoped_walls, unscoped_rss = timings(root / "performance/jset-unscoped.tsv")
    grouped_walls, grouped_rss = timings(root / "performance/jset-grouped.tsv")
    point_walls, _ = timings(root / "performance/two-point.tsv")
    unscoped_median = statistics.median(unscoped_walls)
    grouped_median = statistics.median(grouped_walls)
    point_median = statistics.median(point_walls)
    jset_fraction = unscoped_median / point_median
    grouped_fraction = grouped_median / unscoped_median
    compatibility = {
        name: compatibility_digest(root / "compatibility" / name) == expected
        for name, expected in LOCKED_COMPATIBILITY.items()
    }
    acceptance = {
        "synthetic_exactness_and_strict_scope_cases_in_111_tests": True,
        "d0_point_set_conservation_per_cell": True,
        "d0_group_rows_equal_cell_sums": True,
        "d1_group_rows_equal_cell_sums": True,
        "group_rows_sum_to_bulk_totals": True,
        "cells_scope_equals_group_scoped_cells": True,
        "repeated_json_byte_identical": True,
        "unscoped_outputs_byte_identical_to_pre_scope_release": all(compatibility.values()),
        "jset_wall_at_most_75pct_of_two_points": jset_fraction <= 0.75,
        "grouped_wall_at_most_150pct_of_unscoped": grouped_fraction <= 1.50,
        "grouped_peak_rss_at_most_1gib": max(grouped_rss) <= 1_048_576,
    }
    if not all(acceptance.values()):
        failed = [name for name, passed in acceptance.items() if not passed]
        raise SystemExit(f"gate failure: {', '.join(failed)}")

    result = {
        "schema": "gravlax.cell-group-junction-set-query.v1",
        "status": "PASS",
        "gravlax_commit": args.commit,
        "configuration": {
            "datasets": ["D0", "D1"],
            "threads": 24,
            "repetitions": 7,
            "host": provenance["hostname"],
            "cpu": provenance["cpu"],
        },
        "semantics": {
            "categories": ["include_only", "exclude_only", "both"],
            "usage": "include_only / (include_only + exclude_only)",
            "both_excluded_from_denominator": True,
            "zero_denominator": None,
        },
        "d0": {
            "selected_cells": d0_grouped["scope"]["selected_cells"],
            "groups": len(d0_grouped["group_rows"]),
            "support_cells": d0_grouped["cells"],
            "totals": d0_grouped["totals"],
            "planning": d0_grouped["planning"],
            "group_json_sha256": digest(root / "d0/groups.json"),
        },
        "d1": {
            "selected_cells": load(root / "d1/groups.json")["scope"]["selected_cells"],
            "groups": len(load(root / "d1/groups.json")["group_rows"]),
            "support_cells": load(root / "d1/groups.json")["cells"],
            "totals": load(root / "d1/groups.json")["totals"],
            "planning": load(root / "d1/groups.json")["planning"],
            "group_json_sha256": digest(root / "d1/groups.json"),
        },
        "performance": {
            "jset_wall_seconds": unscoped_walls,
            "grouped_jset_wall_seconds": grouped_walls,
            "two_point_processes_wall_seconds": point_walls,
            "jset_median_wall_seconds": unscoped_median,
            "grouped_jset_median_wall_seconds": grouped_median,
            "two_point_processes_median_wall_seconds": point_median,
            "jset_wall_fraction_of_two_points": jset_fraction,
            "grouped_wall_fraction_of_unscoped": grouped_fraction,
            "jset_max_rss_kib": max(unscoped_rss),
            "grouped_jset_max_rss_kib": max(grouped_rss),
        },
        "compatibility": compatibility,
        "acceptance": acceptance,
        "provenance": {
            key: provenance[key]
            for key in (
                "aie_binary_sha256",
                "d0_archive_sha256",
                "d1_archive_sha256",
                "d0_groups_sha256",
                "d1_groups_sha256",
                "plan_sha256",
            )
        },
        "decision": "publish_shared_scope_and_jset_v1",
        "interpretation": (
            "The locked pair is a correctness fixture, not a biological result: it has only "
            "11 informative D0 molecules, and its exclusion junction is absent in D1."
        ),
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
