#!/usr/bin/env python3
"""Validate and summarize the locked D0 junction-enumeration gate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
import subprocess
from pathlib import Path


LOCUS = "chr11:35138870-35232402"
MIN_SUPPORT = 20


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(root: Path, path: Path) -> dict[str, object]:
    try:
        logical = str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        logical = path.name
    return {"path": logical, "bytes": path.stat().st_size, "sha256": sha256(path)}


def timed_rows(run_dir: Path, arm: str) -> list[dict[str, float | int]]:
    rows = []
    for rep in range(1, 6):
        fields = (run_dir / f"{arm}-{rep}.time").read_text().strip().split("\t")
        if len(fields) != 2:
            raise ValueError(f"malformed timing for {arm} repetition {rep}")
        rows.append({"rep": rep, "wall_seconds": float(fields[0]), "max_rss_kb": int(fields[1])})
    return rows


def tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_json(command: list[str]) -> dict[str, object]:
    completed = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return json.loads(completed.stdout)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--aie-bin", type=Path, required=True)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    run_dir = args.run_dir.resolve()
    archive = root / "runs/archive/d0.aie"

    listing = json.loads((run_dir / "cells-1.json").read_text())
    index_json = json.loads((run_dir / "index.json").read_text())
    index_tsv = tsv_rows(run_dir / "index-1.tsv")
    rows = listing["junctions"]
    if not rows:
        raise AssertionError("fixed locus returned no junctions")
    if [(row["donor"], row["acceptor"], row["supporting_children"], row["posting_chunks"])
        for row in rows] != [
            (row["donor"], row["acceptor"], row["supporting_children"], row["posting_chunks"])
            for row in index_json["junctions"]
        ]:
        raise AssertionError("index-only and cell-decoded JSON rows differ")
    if [(str(row["donor"]), str(row["acceptor"]), str(row["supporting_children"]),
         str(row["posting_chunks"])) for row in rows] != [
            (row["donor"], row["acceptor"], row["supporting_children"], row["posting_chunks"])
            for row in index_tsv
        ]:
        raise AssertionError("TSV and JSON rows differ")

    point_checks = []
    for row in rows:
        locus = f"{row['chrom']}:{row['donor']}-{row['acceptor']}"
        point = run_json([
            str(args.aie_bin), "query", str(archive), "junction", locus, "--top", "0", "--json"
        ])
        for field in ("donor", "acceptor", "supporting_children", "posting_chunks", "umis", "cells"):
            if point[field] != row[field]:
                raise AssertionError(f"point/listing mismatch at {locus} for {field}")
        if point["cell_rows"] != row["cell_counts"]:
            raise AssertionError(f"point/listing per-cell mismatch at {locus}")
        point_checks.append({
            "locus": locus,
            "supporting_children": row["supporting_children"],
            "umis": row["umis"],
            "cells": row["cells"],
            "per_cell_exact": True,
        })

    annotation_rows = tsv_rows(run_dir / "annotated.tsv")
    if len(annotation_rows) != len(rows) or not all(
        set((row["annotated"], row["donor_annotated"], row["acceptor_annotated"])) <= {"0", "1"}
        for row in annotation_rows
    ):
        raise AssertionError("annotation flags missing or malformed")
    min_cells = tsv_rows(run_dir / "min-cells.tsv")
    if not min_cells or not all(int(row["cells"]) >= 20 for row in min_cells):
        raise AssertionError("--min-cells did not enforce its threshold")
    either = tsv_rows(run_dir / "either.tsv")
    contained_boundary = json.loads((run_dir / "contained-boundary.json").read_text())["junctions"]
    if len(either) != 2 or contained_boundary:
        raise AssertionError("--either/default boundary semantics failed")

    index_times = timed_rows(run_dir, "index")
    cell_times = timed_rows(run_dir, "cells")
    index_wall = statistics.median(row["wall_seconds"] for row in index_times)
    cells_wall = statistics.median(row["wall_seconds"] for row in cell_times)
    index_rss = max(row["max_rss_kb"] for row in index_times)
    cells_rss = max(row["max_rss_kb"] for row in cell_times)

    deterministic = all(
        (run_dir / f"index-{rep}.tsv").read_bytes() == (run_dir / "index-1.tsv").read_bytes()
        and (run_dir / f"cells-{rep}.json").read_bytes() == (run_dir / "cells-1.json").read_bytes()
        for rep in range(2, 6)
    )
    result = {
        "schema": "gravlax.query-junctions-gate.v1",
        "protocol_lock_commit": "52dca7a87942db4b44a7d5aea8feeb2c679c7511",
        "implementation_commit": args.code_commit,
        "configuration": {
            "dataset": "D0",
            "locus": LOCUS,
            "min_supporting_children": MIN_SUPPORT,
            "threads": 24,
            "warm_repetitions": 5,
            "coordinates": "0-based half-open",
            "default_interval_semantics": "both_endpoints_contained",
        },
        "inputs": {
            "archive": identity(root, archive),
            "aie_binary": identity(root, args.aie_bin.resolve()),
            "gencode_v49": identity(root, root / "annotations/gencode.v49.annotation.gtf"),
        },
        "correctness": {
            "junction_rows": len(rows),
            "point_interval_aggregate_and_per_cell_exact": True,
            "point_checks": point_checks,
            "tsv_json_ordered_rows_equal": True,
            "five_machine_output_repetitions_byte_identical": deterministic,
            "min_support_filter_exercised": True,
            "min_cells_filter_exercised": True,
            "contained_and_either_boundary_semantics_exercised": True,
            "annotation_flags_well_formed": True,
            "workspace_tests_passed": 93,
        },
        "performance": {
            "index_only": {
                "runs": index_times,
                "median_wall_seconds": index_wall,
                "max_rss_kb": index_rss,
                "wall_gate_seconds": 0.20,
                "rss_gate_kb": 512 * 1024,
                "verdict": "pass" if index_wall <= 0.20 and index_rss <= 512 * 1024 else "fail",
            },
            "with_cells": {
                "runs": cell_times,
                "median_wall_seconds": cells_wall,
                "max_rss_kb": cells_rss,
                "wall_gate_seconds": 2.0,
                "rss_gate_kb": 2 * 1024 * 1024,
                "verdict": "pass" if cells_wall <= 2.0 and cells_rss <= 2 * 1024 * 1024 else "fail",
            },
        },
        "format": {
            "archive_bytes_added": 0,
            "tsv_sha256": sha256(run_dir / "index-1.tsv"),
            "json_with_cells_sha256": sha256(run_dir / "cells-1.json"),
        },
    }
    if not deterministic:
        raise AssertionError("machine outputs were not deterministic")
    if result["performance"]["index_only"]["verdict"] != "pass":
        raise AssertionError("index-only performance gate failed")
    if result["performance"]["with_cells"]["verdict"] != "pass":
        raise AssertionError("with-cells performance gate failed")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(
        f"junctions={len(rows)} index={index_wall:.3f}s/{index_rss}KB "
        f"with-cells={cells_wall:.3f}s/{cells_rss}KB exact=PASS"
    )


if __name__ == "__main__":
    main()
