#!/usr/bin/env python3
"""Compare every batched region/junction answer with its standalone command."""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path


REGION_SUMMARY = re.compile(
    r"^region (\S+): (\d+) molecules, (\d+) UMIs across (\d+) cells"
)


def run_region(aie: Path, archive: Path, locus: str, top: int) -> dict[str, object]:
    completed = subprocess.run(
        [str(aie), "query", str(archive), "region", locus, "--top", str(top)],
        check=True,
        text=True,
        capture_output=True,
    )
    lines = completed.stdout.splitlines()
    match = REGION_SUMMARY.match(lines[0])
    if match is None:
        raise SystemExit(f"could not parse region summary: {lines[0]!r}")
    cells = [
        {"barcode": fields[0], "umis": int(fields[1])}
        for line in lines[1:]
        if (fields := line.strip().split("\t")) and len(fields) == 2
    ]
    return {
        "molecules": int(match.group(2)),
        "umis": int(match.group(3)),
        "cells": int(match.group(4)),
        "cell_rows": cells,
    }


def run_junction(aie: Path, archive: Path, locus: str, top: int) -> dict[str, object]:
    completed = subprocess.run(
        [
            str(aie), "query", str(archive), "junction", locus,
            "--top", str(top), "--json",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    row = json.loads(completed.stdout)
    return {
        "umis": row["umis"],
        "cells": row["cells"],
        "supporting_children": row["supporting_children"],
        "posting_chunks": row["posting_chunks"],
        "cell_rows": row["cell_rows"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aie", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    with args.plan.open(newline="") as handle:
        plan = list(csv.DictReader(handle, delimiter="\t"))
    batch = json.loads(args.batch.read_text())
    if batch.get("schema") != "gravlax.query.batch.v1":
        raise SystemExit("unexpected batch schema")
    observed = batch["queries"]
    if [row["id"] for row in observed] != [row["id"] for row in plan]:
        raise SystemExit("batch output does not preserve plan order")

    for spec, actual in zip(plan, observed, strict=True):
        if spec["kind"] == "region":
            expected = run_region(args.aie, args.archive, spec["locus"], args.top)
            keys = ("molecules", "umis", "cells", "cell_rows")
        else:
            expected = run_junction(args.aie, args.archive, spec["locus"], args.top)
            keys = ("umis", "cells", "supporting_children", "posting_chunks", "cell_rows")
        for key in keys:
            if actual[key] != expected[key]:
                raise SystemExit(f"{spec['id']} differs from standalone {spec['kind']} at {key}")

    args.out.write_text(json.dumps({
        "schema": "gravlax.batch-correctness.v1",
        "queries_checked": len(plan),
        "aggregate_counts_exact": True,
        "returned_per_cell_rows_exact": True,
        "plan_order_preserved": True,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
