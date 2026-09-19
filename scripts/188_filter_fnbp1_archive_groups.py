#!/usr/bin/env python3
"""Restrict frozen H5 labels to the cells present in the targeted D5 archive."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_groups(path: Path) -> dict[str, str]:
    groups: dict[str, str] = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split("\t")
        if len(fields) != 2:
            raise SystemExit(f"{path}:{number}: expected barcode<TAB>group")
        barcode, group = fields
        if barcode in groups:
            raise SystemExit(f"duplicate barcode in {path}: {barcode}")
        groups[barcode] = group
    return groups


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-map", type=Path, required=True)
    parser.add_argument("--archive-cells", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    groups = read_groups(args.full_map)
    archive = json.loads(args.archive_cells.read_text())
    if archive.get("schema") != "gravlax.query.region.v1":
        raise SystemExit("archive-cell enumeration is not a region-query result")
    if archive.get("cell_rows_truncated"):
        raise SystemExit("archive-cell enumeration was truncated")
    rows = archive["cell_rows"]
    archive_barcodes = [row["barcode"] for row in rows]
    if len(archive_barcodes) != archive["scope"]["selected_cells"]:
        raise SystemExit("whole-chromosome region query did not enumerate the full archive dictionary")
    if len(set(archive_barcodes)) != len(archive_barcodes):
        raise SystemExit("duplicate archive barcodes")
    missing = sorted(set(archive_barcodes) - groups.keys())
    if missing:
        raise SystemExit(f"{len(missing)} archive barcodes lack frozen H5 labels; first={missing[0]}")

    selected = set(archive_barcodes)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as handle:
        for barcode, group in groups.items():
            if barcode in selected:
                handle.write(f"{barcode}\t{group}\n")

    full_counts = Counter(groups.values())
    selected_counts = Counter(groups[barcode] for barcode in archive_barcodes)
    summary = {
        "schema": "gravlax.fnbp1-broad-cell-archive-scope.v1",
        "status": "FROZEN_LABELS_RESTRICTED_TO_EXISTING_ARCHIVE_DICTIONARY",
        "full_H5_cells": len(groups),
        "archive_cells": len(archive_barcodes),
        "full_H5_group_sizes": dict(sorted(full_counts.items())),
        "archive_selected_group_sizes": dict(sorted(selected_counts.items())),
        "missing_from_targeted_archive_by_group": {
            group: full_counts[group] - selected_counts[group] for group in sorted(full_counts)
        },
        "rate_denominator": "full_H5_group_sizes",
        "full_group_map_sha256": sha256(args.full_map),
        "archive_group_map": args.out.name,
        "archive_group_map_sha256": sha256(args.out),
        "archive_cell_enumeration_sha256": sha256(args.archive_cells),
    }
    summary_path = args.out.with_name("archive-scope-summary.json")
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
