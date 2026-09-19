#!/usr/bin/env python3
"""Summarize the registered eight-donor fixed-event validation."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from statistics import median


EVENT_ID = "cassette:chr9:129908999-129923843:129924026-129924959:129908999-129924959"
JUNCTIONS = {
    "include_left": ("chr9", 129909000, 129923843),
    "include_right": ("chr9", 129924027, 129924959),
    "skip": ("chr9", 129909000, 129924959),
}
RESOLVED_OTHER = {"astrocyte", "neuronal", "oligodendrocyte_OPC", "vascular"}
CATEGORIES = ("include_only", "exclude_only", "both")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def groups(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            barcode, group = line.rstrip().split("\t")
            if barcode in result:
                raise SystemExit(f"duplicate group barcode in {path}:{number}")
            result[barcode] = group
    return result


def star_junction_audit(directory: Path, frozen_groups: dict[str, str]) -> dict[str, object]:
    raw = directory / "Solo.out/SJ/raw"
    selected_rows: dict[int, str] = {}
    by_coordinate = {coordinate: name for name, coordinate in JUNCTIONS.items()}
    with (raw / "features.tsv").open() as handle:
        for index, line in enumerate(handle, 1):
            fields = line.split()
            coordinate = (fields[0], int(fields[1]), int(fields[2]))
            if coordinate in by_coordinate:
                selected_rows[index] = by_coordinate[coordinate]
    feature_present = {
        name: name in selected_rows.values() for name in JUNCTIONS
    }

    selected_columns: dict[int, str] = {}
    with (raw / "barcodes.tsv").open() as handle:
        for index, line in enumerate(handle, 1):
            barcode = line.strip().removesuffix("-1")
            if barcode in frozen_groups:
                selected_columns[index] = frozen_groups[barcode]
    if len(selected_columns) != len(frozen_groups):
        raise SystemExit(
            f"STAR whitelist contains {len(selected_columns)}/{len(frozen_groups)} frozen barcodes"
        )

    totals = {name: 0 for name in JUNCTIONS}
    by_group: dict[str, dict[str, int]] = {
        group: {name: 0 for name in JUNCTIONS} for group in sorted(set(frozen_groups.values()))
    }
    dimensions_seen = False
    with (raw / "matrix.mtx").open() as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.split()
            if not dimensions_seen:
                dimensions_seen = True
                continue
            row, column, value = map(int, fields)
            if row in selected_rows and column in selected_columns:
                name = selected_rows[row]
                group = selected_columns[column]
                totals[name] += value
                by_group[group][name] += value
    return {
        "called_nucleus_UMIs": totals,
        "by_group": by_group,
        "feature_present": feature_present,
    }


def add_counts(rows: list[dict[str, object]]) -> dict[str, object]:
    result = {name: sum(int(row[name]) for row in rows) for name in CATEGORIES}
    result["informative_umis"] = result["include_only"] + result["exclude_only"]
    result["cells"] = sum(int(row["cells"]) for row in rows)
    result["usage_fraction"] = (
        result["include_only"] / result["informative_umis"]
        if result["informative_umis"]
        else None
    )
    return result


def sign_p(negative: int, positive: int) -> float:
    n = negative + positive
    tail = min(negative, positive)
    return min(1.0, 2.0 * sum(math.comb(n, value) for value in range(tail + 1)) / (2**n))


def parse_time(path: Path) -> dict[str, object]:
    result: dict[str, object] = {}
    for line in path.read_text().splitlines():
        if "Maximum resident set size (kbytes):" in line:
            result["peak_rss_kib"] = int(line.rsplit(":", 1)[1])
        elif "Elapsed (wall clock) time" in line:
            value = line.rsplit(": ", 1)[1].strip()
            parts = value.split(":")
            seconds = 0.0
            for part in parts:
                seconds = seconds * 60 + float(part)
            result["wall_seconds"] = seconds
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--feasibility", type=Path, required=True)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-tsv", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    feasibility = json.loads(args.feasibility.read_text())
    acquisition = json.loads(args.acquisition.read_text())
    if feasibility["status"] != "PASS" or feasibility["event_evidence_inspected"]:
        raise SystemExit("event-blind feasibility did not authorize this run")
    if acquisition["status"] != "PASS" or acquisition["event_evidence_inspected"]:
        raise SystemExit("FASTQ acquisition is incomplete or was not event-blind")

    donor_rows: list[dict[str, object]] = []
    donor_results: dict[str, dict[str, object]] = {}
    for donor in "ABCDEFGH":
        donor_dir = args.run / f"donor-{donor}"
        query = json.loads((donor_dir / "jset-by-group.json").read_text())
        if query["schema"] != "gravlax.query.jset.v1":
            raise SystemExit(f"unexpected jset schema for donor {donor}")
        full_groups = groups(
            args.feasibility.parent / "groups" / f"donor-{donor}.tsv"
        )
        full_sizes = feasibility["donors"][donor]["groups"]
        query_rows = {row["group"]: row for row in query["group_rows"]}
        normalized_rows: list[dict[str, object]] = []
        for group in sorted(full_sizes):
            row = query_rows.get(
                group,
                {
                    "group": group,
                    "selected_cells": 0,
                    "cells": 0,
                    "include_only": 0,
                    "exclude_only": 0,
                    "both": 0,
                    "informative_umis": 0,
                    "usage_fraction": None,
                },
            ).copy()
            row["full_group_nuclei"] = int(full_sizes[group])
            normalized_rows.append(row)
            donor_rows.append({"donor": donor, **row})
        summed = add_counts(normalized_rows)
        for field in (*CATEGORIES, "informative_umis"):
            if int(summed[field]) != int(query["totals"][field]):
                raise SystemExit(f"donor {donor} group sum differs from bulk {field}")

        micro = add_counts([row for row in normalized_rows if row["group"] == "microglia_immune"])
        other = add_counts([row for row in normalized_rows if row["group"] in RESOLVED_OTHER])
        eligible_bulk = int(summed["informative_umis"]) >= 10
        eligible_micro = int(micro["informative_umis"]) >= 4 and int(other["informative_umis"]) >= 10
        delta = (
            float(micro["usage_fraction"]) - float(other["usage_fraction"])
            if eligible_micro
            else None
        )
        log = (donor_dir / "star/Log.final.out").read_text()
        input_reads_match = re.search(r"Number of input reads\s*\|\s*(\d+)", log)
        donor_results[donor] = {
            "published_QC_nuclei": int(feasibility["donors"][donor]["nuclei"]),
            "input_read_pairs": int(input_reads_match.group(1)) if input_reads_match else None,
            "bulk": summed,
            "microglia_immune": micro,
            "other_resolved": other,
            "eligible_primary": eligible_bulk,
            "eligible_microglia": eligible_micro,
            "microglia_minus_other_usage": delta,
            "group_rows": normalized_rows,
            "STAR_SJ": star_junction_audit(donor_dir / "star", full_groups),
            "resources": {
                "alignment": parse_time(donor_dir / "alignment.time.txt"),
                "ingest": parse_time(donor_dir / "ingest.time.txt"),
                "query": parse_time(donor_dir / "query.time.txt"),
                "target_bam_bytes": (donor_dir / "fnbp1.target.bam").stat().st_size,
                "archive_bytes": (donor_dir / "fnbp1.called.aie").stat().st_size,
            },
            "digests": {
                "target_bam_sha256": sha256(donor_dir / "fnbp1.target.bam"),
                "archive_sha256": sha256(donor_dir / "fnbp1.called.aie"),
                "jset_sha256": sha256(donor_dir / "jset-by-group.json"),
            },
        }

    eligible_primary = [row for row in donor_results.values() if row["eligible_primary"]]
    primary_usages = [float(row["bulk"]["usage_fraction"]) for row in eligible_primary]
    if len(eligible_primary) < 4:
        primary_status = "STOP_UNDERPOWERED"
    elif sum(value > 0.5 for value in primary_usages) / len(primary_usages) >= 0.75 and median(primary_usages) >= 0.75:
        primary_status = "PASS"
    else:
        primary_status = "FAIL"
    primary = {
        "status": primary_status,
        "eligible_donors": len(eligible_primary),
        "donors_usage_above_0.50": sum(value > 0.5 for value in primary_usages),
        "fraction_usage_above_0.50": (
            sum(value > 0.5 for value in primary_usages) / len(primary_usages)
            if primary_usages else None
        ),
        "median_donor_usage": median(primary_usages) if primary_usages else None,
    }

    eligible_micro = [row for row in donor_results.values() if row["eligible_microglia"]]
    deltas = [float(row["microglia_minus_other_usage"]) for row in eligible_micro]
    nonzero = [value for value in deltas if value != 0]
    negative = sum(value < 0 for value in deltas)
    positive = sum(value > 0 for value in deltas)
    ties = len(deltas) - len(nonzero)
    if len(eligible_micro) < 6 or len(nonzero) < 6:
        micro_status = "STOP_DESCRIPTIVE"
        p_value = None
    else:
        p_value = sign_p(negative, positive)
        micro_status = (
            "PASS"
            if negative / len(eligible_micro) >= 0.75 and p_value <= 0.05
            else "FAIL"
        )
    micro = {
        "status": micro_status,
        "eligible_donors": len(eligible_micro),
        "nonzero_differences": len(nonzero),
        "negative_differences": negative,
        "positive_differences": positive,
        "ties": ties,
        "fraction_negative_among_eligible": negative / len(eligible_micro) if eligible_micro else None,
        "two_sided_exact_sign_test_p": p_value,
        "stop_reason": (
            "fewer than six eligible donors or fewer than six nonzero paired differences"
            if micro_status == "STOP_DESCRIPTIVE" else None
        ),
    }

    result = {
        "schema": "gravlax.fnbp1-multidonor-sez-validation.v1",
        "status": primary_status,
        "locked_event": EVENT_ID,
        "unit_of_replication": "donor",
        "registration": {
            "gate_commit": "95b79ad9190063d3a6702f8bdfab7d3cf9cd1cb6",
            "gate_commit_time": "2026-09-01T11:01:55-04:00",
            "feasibility_commit": "9683cb9ae0b9d9b220ed5bd41246e14253652983",
            "raw_acquisition_started_after_both_commits": True,
        },
        "inputs": {
            "manifest_sha256": sha256(args.manifest),
            "feasibility_sha256": sha256(args.feasibility),
            "acquisition_sha256": sha256(args.acquisition),
            "executable_sha256": sha256(args.executable),
        },
        "primary_brain_replication": primary,
        "microglia_localization": micro,
        "donors": donor_results,
        "interpretation_limits": [
            "no mechanism",
            "no causal age effect",
            "no generalization outside adult human SEZ",
            "cells and molecules are not replicates",
            "STAR-SJ is an independent reduction of the same annotation-free alignments, not an independent mapper",
        ],
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    with args.out_tsv.open("w") as handle:
        fields = [
            "donor", "group", "full_group_nuclei", "selected_cells", "cells",
            "include_only", "exclude_only", "both", "informative_umis", "usage_fraction",
        ]
        handle.write("\t".join(fields) + "\n")
        for row in donor_rows:
            handle.write(
                "\t".join(
                    "NA" if row.get(field) is None else str(row.get(field))
                    for field in fields
                )
                + "\n"
            )
    print(json.dumps({
        "status": result["status"],
        "primary_brain_replication": primary,
        "microglia_localization": micro,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
