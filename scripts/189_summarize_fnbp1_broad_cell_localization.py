#!/usr/bin/env python3
"""Summarize the frozen-label, fixed-event D5 broad-cell localization."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path


EVENT_ID = "cassette:chr9:129908999-129923843:129924026-129924959:129908999-129924959"
GROUP_ORDER = [
    "neuronal",
    "oligodendrocyte_OPC",
    "astrocyte",
    "microglia_immune",
    "vascular",
    "unresolved",
    "ambiguous_or_doublet",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def wilson(successes: int, total: int) -> list[float] | None:
    if total == 0:
        return None
    z = 1.959963984540054
    z2 = z * z
    p = successes / total
    scale = 1.0 + z2 / total
    center = (p + z2 / (2.0 * total)) / scale
    half = z * math.sqrt(p * (1.0 - p) / total + z2 / (4.0 * total * total)) / scale
    return [max(0.0, center - half), min(1.0, center + half)]


def timing(path: Path) -> dict[str, float | int]:
    text = path.read_text()
    elapsed = re.search(r"Elapsed \(wall clock\) time.*: ([0-9:.]+)", text)
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if not elapsed or not rss:
        raise SystemExit(f"could not parse GNU time output: {path}")
    pieces = [float(piece) for piece in elapsed.group(1).split(":")]
    seconds = pieces[-1]
    if len(pieces) >= 2:
        seconds += 60.0 * pieces[-2]
    if len(pieces) == 3:
        seconds += 3600.0 * pieces[0]
    return {"wall_seconds": seconds, "peak_rss_kib": int(rss.group(1))}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--confirmation-summary", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-tsv", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    groups = json.loads((args.run / "groups-summary.json").read_text())
    scope = json.loads((args.run / "archive-scope-summary.json").read_text())
    query = json.loads((args.run / "events-by-broad-cell.json").read_text())
    confirmation = json.loads(args.confirmation_summary.read_text())
    if groups["manifest_sha256"] != sha256(args.manifest):
        raise SystemExit("group map was not built from the supplied locked manifest")
    if groups["group_map_sha256"] != scope["full_group_map_sha256"]:
        raise SystemExit("full group-map digest mismatch")
    if query["scope"]["selected_cells"] != scope["archive_cells"]:
        raise SystemExit("query scope does not equal the enumerated archive dictionary")

    matches = [event for event in query["events"] if event["id"] == EVENT_ID]
    if len(matches) != 1:
        raise SystemExit(f"expected one locked event, found {len(matches)}")
    event = matches[0]
    expected_totals = confirmation["Cell_Ranger_alignment_archive"]["totals"]
    for key in ("include_only", "exclude_only", "both", "informative_umis"):
        if event["totals"][key] != expected_totals[key]:
            raise SystemExit(f"bulk total mismatch for {key}")

    query_rows = {row["group"]: row for row in event["group_rows"]}
    if set(query_rows) != set(scope["full_H5_group_sizes"]):
        raise SystemExit("query groups do not match frozen H5 groups")
    rows = []
    for group in GROUP_ORDER:
        row = query_rows[group]
        full_cells = int(scope["full_H5_group_sizes"][group])
        archive_cells = int(scope["archive_selected_group_sizes"][group])
        include = int(row["include_only"])
        exclude = int(row["exclude_only"])
        both = int(row["both"])
        informative = int(row["informative_umis"])
        event_cells = int(row["cells"])
        if row["selected_cells"] != archive_cells:
            raise SystemExit(f"archive-selected cell mismatch for {group}")
        rows.append(
            {
                "group": group,
                "full_H5_nuclei": full_cells,
                "archive_selected_nuclei": archive_cells,
                "event_bearing_cells": event_cells,
                "event_bearing_cells_per_1000_full_H5_nuclei": 1000.0 * event_cells / full_cells,
                "event_bearing_fraction_wilson_95": wilson(event_cells, full_cells),
                "include_only": include,
                "exclude_only": exclude,
                "both": both,
                "informative_umis": informative,
                "inclusion_usage": None if informative == 0 else include / informative,
                "inclusion_usage_wilson_95": wilson(include, informative),
            }
        )

    totals = event["totals"]
    result = {
        "schema": "gravlax.fnbp1-broad-cell-localization.v1",
        "status": "DESCRIPTIVE_SINGLE_DONOR_LOCALIZATION_COMPLETE",
        "dataset": manifest["dataset"],
        "locked_event": EVENT_ID,
        "label_lock": {
            "manifest": "experiments/manifests/fnbp1-broad-cell-localization.json",
            "manifest_sha256": sha256(args.manifest),
            "group_map_sha256": groups["group_map_sha256"],
            "archive_group_map_sha256": scope["archive_group_map_sha256"],
            "labels_frozen_before_event_query": True,
            "FNBP1_used_for_labels": False,
        },
        "scope": {
            "full_H5_nuclei": scope["full_H5_cells"],
            "targeted_archive_nuclei": scope["archive_cells"],
            "rate_denominator": "full H5-labelled group size",
        },
        "bulk_exact_reproduction": {
            "include_only": totals["include_only"],
            "exclude_only": totals["exclude_only"],
            "both": totals["both"],
            "informative_umis": totals["informative_umis"],
            "inclusion_usage": totals["usage_fraction"],
        },
        "groups": rows,
        "execution": {
            "query": timing(args.run / "query.time.txt"),
            "executable_sha256": sha256(args.executable),
            "raw_query_output_sha256": sha256(args.run / "events-by-broad-cell.json"),
        },
        "interpretation": {
            "descriptive_observation": "both exclusion-only molecules occur in the broad microglia/immune group, which also contains two of the 72 inclusion-only molecules",
            "localization": "event-bearing cells occur in every broad compartment, with the largest absolute count in oligodendrocyte/OPC nuclei",
            "power": "only 74 informative molecules occur in one donor; microglia/immune has four informative molecules and a 0.50 inclusion estimate with a wide Wilson interval",
            "label_limit": "the conservative frozen marker rule leaves 3,533 of 10,261 nuclei unresolved and marks 303 ambiguous/doublet; labels are broad marker-rule assignments, not externally curated truth",
            "claims_forbidden": manifest["reporting"]["claims_forbidden"],
        },
    }
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    with args.out_tsv.open("w") as handle:
        handle.write(
            "group\tfull_H5_nuclei\tarchive_selected_nuclei\tevent_bearing_cells\t"
            "event_bearing_cells_per_1000_full_H5_nuclei\tinclude_only\texclude_only\tboth\t"
            "informative_umis\tinclusion_usage\tinclusion_usage_wilson_95_low\t"
            "inclusion_usage_wilson_95_high\n"
        )
        for row in rows:
            interval = row["inclusion_usage_wilson_95"] or [None, None]
            values = [
                row["group"], row["full_H5_nuclei"], row["archive_selected_nuclei"],
                row["event_bearing_cells"], f"{row['event_bearing_cells_per_1000_full_H5_nuclei']:.6f}",
                row["include_only"], row["exclude_only"], row["both"], row["informative_umis"],
                "NA" if row["inclusion_usage"] is None else f"{row['inclusion_usage']:.6f}",
                "NA" if interval[0] is None else f"{interval[0]:.6f}",
                "NA" if interval[1] is None else f"{interval[1]:.6f}",
            ]
            handle.write("\t".join(map(str, values)) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
