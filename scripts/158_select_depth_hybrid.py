#!/usr/bin/env python3
"""Select the locked D0 monotone depth-hybrid transition scale."""
from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path


def load_mode(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    modes = payload.get("modes", [])
    if len(modes) != 1 or modes[0].get("name") != "depth-hybrid":
        raise SystemExit(f"expected one depth-hybrid mode in {path}")
    return modes[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    grouped: dict[float, list[dict[str, object]]] = defaultdict(list)
    with (args.root / "runs.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["model"] == "depth-hybrid" and row["control"] == "real":
                grouped[float(row["depth_scale"])].append(
                    load_mode(Path(row["out"]) / "metrics.json")
                )

    records = []
    for scale, modes in sorted(grouped.items()):
        if len(modes) != 3:
            raise SystemExit(f"scale {scale} has {len(modes)} seeds")
        records.append(
            {
                "depth_scale": scale,
                "seeds": len(modes),
                "n_by_seed": [int(mode["n"]) for mode in modes],
                "mean_negative_log_loss": statistics.fmean(
                    float(mode["negative_log_loss"]) for mode in modes
                ),
                "mean_multiclass_brier": statistics.fmean(
                    float(mode["multiclass_brier"]) for mode in modes
                ),
                "mean_top1_percent": statistics.fmean(
                    float(mode["top1_percent"]) for mode in modes
                ),
            }
        )
    if len(records) != 10:
        raise SystemExit(f"expected 10 scales, found {len(records)}")
    selected = min(records, key=lambda row: (row["mean_negative_log_loss"], -row["depth_scale"]))
    args.out.write_text(
        json.dumps(
            {
                "status": "D0_development_selection_only",
                "analysis": "monotone_depth_gated_hybrid",
                "selection_rule": "minimum mean NLL; larger scale exact tie-break",
                "configuration_count": len(records),
                "selected": selected,
                "configurations": records,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
