#!/usr/bin/env python3
"""Select the locked D0 convex-EM configuration and its tuned no-group comparator."""
from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path


def load_mode(path: Path) -> dict[str, float]:
    payload = json.loads(path.read_text())
    modes = payload.get("modes", [])
    if len(modes) != 1 or modes[0].get("name") != "convex":
        raise SystemExit(f"expected one convex mode in {path}")
    return modes[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    grouped: dict[tuple[float, float, float], list[dict[str, float]]] = defaultdict(list)
    with (args.root / "runs.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["control"] != "real":
                continue
            key = (
                float(row["cell_weight"]),
                float(row["group_weight"]),
                float(row["group_prior"]),
            )
            grouped[key].append(load_mode(Path(row["out"]) / "metrics.json"))

    expected_seeds = 3
    records = []
    for (cell_weight, group_weight, group_prior), modes in sorted(grouped.items()):
        if len(modes) != expected_seeds:
            raise SystemExit(
                f"configuration {(cell_weight, group_weight, group_prior)} has {len(modes)} seeds"
            )
        records.append(
            {
                "cell_weight": cell_weight,
                "group_weight": group_weight,
                "global_weight": round(1.0 - cell_weight - group_weight, 12),
                "group_prior": group_prior,
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

    if len(records) != 68:
        raise SystemExit(f"expected 68 configurations, found {len(records)}")

    def rank(record: dict[str, float]) -> tuple[float, float, float, float]:
        return (
            record["mean_negative_log_loss"],
            record["group_weight"],
            record["cell_weight"],
            record["group_prior"],
        )

    selected = min(records, key=rank)
    no_group = min((record for record in records if record["group_weight"] == 0.0), key=rank)
    group_increment = (
        no_group["mean_negative_log_loss"] - selected["mean_negative_log_loss"]
    ) / no_group["mean_negative_log_loss"]
    result = {
        "status": "development_selection_only",
        "analysis": "candidate_normalized_convex_em_D0_grid",
        "selection_dataset": "D0",
        "selection_rule": "minimum mean NLL; lower group, cell, then prior tie-break",
        "configuration_count": len(records),
        "selected": selected,
        "tuned_no_group": no_group,
        "relative_nll_group_increment_over_tuned_no_group": group_increment,
        "configurations": records,
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
