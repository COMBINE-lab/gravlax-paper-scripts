#!/usr/bin/env python3
"""Select the locked D0 posterior-mean Dirichlet-proxy configuration."""
from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path


def load_mode(path: Path, expected: str) -> dict[str, object]:
    payload = json.loads(path.read_text())
    modes = payload.get("modes", [])
    if len(modes) != 1 or modes[0].get("name") != expected:
        raise SystemExit(f"expected one {expected} mode in {path}")
    return modes[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    grouped: dict[tuple[float, float], list[dict[str, object]]] = defaultdict(list)
    with (args.root / "runs.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["model"] != "dirichlet-proxy" or row["control"] != "real":
                continue
            key = (float(row["cell_prior"]), float(row["group_prior"]))
            grouped[key].append(load_mode(Path(row["out"]) / "metrics.json", "dirichlet-proxy"))

    records = []
    for (cell_prior, group_prior), modes in sorted(grouped.items()):
        if len(modes) != 3:
            raise SystemExit(f"configuration {(cell_prior, group_prior)} has {len(modes)} seeds")
        records.append(
            {
                "cell_prior": cell_prior,
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
    if len(records) != 25:
        raise SystemExit(f"expected 25 configurations, found {len(records)}")

    selected = min(
        records,
        key=lambda row: (
            row["mean_negative_log_loss"],
            row["cell_prior"],
            row["group_prior"],
        ),
    )
    result = {
        "status": "D0_development_selection_only",
        "analysis": "posterior_mean_dirichlet_proxy",
        "selection_dataset": "D0",
        "selection_rule": "minimum mean NLL; lower cell prior, then group prior tie-break",
        "configuration_count": len(records),
        "selected": selected,
        "configurations": records,
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
