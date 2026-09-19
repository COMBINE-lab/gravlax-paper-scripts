#!/usr/bin/env python3
"""Apply the locked joint D0/D4 selection rule for hybrid promotion."""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


METRICS = ("negative_log_loss", "multiclass_brier", "top1_percent")
DATASETS = ("D0", "D4")
SEEDS = (7, 17, 29)


def load_modes(path: Path) -> dict[str, dict[str, object]]:
    payload = json.loads(path.read_text())
    return {str(mode["name"]): mode for mode in payload["modes"]}


def weighted(rows: list[dict[str, object]], metric: str) -> float:
    total = sum(int(row["n"]) for row in rows)
    return sum(int(row["n"]) * float(row[metric]) for row in rows) / total


def aggregate(rows: list[dict[str, object]]) -> dict[str, object]:
    if len(rows) != len(SEEDS):
        raise SystemExit(f"expected {len(SEEDS)} seed rows, found {len(rows)}")
    return {
        "targets": sum(int(row["n"]) for row in rows),
        **{metric: weighted(rows, metric) for metric in METRICS},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    baselines: dict[str, dict[str, list[dict[str, object]]]] = {
        dataset: defaultdict(list) for dataset in DATASETS
    }
    grid: dict[tuple[float, float, str], list[dict[str, object]]] = defaultdict(list)
    with (args.root / "runs.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            dataset = row["dataset"]
            modes = load_modes(Path(row["out"]) / "metrics.json")
            if row["run_type"] == "baseline":
                for name in ("pooled", "dirichlet-proxy", "depth-hybrid"):
                    if name not in modes:
                        raise SystemExit(f"missing {name} in {row['out']}")
                    baselines[dataset][name].append(modes[name])
            else:
                if tuple(modes) != ("depth-hybrid",):
                    raise SystemExit(f"expected hybrid-only metrics in {row['out']}")
                key = (float(row["depth_scale"]), float(row["depth_power"]), dataset)
                grid[key].append(modes["depth-hybrid"])

    reference = {
        dataset: {
            name: aggregate(baselines[dataset][name])
            for name in ("pooled", "dirichlet-proxy")
        }
        for dataset in DATASETS
    }
    configurations = []
    for scale in (2, 4, 8, 16, 32, 64, 128, 256):
        for power in (0.5, 1, 2, 4, 8):
            by_dataset = {
                dataset: aggregate(grid[(float(scale), float(power), dataset)])
                for dataset in DATASETS
            }
            criteria: dict[str, bool] = {}
            for dataset in DATASETS:
                candidate = by_dataset[dataset]
                refs = reference[dataset]
                criteria[f"{dataset}_NLL_below_pooled"] = (
                    float(candidate["negative_log_loss"])
                    < float(refs["pooled"]["negative_log_loss"])
                )
                criteria[f"{dataset}_NLL_below_proxy"] = (
                    float(candidate["negative_log_loss"])
                    < float(refs["dirichlet-proxy"]["negative_log_loss"])
                )
                criteria[f"{dataset}_Brier_at_most_pooled"] = (
                    float(candidate["multiclass_brier"])
                    <= float(refs["pooled"]["multiclass_brier"])
                )
                criteria[f"{dataset}_Brier_at_most_proxy"] = (
                    float(candidate["multiclass_brier"])
                    <= float(refs["dirichlet-proxy"]["multiclass_brier"])
                )
                best_top1 = max(
                    float(refs["pooled"]["top1_percent"]),
                    float(refs["dirichlet-proxy"]["top1_percent"]),
                )
                criteria[f"{dataset}_top1_drop_at_most_0_10pp"] = (
                    best_top1 - float(candidate["top1_percent"]) <= 0.10
                )
            targets = sum(int(by_dataset[d]["targets"]) for d in DATASETS)
            joint_nll = sum(
                int(by_dataset[d]["targets"])
                * float(by_dataset[d]["negative_log_loss"])
                for d in DATASETS
            ) / targets
            configurations.append(
                {
                    "depth_scale": scale,
                    "depth_power": power,
                    "joint_target_weighted_negative_log_loss": joint_nll,
                    "by_dataset": by_dataset,
                    "criteria": criteria,
                    "eligible": all(criteria.values()),
                }
            )

    eligible = [row for row in configurations if row["eligible"]]
    selected = min(
        eligible,
        key=lambda row: (
            row["joint_target_weighted_negative_log_loss"],
            row["depth_power"],
            row["depth_scale"],
        ),
        default=None,
    )
    args.out.write_text(
        json.dumps(
            {
                "status": "ELIGIBLE_MODEL_SELECTED" if selected else "NO_ELIGIBLE_MODEL",
                "selection_rule": (
                    "joint minimum target-weighted NLL among candidates satisfying every "
                    "per-dataset pooled/proxy NLL, Brier, and top-1 criterion"
                ),
                "reference": reference,
                "eligible_configuration_count": len(eligible),
                "selected": selected,
                "configurations": configurations,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()

