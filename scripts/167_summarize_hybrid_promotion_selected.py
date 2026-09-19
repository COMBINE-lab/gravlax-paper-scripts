#!/usr/bin/env python3
"""Summarize paired D0/D4 diagnostics for the selected hybrid candidate."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


DATASETS = ("D0", "D4")
REFERENCES = ("pooled", "dirichlet-proxy")


def weighted(values: list[tuple[int, float]]) -> float:
    total = sum(n for n, _ in values)
    return sum(n * value for n, value in values) / total


def comparison(payload: dict[str, object], reference: str) -> dict[str, object]:
    for item in payload["paired_comparisons"]:
        if item["candidate"] == "depth-hybrid" and item["reference"] == reference:
            return item
    raise SystemExit(f"missing depth-hybrid versus {reference} paired comparison")


def aggregate_comparison(rows: list[dict[str, object]], reference: str) -> dict[str, object]:
    pairs = [comparison(row, reference) for row in rows]
    total = sum(int(pair["n"]) for pair in pairs)
    nll = weighted([(int(pair["n"]), float(pair["mean_negative_log_loss_difference"])) for pair in pairs])
    brier = weighted([(int(pair["n"]), float(pair["mean_brier_difference"])) for pair in pairs])
    # Seeds reuse cells and partially overlapping masked targets. Treat their clustered estimates
    # as perfectly positively correlated: sum(weight * SE), a conservative bound versus an
    # independence assumption.
    nll_se = sum(int(pair["n"]) * float(pair["negative_log_loss_clustered_se"]) for pair in pairs) / total
    brier_se = sum(int(pair["n"]) * float(pair["brier_clustered_se"]) for pair in pairs) / total
    return {
        "candidate": "depth-hybrid",
        "reference": reference,
        "n": total,
        "seed_combination": "target_weighted_mean_with_perfect_positive_SE_correlation",
        "mean_negative_log_loss_difference": nll,
        "negative_log_loss_conservative_se": nll_se,
        "negative_log_loss_ci95": [nll - 1.96 * nll_se, nll + 1.96 * nll_se],
        "mean_brier_difference": brier,
        "brier_conservative_se": brier_se,
        "brier_ci95": [brier - 1.96 * brier_se, brier + 1.96 * brier_se],
        "top1_percentage_point_difference": weighted(
            [(int(pair["n"]), float(pair["top1_percentage_point_difference"])) for pair in pairs]
        ),
        "by_seed": pairs,
    }


def aggregate_modes(rows: list[dict[str, object]]) -> dict[str, dict[str, float]]:
    result = {}
    for name in ("pooled", "dirichlet-proxy", "depth-hybrid"):
        modes = [next(mode for mode in row["modes"] if mode["name"] == name) for row in rows]
        result[name] = {
            metric: weighted([(int(mode["n"]), float(mode[metric])) for mode in modes])
            for metric in ("negative_log_loss", "multiclass_brier", "top1_percent")
        }
    return result


def elapsed_seconds(value: str) -> float:
    fields = [float(field) for field in value.split(":")]
    return sum(field * 60**power for power, field in enumerate(reversed(fields)))


def performance(root: Path) -> dict[str, float]:
    walls, rss = [], []
    for path in root.glob("D*/seed-*/time.txt"):
        text = path.read_text()
        wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)", text)
        memory = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
        if wall is None or memory is None:
            raise SystemExit(f"incomplete timing: {path}")
        walls.append(elapsed_seconds(wall.group(1)))
        rss.append(int(memory.group(1)))
    if len(walls) != 6:
        raise SystemExit(f"expected six timed runs, found {len(walls)}")
    return {"maximum_wall_seconds": max(walls), "maximum_peak_rss_kib": max(rss)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    selection = json.loads(args.selection.read_text())
    if selection["status"] != "ELIGIBLE_MODEL_SELECTED":
        raise SystemExit("selected audit requires an eligible development model")

    datasets = {}
    criteria = {}
    for dataset in DATASETS:
        paths = sorted((args.root / dataset).glob("seed-*/metrics.json"))
        rows = [json.loads(path.read_text()) for path in paths]
        rows.sort(key=lambda row: (7, 17, 29).index(int(row["seed"])))
        if [int(row["seed"]) for row in rows] != [7, 17, 29]:
            raise SystemExit(f"expected seeds 7, 17, 29 for {dataset}")
        comparisons = {ref: aggregate_comparison(rows, ref) for ref in REFERENCES}
        modes = aggregate_modes(rows)
        datasets[dataset] = {"modes": modes, "paired_comparisons": comparisons}
        for ref in REFERENCES:
            item = comparisons[ref]
            criteria[f"{dataset}_NLL_upper_CI_below_zero_vs_{ref}"] = (
                float(item["negative_log_loss_ci95"][1]) < 0
            )
            criteria[f"{dataset}_Brier_upper_CI_at_most_zero_vs_{ref}"] = (
                float(item["brier_ci95"][1]) <= 0
            )
            criteria[f"{dataset}_top1_drop_at_most_0_10pp_vs_{ref}"] = (
                float(item["top1_percentage_point_difference"]) >= -0.10
            )
    perf = performance(args.root)
    criteria["peak_RSS_below_4_GiB"] = perf["maximum_peak_rss_kib"] < 4 * 1024 * 1024
    result = {
        "status": "selected_D0_D4_development_audit_complete",
        "selected": selection["selected"],
        "datasets": datasets,
        "acceptance": criteria,
        "verdict": "DEVELOPMENT_PASS" if all(criteria.values()) else "DEVELOPMENT_FAIL",
        "next_step": "freeze_and_open_D3" if all(criteria.values()) else "retain_pooled_default",
        "performance": perf,
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
