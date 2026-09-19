#!/usr/bin/env python3
"""Summarize the registered hierarchical-EM gate and apply its frozen criteria."""
from __future__ import annotations

import argparse
import json
import re
import statistics
from pathlib import Path


def parse_time(path: Path) -> dict[str, float]:
    text = path.read_text()
    elapsed = re.search(r"Elapsed \(wall clock\) time .*: (.+)", text).group(1).strip()
    parts = [float(value) for value in elapsed.split(":")]
    wall = sum(value * 60**index for index, value in enumerate(reversed(parts)))
    rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", text).group(1))
    return {"wall_seconds": wall, "max_rss_kb": rss}


def read_run(path: Path) -> dict[str, object]:
    data = json.loads((path / "metrics.json").read_text())
    return {
        "seed": int(data["seed"]),
        "time": parse_time(path / "time.txt"),
        "modes": {item["name"]: item for item in data["modes"]},
    }


def read_protocol(path: Path) -> dict[str, str]:
    rows = [line.rstrip("\n").split("\t", 1) for line in path.open()]
    if not rows or rows[0] != ["key", "value"]:
        raise SystemExit(f"malformed protocol: {path}")
    return dict(rows[1:])


def aggregate(runs: list[dict[str, object]], mode: str) -> dict[str, float]:
    rows = [run["modes"][mode] for run in runs]
    total_n = sum(int(row["n"]) for row in rows)
    return {
        "n_across_seeds": total_n,
        "top1_percent": sum(float(row["top1_percent"]) * int(row["n"]) for row in rows) / total_n,
        "expected_accuracy_percent": sum(
            float(row["expected_accuracy_percent"]) * int(row["n"]) for row in rows
        ) / total_n,
        "negative_log_loss": sum(
            float(row["negative_log_loss"]) * int(row["n"]) for row in rows
        ) / total_n,
        "multiclass_brier": sum(
            float(row["multiclass_brier"]) * int(row["n"]) for row in rows
        ) / total_n,
    }


def load_control(root: Path, dataset: str, control: str) -> list[dict[str, object]]:
    return [read_run(root / "selected" / dataset / control / f"seed-{seed}") for seed in (7, 17, 29)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gate_dir", type=Path)
    parser.add_argument("--packed-d1", type=Path, required=True)
    parser.add_argument("--groups-json", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    selection = json.loads((args.gate_dir / "selection.json").read_text())
    groups = json.loads(args.groups_json.read_text())
    datasets: dict[str, object] = {}
    collapse_max = 0.0
    all_target_counts_equal = True
    for dataset in ("D0", "D1", "D2_prime"):
        controls = {name: load_control(args.gate_dir, dataset, name) for name in ("real", "shuffled", "collapse")}
        target_counts = {
            f"{control}:seed-{run['seed']}": run["modes"]["pooled"]["n"]
            for control, runs in controls.items()
            for run in runs
        }
        all_target_counts_equal &= all(
            len(
                {
                    next(run for run in controls[control] if run["seed"] == seed)["modes"][
                        "pooled"
                    ]["n"]
                    for control in controls
                }
            )
            == 1
            for seed in (7, 17, 29)
        )
        summaries = {
            control: {mode: aggregate(runs, mode) for mode in ("pooled", "blend", "group", "hierarchical")}
            for control, runs in controls.items()
        }
        for left, right in (("negative_log_loss", "negative_log_loss"), ("multiclass_brier", "multiclass_brier"), ("top1_percent", "top1_percent"), ("expected_accuracy_percent", "expected_accuracy_percent")):
            collapse_max = max(
                collapse_max,
                abs(summaries["collapse"]["group"][left] - summaries["collapse"]["pooled"][right]),
            )
        baseline_name = min(
            ("pooled", "blend"),
            key=lambda mode: summaries["real"][mode]["negative_log_loss"],
        )
        wall = [run["time"]["wall_seconds"] for run in controls["real"]]
        rss = [run["time"]["max_rss_kb"] for run in controls["real"]]
        datasets[dataset] = {
            "target_counts": target_counts,
            "controls": summaries,
            "best_nonhierarchical_mode": baseline_name,
            "best_nonhierarchical": summaries["real"][baseline_name],
            "real_hierarchical": summaries["real"]["hierarchical"],
            "shuffled_hierarchical": summaries["shuffled"]["hierarchical"],
            "timing_real": {
                "wall_seconds": wall,
                "wall_median_seconds": statistics.median(wall),
                "max_rss_kb": rss,
                "max_rss_kb_max": max(rss),
            },
        }

    confirmations = [datasets["D1"], datasets["D2_prime"]]
    baseline_loss_sum = sum(
        item["best_nonhierarchical"]["negative_log_loss"] * item["best_nonhierarchical"]["n_across_seeds"]
        for item in confirmations
    )
    hierarchy_loss_sum = sum(
        item["real_hierarchical"]["negative_log_loss"] * item["real_hierarchical"]["n_across_seeds"]
        for item in confirmations
    )
    confirmation_n = sum(item["real_hierarchical"]["n_across_seeds"] for item in confirmations)
    combined_baseline = baseline_loss_sum / confirmation_n
    combined_hierarchy = hierarchy_loss_sum / confirmation_n
    combined_improvement = (combined_baseline - combined_hierarchy) / combined_baseline

    per_dataset_top1_ok = all(
        item["real_hierarchical"]["top1_percent"]
        >= item["best_nonhierarchical"]["top1_percent"] - 0.5
        for item in confirmations
    )
    d2 = datasets["D2_prime"]
    d2_top1_gain = d2["real_hierarchical"]["top1_percent"] - d2["best_nonhierarchical"]["top1_percent"]
    d2_loss_gain = (
        d2["best_nonhierarchical"]["negative_log_loss"] - d2["real_hierarchical"]["negative_log_loss"]
    ) / d2["best_nonhierarchical"]["negative_log_loss"]
    real_gain = sum(
        (item["best_nonhierarchical"]["negative_log_loss"] - item["real_hierarchical"]["negative_log_loss"])
        * item["real_hierarchical"]["n_across_seeds"]
        for item in confirmations
    )
    beyond_shuffle = sum(
        (item["shuffled_hierarchical"]["negative_log_loss"] - item["real_hierarchical"]["negative_log_loss"])
        * item["real_hierarchical"]["n_across_seeds"]
        for item in confirmations
    )
    real_fraction_beyond_shuffle = beyond_shuffle / real_gain if real_gain > 0 else None

    packed = json.loads(args.packed_d1.read_text())
    packed_wall = float(packed["arms"]["masked"]["wall_median_seconds"])
    d1_wall_ratio = datasets["D1"]["timing_real"]["wall_median_seconds"] / packed_wall
    checks = {
        "target_counts_identical_across_controls": all_target_counts_equal,
        "combined_confirmation_log_loss_improvement_at_least_5pct": combined_improvement >= 0.05,
        "per_confirmation_top1_drop_at_most_0_5pp": per_dataset_top1_ok,
        "d2_sparse_alternative": d2_top1_gain >= 1.0 or d2_loss_gain >= 0.10,
        "real_groups_explain_half_gain_beyond_shuffle": (
            real_fraction_beyond_shuffle is not None and real_fraction_beyond_shuffle >= 0.50
        ),
        "one_group_collapse_matches_pooled": collapse_max <= 1e-10,
        "d1_wall_at_most_1_5x_packed": d1_wall_ratio <= 1.5,
        "d1_peak_rss_below_6_gib": datasets["D1"]["timing_real"]["max_rss_kb_max"] < 6 * 1024 * 1024,
    }
    result = {
        "status": "complete",
        "analysis": "registered_hierarchical_em_gate",
        "protocol": read_protocol(args.gate_dir / "protocol.tsv"),
        "selection": selection,
        "group_construction": groups,
        "datasets": datasets,
        "confirmation": {
            "combined_baseline_negative_log_loss": combined_baseline,
            "combined_hierarchical_negative_log_loss": combined_hierarchy,
            "relative_log_loss_improvement": combined_improvement,
            "d2_top1_gain_percentage_points": d2_top1_gain,
            "d2_relative_log_loss_improvement": d2_loss_gain,
            "real_fraction_of_gain_beyond_shuffle": real_fraction_beyond_shuffle,
            "collapse_max_metric_absolute_difference": collapse_max,
            "d1_wall_ratio_vs_packed": d1_wall_ratio,
        },
        "gate": {"checks": checks, "verdict": "PASS" if all(checks.values()) else "STOP"},
    }
    rendered = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.out:
        args.out.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
