#!/usr/bin/env python3
"""Summarize one locked depth-hybrid confirmation dataset."""
from __future__ import annotations

import argparse
import json
import re
import statistics
from pathlib import Path

import yaml


METRICS = (
    "negative_log_loss",
    "multiclass_brier",
    "top1_percent",
    "expected_accuracy_percent",
)
DEPTH_NAMES = ("0-1", "1-4", "4-16", "16+")


def load_mode(path: Path, expected: str) -> dict[str, object]:
    payload = json.loads(path.read_text())
    modes = payload.get("modes", [])
    if len(modes) != 1 or modes[0].get("name") != expected:
        raise SystemExit(f"expected one {expected} mode in {path}")
    return {"seed": int(payload["seed"]), **modes[0]}


def summarize(root: Path, expected: str) -> dict[str, object]:
    rows = sorted(
        (load_mode(path, expected) for path in root.glob("seed-*/metrics.json")),
        key=lambda row: int(row["seed"]),
    )
    if [int(row["seed"]) for row in rows] != [7, 17, 29]:
        raise SystemExit(f"expected seeds 7, 17, 29 under {root}")
    depth = []
    for name in DEPTH_NAMES:
        strata = []
        for row in rows:
            by_name = {stratum["name"]: stratum for stratum in row["evidence_depth_strata"]}
            if tuple(by_name) != DEPTH_NAMES:
                raise SystemExit(f"unexpected depth strata under {root}")
            strata.append(by_name[name])
        depth.append(
            {
                "name": name,
                "n_by_seed": [int(item["n"]) for item in strata],
                **{
                    f"mean_{metric}": statistics.fmean(float(item[metric]) for item in strata)
                    for metric in METRICS
                },
            }
        )
    return {
        "means": {
            metric: statistics.fmean(float(row[metric]) for row in rows)
            for metric in METRICS
        },
        "by_seed": [
            {"seed": row["seed"], "n": row["n"], **{metric: row[metric] for metric in METRICS}}
            for row in rows
        ],
        "evidence_depth_strata": depth,
    }


def parse_protocol(path: Path) -> dict[str, str]:
    rows = [line.rstrip("\n").split("\t", 1) for line in path.open()]
    return dict(rows[1:])


def elapsed_seconds(value: str) -> float:
    fields = [float(field) for field in value.split(":")]
    return sum(field * 60**power for power, field in enumerate(reversed(fields)))


def timings(root: Path) -> dict[str, object]:
    walls, rss = [], []
    for path in root.glob("**/time.txt"):
        text = path.read_text()
        wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)", text)
        memory = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
        if wall is None or memory is None:
            raise SystemExit(f"incomplete timing: {path}")
        walls.append(elapsed_seconds(wall.group(1)))
        rss.append(int(memory.group(1)))
    if len(walls) != 15:
        raise SystemExit(f"expected 15 timed runs, found {len(walls)}")
    return {
        "runs": len(walls),
        "wall_seconds": {
            "minimum": min(walls),
            "median": statistics.median(walls),
            "maximum": max(walls),
        },
        "peak_rss_kib": {
            "minimum": min(rss),
            "median": statistics.median(rss),
            "maximum": max(rss),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest = yaml.safe_load(args.manifest.read_text())
    protocol = parse_protocol(args.root / "protocol.tsv")
    dataset = protocol["dataset"]
    acceptance = manifest["acceptance"]["D4" if dataset == "D4" else "D3_replication"]

    hybrid_real = summarize(args.root / "depth-hybrid/real", "depth-hybrid")
    hybrid_shuffled = summarize(args.root / "depth-hybrid/shuffled", "depth-hybrid")
    proxy_real = summarize(args.root / "dirichlet-proxy/real", "dirichlet-proxy")
    proxy_shuffled = summarize(args.root / "dirichlet-proxy/shuffled", "dirichlet-proxy")
    fixed_real = summarize(args.root / "fixed-convex/real", "convex")
    arms = (hybrid_real, hybrid_shuffled, proxy_real, proxy_shuffled, fixed_real)
    populations = [[row["n"] for row in arm["by_seed"]] for arm in arms]
    if any(population != populations[0] for population in populations[1:]):
        raise SystemExit("factorial arms are not paired on identical populations")

    hybrid_nll = float(hybrid_real["means"]["negative_log_loss"])
    proxy_nll = float(proxy_real["means"]["negative_log_loss"])
    shuffled_nll = float(hybrid_shuffled["means"]["negative_log_loss"])
    model_gain = (proxy_nll - hybrid_nll) / proxy_nll
    group_gain = (shuffled_nll - hybrid_nll) / shuffled_nll
    model_wins = [
        float(hybrid["negative_log_loss"]) < float(proxy["negative_log_loss"])
        for hybrid, proxy in zip(hybrid_real["by_seed"], proxy_real["by_seed"], strict=True)
    ]
    group_wins = [
        float(real["negative_log_loss"]) < float(shuffled["negative_log_loss"])
        for real, shuffled in zip(
            hybrid_real["by_seed"], hybrid_shuffled["by_seed"], strict=True
        )
    ]
    depth_comparison = []
    for index, hybrid in enumerate(hybrid_real["evidence_depth_strata"]):
        reference = (
            fixed_real["evidence_depth_strata"][index]
            if index < 3
            else proxy_real["evidence_depth_strata"][index]
        )
        hybrid_depth = float(hybrid["mean_negative_log_loss"])
        reference_depth = float(reference["mean_negative_log_loss"])
        depth_comparison.append(
            {
                "name": hybrid["name"],
                "n_by_seed": hybrid["n_by_seed"],
                "reference": "fixed-convex" if index < 3 else "Dirichlet-proxy",
                "hybrid_mean_negative_log_loss": hybrid_depth,
                "reference_mean_negative_log_loss": reference_depth,
                "relative_nll_change_vs_reference":
                    (hybrid_depth - reference_depth) / reference_depth,
                "absolute_brier_change_vs_reference":
                    float(hybrid["mean_multiclass_brier"])
                    - float(reference["mean_multiclass_brier"]),
                "top1_percentage_point_change_vs_reference":
                    float(hybrid["mean_top1_percent"])
                    - float(reference["mean_top1_percent"]),
            }
        )
    proxy_depth_audit = []
    for hybrid, proxy in zip(
        hybrid_real["evidence_depth_strata"],
        proxy_real["evidence_depth_strata"],
        strict=True,
    ):
        proxy_depth_audit.append(
            {
                "name": hybrid["name"],
                "n_by_seed": hybrid["n_by_seed"],
                "absolute_nll_change_hybrid_minus_proxy":
                    float(hybrid["mean_negative_log_loss"])
                    - float(proxy["mean_negative_log_loss"]),
                "absolute_brier_change_hybrid_minus_proxy":
                    float(hybrid["mean_multiclass_brier"])
                    - float(proxy["mean_multiclass_brier"]),
                "top1_percentage_point_change_hybrid_minus_proxy":
                    float(hybrid["mean_top1_percent"])
                    - float(proxy["mean_top1_percent"]),
            }
        )
    performance = timings(args.root)
    if dataset == "D4":
        criteria = {
            "model_form_hybrid_beats_proxy_on_all_seeds": all(model_wins),
            "model_form_mean_relative_nll_improvement_at_least_0_2_percent":
                model_gain >= float(acceptance["model_form_mean_relative_nll_improvement_min"]),
            "group_information_real_beats_shuffled_on_all_seeds": all(group_wins),
            "group_information_mean_relative_nll_improvement_at_least_0_2_percent":
                group_gain >= float(acceptance["group_information_mean_relative_nll_improvement_min"]),
        }
    else:
        criteria = {
            "model_form_mean_relative_nll_improvement_positive":
                model_gain > float(acceptance["model_form_mean_relative_nll_improvement_min"]),
            "model_form_hybrid_beats_proxy_on_at_least_2_seeds":
                sum(model_wins) >= int(acceptance["model_form_hybrid_beats_proxy_seeds_min"]),
            "group_information_mean_relative_nll_improvement_positive":
                group_gain > float(acceptance["group_information_mean_relative_nll_improvement_min"]),
            "group_information_real_beats_shuffled_on_at_least_2_seeds":
                sum(group_wins) >= int(acceptance["group_information_real_beats_shuffled_seeds_min"]),
        }
    criteria.update(
        {
            "model_form_mean_brier_does_not_increase":
                float(hybrid_real["means"]["multiclass_brier"])
                <= float(proxy_real["means"]["multiclass_brier"]),
            "group_information_mean_brier_does_not_increase":
                float(hybrid_real["means"]["multiclass_brier"])
                <= float(hybrid_shuffled["means"]["multiclass_brier"]),
            "top1_drop_vs_proxy_at_most_0_10_percentage_point":
                float(proxy_real["means"]["top1_percent"])
                - float(hybrid_real["means"]["top1_percent"])
                <= float(acceptance["top1_drop_vs_proxy_percentage_points_max"]),
            "each_shallow_stratum_regression_vs_fixed_at_most_0_25_percent": all(
                float(row["relative_nll_change_vs_reference"])
                <= float(acceptance["each_shallow_stratum_relative_nll_regression_vs_fixed_max"])
                for row in depth_comparison[:3]
            ),
            "high_stratum_regression_vs_proxy_at_most_0_25_percent":
                float(depth_comparison[3]["relative_nll_change_vs_reference"])
                <= float(acceptance["high_stratum_relative_nll_regression_vs_proxy_max"]),
            "peak_rss_below_4_GiB":
                int(performance["peak_rss_kib"]["maximum"])
                < float(acceptance["peak_rss_gib_max"]) * 1024 * 1024,
        }
    )
    if all(criteria.values()):
        verdict = "PASS"
    elif model_gain > 0 and group_gain > 0:
        verdict = "MARGINAL"
    else:
        verdict = "STOP"

    result = {
        "status": "untouched_confirmation_complete" if dataset == "D4" else "held_back_replication_complete",
        "dataset": dataset,
        "protocol": protocol,
        "arms": {
            "hybrid_real": hybrid_real,
            "hybrid_shuffled": hybrid_shuffled,
            "Dirichlet_proxy_real": proxy_real,
            "Dirichlet_proxy_shuffled": proxy_shuffled,
            "fixed_convex_real": fixed_real,
        },
        "factorial_contrasts": {
            "model_form_relative_nll_improvement_hybrid_vs_proxy_real": model_gain,
            "model_form_hybrid_better_by_seed": model_wins,
            "group_information_relative_nll_improvement_real_vs_shuffled_hybrid": group_gain,
            "group_information_real_better_by_seed": group_wins,
            "model_form_absolute_brier_change_hybrid_minus_proxy":
                float(hybrid_real["means"]["multiclass_brier"])
                - float(proxy_real["means"]["multiclass_brier"]),
            "model_form_brier_change_by_seed": [
                float(hybrid["multiclass_brier"]) - float(proxy["multiclass_brier"])
                for hybrid, proxy in zip(
                    hybrid_real["by_seed"], proxy_real["by_seed"], strict=True
                )
            ],
            "model_form_top1_percentage_point_change_hybrid_minus_proxy":
                float(hybrid_real["means"]["top1_percent"])
                - float(proxy_real["means"]["top1_percent"]),
            "group_information_absolute_brier_change_real_minus_shuffled":
                float(hybrid_real["means"]["multiclass_brier"])
                - float(hybrid_shuffled["means"]["multiclass_brier"]),
            "group_information_top1_percentage_point_change_real_minus_shuffled":
                float(hybrid_real["means"]["top1_percent"])
                - float(hybrid_shuffled["means"]["top1_percent"]),
            "proxy_group_information_relative_nll_improvement_real_vs_shuffled":
                (
                    float(proxy_shuffled["means"]["negative_log_loss"])
                    - float(proxy_real["means"]["negative_log_loss"])
                ) / float(proxy_shuffled["means"]["negative_log_loss"]),
            "difference_in_differences_NLL":
                (
                    float(proxy_real["means"]["negative_log_loss"])
                    - float(hybrid_real["means"]["negative_log_loss"])
                )
                - (
                    float(proxy_shuffled["means"]["negative_log_loss"])
                    - float(hybrid_shuffled["means"]["negative_log_loss"])
                ),
        },
        "evidence_depth_comparison": depth_comparison,
        "model_form_depth_residual_audit": proxy_depth_audit,
        "acceptance": criteria,
        "verdict": verdict,
        "next_step": (
            "open_D3_replication" if dataset == "D4" and verdict == "PASS"
            else "evaluate_replication_and_close_EM_branch" if dataset == "D3"
            else "close_EM_branch_without_D3"
        ),
        "full_hierarchical_Dirichlet": {
            "licensed": False,
            "structured_depth_residual_observed":
                any(float(row["absolute_brier_change_hybrid_minus_proxy"]) < 0 for row in proxy_depth_audit)
                and any(float(row["absolute_brier_change_hybrid_minus_proxy"]) > 0 for row in proxy_depth_audit),
            "replicated_real_group_signal_available": dataset == "D3" and verdict == "PASS",
            "decision":
                "not_licensed_D4_did_not_open_replication"
                if dataset == "D4" and verdict != "PASS"
                else "not_licensed_by_this_single_dataset_summary",
        },
        "production_emission": "pooled",
        "performance": performance,
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
