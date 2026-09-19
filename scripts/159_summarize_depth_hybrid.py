#!/usr/bin/env python3
"""Evaluate the locked monotone depth-hybrid D0 gate."""
from __future__ import annotations

import argparse
import json
import re
import statistics
from pathlib import Path


METRICS = ("negative_log_loss", "multiclass_brier", "top1_percent")
DEPTH_NAMES = ("0-1", "1-4", "4-16", "16+")


def load_mode(path: Path, expected: str) -> dict[str, object]:
    payload = json.loads(path.read_text())
    modes = payload.get("modes", [])
    if len(modes) != 1 or modes[0].get("name") != expected:
        raise SystemExit(f"expected one {expected} mode in {path}")
    return {"seed": int(payload["seed"]), **modes[0]}


def summarize(paths: list[Path], expected: str) -> dict[str, object]:
    rows = sorted((load_mode(path, expected) for path in paths), key=lambda row: int(row["seed"]))
    if [row["seed"] for row in rows] != [7, 17, 29]:
        raise SystemExit(f"expected {expected} seeds 7, 17, and 29")
    depth = []
    for name in DEPTH_NAMES:
        strata = []
        for row in rows:
            by_name = {item["name"]: item for item in row["evidence_depth_strata"]}
            if tuple(by_name) != DEPTH_NAMES:
                raise SystemExit(f"unexpected depth strata in {expected}, seed {row['seed']}")
            strata.append(by_name[name])
        depth.append(
            {
                "name": name,
                "n_by_seed": [int(item["n"]) for item in strata],
                "by_seed": [
                    {
                        "seed": row["seed"],
                        "n": item["n"],
                        **{metric: item[metric] for metric in METRICS},
                    }
                    for row, item in zip(rows, strata, strict=True)
                ],
                **{
                    f"mean_{metric}": statistics.fmean(float(item[metric]) for item in strata)
                    for metric in METRICS
                },
            }
        )
    return {
        "means": {metric: statistics.fmean(float(row[metric]) for row in rows) for metric in METRICS},
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


def timing(root: Path) -> dict[str, object]:
    walls, rss = [], []
    for path in root.glob("**/time.txt"):
        text = path.read_text()
        wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)", text)
        memory = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
        if wall is None or memory is None:
            raise SystemExit(f"incomplete timing file: {path}")
        walls.append(elapsed_seconds(wall.group(1)))
        rss.append(int(memory.group(1)))
    if len(walls) != 45:
        raise SystemExit(f"expected 45 timings, found {len(walls)}")
    return {
        "runs": len(walls),
        "wall_seconds": {"minimum": min(walls), "median": statistics.median(walls), "maximum": max(walls)},
        "peak_rss_kib": {"minimum": min(rss), "median": statistics.median(rss), "maximum": max(rss)},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    selection = json.loads((args.root / "selection.json").read_text())
    scale = float(selection["selected"]["depth_scale"])
    scale_label = f"{scale:g}"
    hybrid_real = summarize(sorted((args.root / f"grid/scale-{scale_label}").glob("seed-*/metrics.json")), "depth-hybrid")
    hybrid_shuffled = summarize(sorted((args.root / "selected-shuffled").glob("seed-*/metrics.json")), "depth-hybrid")
    fixed_real = summarize(sorted((args.root / "constituents/fixed-convex/real").glob("seed-*/metrics.json")), "convex")
    fixed_shuffled = summarize(sorted((args.root / "constituents/fixed-convex/shuffled").glob("seed-*/metrics.json")), "convex")
    proxy_real = summarize(sorted((args.root / "constituents/dirichlet-proxy/real").glob("seed-*/metrics.json")), "dirichlet-proxy")
    proxy_shuffled = summarize(sorted((args.root / "constituents/dirichlet-proxy/shuffled").glob("seed-*/metrics.json")), "dirichlet-proxy")

    for arm in (hybrid_real, fixed_real, proxy_real):
        if [row["n_by_seed"] for row in arm["evidence_depth_strata"]] != [
            row["n_by_seed"] for row in hybrid_real["evidence_depth_strata"]
        ]:
            raise SystemExit("real-arm depth populations are not paired")

    oracle_by_seed = []
    for seed_index, seed in enumerate((7, 17, 29)):
        total = 0.0
        n = int(hybrid_real["by_seed"][seed_index]["n"])
        for depth_index in range(4):
            source = fixed_real if depth_index < 3 else proxy_real
            item = source["evidence_depth_strata"][depth_index]
            by_seed = item["by_seed"][seed_index]
            total += int(by_seed["n"]) * float(by_seed["negative_log_loss"])
        oracle_by_seed.append({"seed": seed, "n": n, "negative_log_loss": total / n})
    oracle_nll = statistics.fmean(row["negative_log_loss"] for row in oracle_by_seed)

    hybrid_nll = float(hybrid_real["means"]["negative_log_loss"])
    proxy_nll = float(proxy_real["means"]["negative_log_loss"])
    hybrid_shuffled_nll = float(hybrid_shuffled["means"]["negative_log_loss"])
    proxy_shuffled_nll = float(proxy_shuffled["means"]["negative_log_loss"])
    gain = proxy_nll - hybrid_nll
    oracle_gain = proxy_nll - oracle_nll
    shuffle_fraction = (proxy_shuffled_nll - hybrid_shuffled_nll) / gain if gain > 0 else None

    depth_comparison = []
    for index, hybrid in enumerate(hybrid_real["evidence_depth_strata"]):
        reference = fixed_real["evidence_depth_strata"][index] if index < 3 else proxy_real["evidence_depth_strata"][index]
        reference_nll = float(reference["mean_negative_log_loss"])
        hybrid_depth_nll = float(hybrid["mean_negative_log_loss"])
        depth_comparison.append(
            {
                "name": hybrid["name"],
                "n_by_seed": hybrid["n_by_seed"],
                "reference": "fixed-convex" if index < 3 else "dirichlet-proxy",
                "hybrid_mean_negative_log_loss": hybrid_depth_nll,
                "reference_mean_negative_log_loss": reference_nll,
                "relative_nll_change_vs_reference": (hybrid_depth_nll - reference_nll) / reference_nll,
            }
        )

    real_beats_shuffled = [
        float(real["negative_log_loss"]) < float(shuffled["negative_log_loss"])
        for real, shuffled in zip(hybrid_real["by_seed"], hybrid_shuffled["by_seed"], strict=True)
    ]
    criteria = {
        "mean_nll_improvement_over_proxy_at_least_0_2_percent": gain / proxy_nll >= 0.002,
        "hard_stratum_oracle_gain_recovered_at_least_half": oracle_gain > 0 and gain / oracle_gain >= 0.50,
        "each_shallow_stratum_regression_vs_fixed_at_most_0_25_percent": all(
            float(row["relative_nll_change_vs_reference"]) <= 0.0025 for row in depth_comparison[:3]
        ),
        "high_stratum_regression_vs_proxy_at_most_0_25_percent": float(depth_comparison[3]["relative_nll_change_vs_reference"]) <= 0.0025,
        "real_hybrid_beats_shuffled_on_all_seeds": all(real_beats_shuffled),
        "shuffle_explained_fraction_at_most_0_50": shuffle_fraction is not None and shuffle_fraction <= 0.50,
        "mean_brier_does_not_increase_vs_proxy": float(hybrid_real["means"]["multiclass_brier"]) <= float(proxy_real["means"]["multiclass_brier"]),
        "top1_drop_vs_best_constituent_at_most_0_10_percentage_point": max(float(fixed_real["means"]["top1_percent"]), float(proxy_real["means"]["top1_percent"])) - float(hybrid_real["means"]["top1_percent"]) <= 0.10,
    }
    if all(criteria.values()):
        verdict = "PASS"
        recommendation = "lock_untouched_D4_confirmation"
    elif gain > 0 and all(real_beats_shuffled):
        verdict = "MARGINAL"
        recommendation = "review_failed_acceptance_criterion_before_any_confirmation"
    else:
        verdict = "STOP"
        recommendation = "retain_constituent_models_as_diagnostics_only"

    result = {
        "status": "D0_development_gate_complete",
        "analysis": "monotone_evidence_depth_hybrid",
        "protocol": parse_protocol(args.root / "protocol.tsv"),
        "selection": selection,
        "D0": {
            "hybrid_real": hybrid_real,
            "hybrid_shuffled": hybrid_shuffled,
            "fixed_convex_real": fixed_real,
            "fixed_convex_shuffled": fixed_shuffled,
            "Dirichlet_proxy_real": proxy_real,
            "Dirichlet_proxy_shuffled": proxy_shuffled,
            "hard_stratum_oracle": {"mean_negative_log_loss": oracle_nll, "by_seed": oracle_by_seed},
            "mean_relative_nll_improvement_over_proxy": gain / proxy_nll,
            "hard_stratum_oracle_gain_recovered": gain / oracle_gain if oracle_gain > 0 else None,
            "shuffle_explained_fraction_of_hybrid_gain": shuffle_fraction,
            "real_hybrid_better_than_shuffled_by_seed": real_beats_shuffled,
            "evidence_depth_comparison": depth_comparison,
            "posthoc_factorial_contrasts": {
                "status": "descriptive_not_part_of_locked_verdict",
                "hybrid_gain_with_real_groups_relative_to_proxy": gain / proxy_nll,
                "hybrid_gain_with_shuffled_groups_relative_to_proxy":
                    (proxy_shuffled_nll - hybrid_shuffled_nll) / proxy_shuffled_nll,
                "real_group_gain_with_hybrid_relative_to_shuffled":
                    (hybrid_shuffled_nll - hybrid_nll) / hybrid_shuffled_nll,
                "real_group_gain_with_proxy_relative_to_shuffled":
                    (proxy_shuffled_nll - proxy_nll) / proxy_shuffled_nll,
                "interpretation": "depth-model and real-group main effects are distinct; the locked ratio compares those axes and remains binding",
            },
        },
        "acceptance": criteria,
        "verdict": verdict,
        "recommendation": recommendation,
        "confirmation_status": "none_D0_is_development_only",
        "full_hierarchical_Dirichlet_licensed": False,
        "production_emission": "pooled",
        "performance": timing(args.root),
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
