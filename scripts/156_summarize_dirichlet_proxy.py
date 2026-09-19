#!/usr/bin/env python3
"""Evaluate the locked D0 Dirichlet-proxy screen against matched controls."""
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
        raise SystemExit(f"expected {expected} metrics for seeds 7, 17, and 29")
    depth = []
    for name in DEPTH_NAMES:
        strata = []
        for row in rows:
            by_name = {stratum["name"]: stratum for stratum in row["evidence_depth_strata"]}
            if tuple(by_name) != DEPTH_NAMES:
                raise SystemExit(f"unexpected depth strata in {expected}, seed {row['seed']}")
            strata.append(by_name[name])
        n_by_seed = [int(stratum["n"]) for stratum in strata]
        depth.append(
            {
                "name": name,
                "n_by_seed": n_by_seed,
                "minimum_n": min(n_by_seed),
                "mean_n": statistics.fmean(n_by_seed),
                **{
                    f"mean_{metric}": statistics.fmean(float(stratum[metric]) for stratum in strata)
                    for metric in METRICS
                },
            }
        )
    return {
        "means": {
            metric: statistics.fmean(float(row[metric]) for row in rows) for metric in METRICS
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


def timing(root: Path) -> dict[str, object]:
    walls: list[float] = []
    rss: list[int] = []
    for path in root.glob("**/time.txt"):
        text = path.read_text()
        wall_match = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)", text)
        rss_match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
        if wall_match is None or rss_match is None:
            raise SystemExit(f"incomplete timing file: {path}")
        walls.append(elapsed_seconds(wall_match.group(1)))
        rss.append(int(rss_match.group(1)))
    if len(walls) != 84:
        raise SystemExit(f"expected 84 timings, found {len(walls)}")
    return {
        "runs": len(walls),
        "wall_seconds": {
            "minimum": min(walls),
            "median": statistics.median(walls),
            "maximum": max(walls),
        },
        "peak_rss_kib": {"minimum": min(rss), "median": statistics.median(rss), "maximum": max(rss)},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    selection = json.loads((args.root / "selection.json").read_text())
    chosen = selection["selected"]
    label = f"kh-{float(chosen['cell_prior']):g}.k0-{float(chosen['group_prior']):g}"
    proxy_real = summarize(sorted((args.root / "grid" / label).glob("seed-*/metrics.json")), "dirichlet-proxy")
    proxy_shuffled = summarize(sorted((args.root / "selected-shuffled").glob("seed-*/metrics.json")), "dirichlet-proxy")
    fixed_real = summarize(sorted((args.root / "fixed-convex/real").glob("seed-*/metrics.json")), "convex")
    fixed_shuffled = summarize(sorted((args.root / "fixed-convex/shuffled").glob("seed-*/metrics.json")), "convex")

    proxy_real_nll = float(proxy_real["means"]["negative_log_loss"])
    proxy_shuffled_nll = float(proxy_shuffled["means"]["negative_log_loss"])
    fixed_real_nll = float(fixed_real["means"]["negative_log_loss"])
    fixed_shuffled_nll = float(fixed_shuffled["means"]["negative_log_loss"])
    real_gain = fixed_real_nll - proxy_real_nll
    shuffled_gain = fixed_shuffled_nll - proxy_shuffled_nll
    shuffle_fraction = shuffled_gain / real_gain if real_gain > 0 else None

    real_beats_shuffled_by_seed = [
        float(real["negative_log_loss"]) < float(shuffled["negative_log_loss"])
        for real, shuffled in zip(proxy_real["by_seed"], proxy_shuffled["by_seed"], strict=True)
    ]
    depth_comparison = []
    for proxy, fixed in zip(
        proxy_real["evidence_depth_strata"], fixed_real["evidence_depth_strata"], strict=True
    ):
        if proxy["name"] != fixed["name"] or proxy["n_by_seed"] != fixed["n_by_seed"]:
            raise SystemExit("proxy and fixed-convex depth strata are not paired")
        fixed_nll = float(fixed["mean_negative_log_loss"])
        proxy_nll = float(proxy["mean_negative_log_loss"])
        depth_comparison.append(
            {
                "name": proxy["name"],
                "n_by_seed": proxy["n_by_seed"],
                "minimum_n": proxy["minimum_n"],
                "proxy_mean_negative_log_loss": proxy_nll,
                "fixed_convex_mean_negative_log_loss": fixed_nll,
                "relative_nll_improvement": (fixed_nll - proxy_nll) / fixed_nll,
            }
        )

    populated = [row for row in depth_comparison if int(row["minimum_n"]) >= 1000]
    criteria = {
        "aggregate_nll_improvement_at_least_0_25_percent": real_gain / fixed_real_nll >= 0.0025,
        "real_proxy_beats_shuffled_on_all_seeds": all(real_beats_shuffled_by_seed),
        "shuffle_explained_fraction_at_most_0_50": shuffle_fraction is not None and shuffle_fraction <= 0.50,
        "mean_brier_does_not_increase": float(proxy_real["means"]["multiclass_brier"]) <= float(fixed_real["means"]["multiclass_brier"]),
        "top1_drop_at_most_0_10_percentage_point": float(fixed_real["means"]["top1_percent"]) - float(proxy_real["means"]["top1_percent"]) <= 0.10,
        "at_least_two_populated_depth_strata_improve": sum(float(row["relative_nll_improvement"]) > 0 for row in populated) >= 2,
        "low_or_medium_depth_stratum_improves_at_least_0_25_percent": any(
            row["name"] in DEPTH_NAMES[:3] and float(row["relative_nll_improvement"]) >= 0.0025
            for row in populated
        ),
        "no_populated_depth_stratum_regresses_over_0_5_percent": all(
            float(row["relative_nll_improvement"]) >= -0.005 for row in populated
        ),
    }
    if all(criteria.values()):
        verdict = "PASS"
        recommendation = "lock_a_separate_full_hierarchical_Dirichlet_implementation_gate"
    elif real_gain > 0 and all(real_beats_shuffled_by_seed) and (shuffle_fraction is None or shuffle_fraction <= 0.50):
        verdict = "MARGINAL"
        recommendation = "do_not_implement_full_model_without_new_signal"
    else:
        verdict = "STOP"
        recommendation = "stop_hierarchy_ladder_and_keep_pooled_production"

    result = {
        "status": "D0_development_screen_complete",
        "analysis": "posterior_mean_Dirichlet_proxy_not_full_variational_model",
        "protocol": parse_protocol(args.root / "protocol.tsv"),
        "selection": selection,
        "D0": {
            "proxy_real": proxy_real,
            "proxy_shuffled": proxy_shuffled,
            "fixed_convex_real": fixed_real,
            "fixed_convex_shuffled": fixed_shuffled,
            "mean_relative_nll_improvement_over_fixed_convex": real_gain / fixed_real_nll,
            "real_proxy_better_than_shuffled_by_seed": real_beats_shuffled_by_seed,
            "shuffle_explained_fraction_of_proxy_gain": shuffle_fraction,
            "top1_drop_percentage_points": float(fixed_real["means"]["top1_percent"]) - float(proxy_real["means"]["top1_percent"]),
            "evidence_depth_comparison": depth_comparison,
        },
        "acceptance": criteria,
        "verdict": verdict,
        "recommendation": recommendation,
        "confirmation_status": "none_D0_is_development_only",
        "production_emission": "pooled",
        "performance": timing(args.root),
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
