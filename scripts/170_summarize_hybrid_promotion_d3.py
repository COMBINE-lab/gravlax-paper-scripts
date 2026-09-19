#!/usr/bin/env python3
"""Summarize the single untouched D3 hybrid-promotion confirmation."""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
from pathlib import Path


def load_audit_helpers() -> object:
    path = Path(__file__).with_name("167_summarize_hybrid_promotion_selected.py")
    spec = importlib.util.spec_from_file_location("hybrid_audit", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


audit = load_audit_helpers()


def hybrid_mode(payload: dict[str, object]) -> dict[str, object]:
    return next(mode for mode in payload["modes"] if mode["name"] == "depth-hybrid")


def weighted(rows: list[dict[str, object]], metric: str) -> float:
    total = sum(int(row["n"]) for row in rows)
    return sum(int(row["n"]) * float(row[metric]) for row in rows) / total


def performance(root: Path) -> dict[str, float]:
    walls, rss = [], []
    for path in root.glob("*/seed-*/time.txt"):
        text = path.read_text()
        wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)", text)
        memory = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
        if wall is None or memory is None:
            raise SystemExit(f"incomplete timing: {path}")
        walls.append(audit.elapsed_seconds(wall.group(1)))
        rss.append(int(memory.group(1)))
    if len(walls) != 6:
        raise SystemExit(f"expected six timed runs, found {len(walls)}")
    return {"maximum_wall_seconds": max(walls), "maximum_peak_rss_kib": max(rss)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    real = [json.loads(path.read_text()) for path in (args.root / "real").glob("seed-*/metrics.json")]
    shuffled = [json.loads(path.read_text()) for path in (args.root / "shuffled").glob("seed-*/metrics.json")]
    real.sort(key=lambda row: (7, 17, 29).index(int(row["seed"])))
    shuffled.sort(key=lambda row: (7, 17, 29).index(int(row["seed"])))
    if [int(row["seed"]) for row in real] != [7, 17, 29] or [int(row["seed"]) for row in shuffled] != [7, 17, 29]:
        raise SystemExit("expected real and shuffled seeds 7, 17, 29")

    modes = audit.aggregate_modes(real)
    comparisons = {ref: audit.aggregate_comparison(real, ref) for ref in ("pooled", "dirichlet-proxy")}
    real_hybrid = [hybrid_mode(row) for row in real]
    shuffled_hybrid = [hybrid_mode(row) for row in shuffled]
    group_wins = [
        float(r["negative_log_loss"]) < float(s["negative_log_loss"])
        for r, s in zip(real_hybrid, shuffled_hybrid, strict=True)
    ]
    real_nll = weighted(real_hybrid, "negative_log_loss")
    shuffled_nll = weighted(shuffled_hybrid, "negative_log_loss")
    group_contrast = {
        "real_mean_negative_log_loss": real_nll,
        "shuffled_mean_negative_log_loss": shuffled_nll,
        "relative_nll_improvement_real_vs_shuffled": (shuffled_nll - real_nll) / shuffled_nll,
        "real_better_by_seed": group_wins,
        "wins": sum(group_wins),
    }
    pooled = comparisons["pooled"]
    proxy = comparisons["dirichlet-proxy"]
    criteria = {
        "NLL_upper_CI_below_zero_vs_pooled": float(pooled["negative_log_loss_ci95"][1]) < 0,
        "Brier_upper_CI_at_most_zero_vs_pooled": float(pooled["brier_ci95"][1]) <= 0,
        "mean_NLL_below_proxy": float(modes["depth-hybrid"]["negative_log_loss"]) < float(modes["dirichlet-proxy"]["negative_log_loss"]),
        "mean_Brier_at_most_proxy": float(modes["depth-hybrid"]["multiclass_brier"]) <= float(modes["dirichlet-proxy"]["multiclass_brier"]),
        "top1_drop_at_most_0_10pp_vs_pooled": float(pooled["top1_percentage_point_difference"]) >= -0.10,
        "top1_drop_at_most_0_10pp_vs_proxy": float(proxy["top1_percentage_point_difference"]) >= -0.10,
        "real_groups_beat_shuffled_on_at_least_2_seeds": sum(group_wins) >= 2,
    }
    perf = performance(args.root)
    criteria["peak_RSS_below_4_GiB"] = perf["maximum_peak_rss_kib"] < 4 * 1024 * 1024
    passed = all(criteria.values())
    result = {
        "status": "untouched_D3_confirmation_complete",
        "candidate": {"depth_scale": 8, "depth_power": 8},
        "modes": modes,
        "paired_comparisons": comparisons,
        "group_information": group_contrast,
        "acceptance": criteria,
        "verdict": "PASS" if passed else "FAIL",
        "decision": "promote_hybrid_default_evaluator_and_implement_opt_in_emission" if passed else "retain_pooled_default_and_label_hybrid_experimental",
        "performance": perf,
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()

