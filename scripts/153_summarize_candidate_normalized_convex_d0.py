#!/usr/bin/env python3
"""Summarize the locked D0 convex grid and selected shuffled-label control."""
from __future__ import annotations

import argparse
import json
import re
import statistics
from pathlib import Path


METRICS = ("negative_log_loss", "multiclass_brier", "top1_percent")


def load_mode(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    modes = payload["modes"]
    if len(modes) != 1 or modes[0]["name"] != "convex":
        raise SystemExit(f"expected one convex mode: {path}")
    return {"seed": int(payload["seed"]), **modes[0]}


def summarize(paths: list[Path]) -> dict[str, object]:
    rows = sorted((load_mode(path) for path in paths), key=lambda row: int(row["seed"]))
    if [row["seed"] for row in rows] != [7, 17, 29]:
        raise SystemExit("expected selected metrics for seeds 7, 17, and 29")
    return {
        "means": {
            metric: statistics.fmean(float(row[metric]) for row in rows) for metric in METRICS
        },
        "by_seed": [
            {
                "seed": row["seed"],
                "n": row["n"],
                **{metric: row[metric] for metric in METRICS},
            }
            for row in rows
        ],
    }


def parse_protocol(path: Path) -> dict[str, str]:
    rows = [line.rstrip("\n").split("\t", 1) for line in path.open()]
    return dict(rows[1:])


def timing(root: Path) -> dict[str, object]:
    walls: list[float] = []
    rss: list[int] = []
    for path in root.glob("grid/*/seed-*/time.txt"):
        text = path.read_text()
        wall_match = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)", text)
        rss_match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
        if wall_match is None or rss_match is None:
            raise SystemExit(f"incomplete timing file: {path}")
        fields = [float(value) for value in wall_match.group(1).split(":")]
        walls.append(sum(value * 60**power for power, value in enumerate(reversed(fields))))
        rss.append(int(rss_match.group(1)))
    if len(walls) != 204:
        raise SystemExit(f"expected 204 grid timings, found {len(walls)}")
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
    parser.add_argument("--hierarchical-result", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    selection = json.loads((args.root / "selection.json").read_text())
    chosen = selection["selected"]
    label = (
        f"cw-{float(chosen['cell_weight']):.2f}."
        f"gw-{float(chosen['group_weight']):.4f}."
        f"k-{float(chosen['group_prior']):g}"
    )
    real = summarize(sorted((args.root / "grid" / label).glob("seed-*/metrics.json")))
    shuffled = summarize(sorted((args.root / "selected-shuffled").glob("seed-*/metrics.json")))
    no_group = selection["tuned_no_group"]
    no_group_nll = float(no_group["mean_negative_log_loss"])
    real_nll = float(real["means"]["negative_log_loss"])
    shuffled_nll = float(shuffled["means"]["negative_log_loss"])
    group_gain = no_group_nll - real_nll
    shuffled_gain = no_group_nll - shuffled_nll

    old = json.loads(args.hierarchical_result.read_text())
    pooled = old["datasets"]["D0"]["controls"]["real"]["pooled"]
    real_beats_shuffled = all(
        float(real_row["negative_log_loss"]) < float(shuffled_row["negative_log_loss"])
        for real_row, shuffled_row in zip(real["by_seed"], shuffled["by_seed"], strict=True)
    )
    result = {
        "status": "D0_development_complete_D4_confirmation_pending",
        "analysis": "candidate_normalized_convex_partial_pooling",
        "protocol": parse_protocol(args.root / "protocol.tsv"),
        "control_execution_scripts_commit": "d83ed08",
        "selection": {
            "selected": chosen,
            "tuned_no_group": no_group,
            "configuration_count": selection["configuration_count"],
            "selection_metric": selection["selection_rule"],
        },
        "D0": {
            "real_selected": real,
            "shuffled_selected": shuffled,
            "pooled_reference": pooled,
            "relative_nll_improvement_over_pooled":
                (float(pooled["negative_log_loss"]) - real_nll)
                / float(pooled["negative_log_loss"]),
            "relative_nll_improvement_over_tuned_no_group": group_gain / no_group_nll,
            "relative_nll_real_over_shuffled": (shuffled_nll - real_nll) / shuffled_nll,
            "shuffle_explained_fraction_of_group_gain": shuffled_gain / group_gain,
            "real_beats_shuffled_on_all_seeds": real_beats_shuffled,
        },
        "performance": timing(args.root),
        "interpretation": {
            "confirmation_status": "none_D0_is_development_only",
            "candidate_normalized_calibration_signal": "encouraging",
            "biological_group_signal": "small_and_most_gain_is_reproduced_by_shuffling",
            "full_dirichlet_escalation_supported_now": False,
            "next_binding_gate": "untouched_D4_real_vs_shuffled_and_tuned_no_group",
            "production_emission": "pooled",
        },
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
