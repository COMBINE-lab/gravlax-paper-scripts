#!/usr/bin/env python3
"""Plot the replicated FYB1 positive control and prospective FNBP1 confirmation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--confirmation", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    screen = json.loads(args.screen.read_text())
    confirmation = json.loads(args.confirmation.read_text())
    fyb1 = screen["selected_events"]["FYB1"]["samples"]
    fnbp1 = screen["selected_events"]["FNBP1"]["samples"]

    plt.rcParams.update(
        {
            "font.size": 8.5,
            "font.family": "sans-serif",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
            "axes.axisbelow": True,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(6.55, 2.35), sharey=True)

    pbmc = ["D0", "D1", "D3", "D4"]
    x = np.arange(len(pbmc))
    for index, dataset in enumerate(pbmc):
        axes[0].plot(
            [index, index],
            [100 * fyb1[dataset]["T_usage"], 100 * fyb1[dataset]["M_usage"]],
            color="#b9c0c5",
            linewidth=1.0,
            zorder=1,
        )
    axes[0].scatter(
        x,
        [100 * fyb1[dataset]["T_usage"] for dataset in pbmc],
        color="#1a6faf",
        label="T cells",
        zorder=3,
    )
    axes[0].scatter(
        x,
        [100 * fyb1[dataset]["M_usage"] for dataset in pbmc],
        color="#c98f1a",
        label="monocytes",
        zorder=3,
    )
    axes[0].set_xticks(x, pbmc)
    axes[0].set_title("a  FYB1 external positive control", loc="left", fontsize=9)
    axes[0].legend(frameon=False, loc="center right", fontsize=7.5)

    order = ["D0", "D1", "D3", "D4", "D2", "D2p"]
    labels = ["D0", "D1", "D3", "D4", "D2\nGBM", "D2′\nbrain"]
    colors = ["#687982"] * 4 + ["#a65c84", "#6d55a3"]
    values = [100 * fnbp1[dataset]["usage"] for dataset in order]
    axes[1].scatter(np.arange(6), values, color=colors, s=27, zorder=3)
    held_out = 100 * confirmation["annotation_free_targeted_STAR_archive"]["totals"][
        "usage_fraction"
    ]
    axes[1].scatter(
        [6], [held_out], marker="*", s=80, color="#168567", edgecolor="white",
        linewidth=0.5, zorder=4, label="locked confirmation",
    )
    axes[1].set_xticks(np.arange(7), labels + ["D5\nbrain"])
    axes[1].set_title("b  FNBP1 prospective replication", loc="left", fontsize=9)
    axes[1].legend(frameon=False, loc="center right", fontsize=7.2)

    for axis in axes:
        axis.set_ylim(-4, 104)
        axis.set_yticks([0, 25, 50, 75, 100])
        axis.yaxis.grid(True, color="#e1e5e7", linewidth=0.7)
        axis.tick_params(axis="both", length=0, labelsize=7.5)
    axes[0].set_ylabel("cassette inclusion (%)")
    fig.tight_layout(w_pad=1.5)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")


if __name__ == "__main__":
    main()
