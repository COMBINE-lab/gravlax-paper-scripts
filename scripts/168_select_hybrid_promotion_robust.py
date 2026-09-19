#!/usr/bin/env python3
"""Select the single locked robustness candidate from the eligible development grid."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("development_selection", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.development_selection.read_text())
    references = source["reference"]
    eligible = [row for row in source["configurations"] if row["eligible"]]
    if not eligible:
        raise SystemExit("original development grid contains no eligible candidates")
    for row in eligible:
        gains = {
            dataset: (
                float(references[dataset]["dirichlet-proxy"]["multiclass_brier"])
                - float(row["by_dataset"][dataset]["multiclass_brier"])
            )
            / float(references[dataset]["dirichlet-proxy"]["multiclass_brier"])
            for dataset in ("D0", "D4")
        }
        row["relative_mean_Brier_improvement_vs_proxy"] = gains
        row["minimum_relative_mean_Brier_improvement_vs_proxy"] = min(gains.values())
    selected = min(
        eligible,
        key=lambda row: (
            -row["minimum_relative_mean_Brier_improvement_vs_proxy"],
            row["joint_target_weighted_negative_log_loss"],
            row["depth_power"],
            row["depth_scale"],
        ),
    )
    args.out.write_text(
        json.dumps(
            {
                "status": "ELIGIBLE_MODEL_SELECTED",
                "stage": "transparent_post_hoc_D0_D4_robustness_selection",
                "selection_rule": (
                    "maximize minimum relative mean Brier improvement versus proxy across D0/D4; "
                    "tie by joint NLL, lower power, lower scale"
                ),
                "reference": references,
                "eligible_configuration_count": len(eligible),
                "selected": selected,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()

