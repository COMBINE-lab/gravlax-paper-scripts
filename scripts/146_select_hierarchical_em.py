#!/usr/bin/env python3
"""Select hierarchical shrinkage once on the registered D0 development grid."""
import argparse
import json
import statistics
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("grid_dir", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for group_alpha in (5, 20, 80):
        for global_alpha in (0, 5, 20):
            values = []
            for seed in (7, 17, 29):
                path = args.grid_dir / f"ag-{group_alpha}.aa-{global_alpha}" / f"seed-{seed}" / "metrics.json"
                data = json.loads(path.read_text())
                mode = next(item for item in data["modes"] if item["name"] == "hierarchical")
                values.append(float(mode["negative_log_loss"]))
            rows.append(
                {
                    "group_alpha": group_alpha,
                    "global_alpha": global_alpha,
                    "negative_log_loss_by_seed": values,
                    "mean_negative_log_loss": statistics.mean(values),
                }
            )
    selected = min(
        rows,
        key=lambda item: (
            item["mean_negative_log_loss"],
            item["group_alpha"],
            item["global_alpha"],
        ),
    )
    args.out.write_text(
        json.dumps(
            {
                "status": "complete",
                "selection_dataset": "D0",
                "selection_metric": "mean_negative_log_loss",
                "grid": rows,
                "selected": selected,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
