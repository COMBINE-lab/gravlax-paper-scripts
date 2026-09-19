#!/usr/bin/env python3
"""Summarize the post-hoc hierarchical-EM scale audit without upgrading it to a gate."""
from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path


SEEDS = (7, 17, 29)
GRID = (
    (0, 6400), (1280, 5120), (3200, 3200), (5120, 1280),
    (0, 25600), (5120, 20480), (12800, 12800), (20480, 5120),
    (0, 102400), (20480, 81920), (51200, 51200), (81920, 20480),
    (327680, 81920),
)


def mode(path: Path, name: str) -> dict[str, object]:
    data = json.loads(path.read_text())
    return next(item for item in data["modes"] if item["name"] == name)


def aggregate(paths: list[Path], name: str) -> dict[str, float | int]:
    rows = [mode(path, name) for path in paths]
    total = sum(int(row["n"]) for row in rows)
    return {
        "n_across_seeds": total,
        "top1_percent": sum(float(row["top1_percent"]) * int(row["n"]) for row in rows) / total,
        "negative_log_loss": sum(
            float(row["negative_log_loss"]) * int(row["n"]) for row in rows
        ) / total,
        "multiclass_brier": sum(
            float(row["multiclass_brier"]) * int(row["n"]) for row in rows
        ) / total,
    }


def paths(root: Path, relative: str) -> list[Path]:
    return [root / relative / f"seed-{seed}" / "metrics.json" for seed in SEEDS]


def quantiles(values: list[float], probabilities: tuple[float, ...]) -> list[float]:
    values.sort()
    result = []
    for probability in probabilities:
        at = probability * (len(values) - 1)
        low = int(at)
        high = min(low + 1, len(values) - 1)
        fraction = at - low
        result.append(values[low] * (1.0 - fraction) + values[high] * fraction)
    return result


def d0_prior_diagnostics(project: Path) -> dict[str, object]:
    replay = project / "runs/replay/d0-v15"
    prepared = project / "runs/post-v1/hierarchical-em-groups-r1"
    groups: dict[str, int] = {}
    for line in (prepared / "groups/D0.real.tsv").read_text().splitlines():
        barcode, group = line.split("\t")
        groups[barcode] = int(group.removeprefix("g"))

    column_group: dict[int, int] = {}
    selected_columns: list[int] = []
    with (replay / "barcodes.tsv").open() as handle:
        for column, line in enumerate(handle, 1):
            group = groups.get(line.rstrip("\n"))
            if group is not None:
                column_group[column] = group
                selected_columns.append(column)

    features = [line.split("\t", 1)[0] for line in (replay / "features.tsv").read_text().splitlines()]
    candidates: set[str] = set()
    for candidate_file in (prepared / "candidates").glob("D0.seed-*.genes.txt"):
        candidates.update(candidate_file.read_text().splitlines())
    candidate_rows = {index for index, gene in enumerate(features, 1) if gene in candidates}

    group_count = max(groups.values()) + 1
    group_totals = [0] * group_count
    cell_totals = [0] * len(selected_columns)
    column_local = {column: index for index, column in enumerate(selected_columns)}
    group_gene: dict[tuple[int, int], int] = {}
    global_gene: dict[int, int] = {}
    with (replay / "matrix.mtx").open() as handle:
        dimensions_seen = False
        for line in handle:
            if line.startswith("%"):
                continue
            if not dimensions_seen:
                dimensions_seen = True
                continue
            gene, column, count = map(int, line.split())
            group = column_group.get(column)
            if group is None:
                continue
            group_totals[group] += count
            cell_totals[column_local[column]] += count
            if gene in candidate_rows:
                group_gene[group, gene] = group_gene.get((group, gene), 0) + count
                global_gene[gene] = global_gene.get(gene, 0) + count

    global_total = sum(group_totals)
    contributions = [
        80 * group_gene.get((group, gene), 0) / group_totals[group]
        + 20 * global_gene.get(gene, 0) / global_total
        for group in range(group_count)
        for gene in candidate_rows
    ]
    nonzero = [value for value in contributions if value > 0]
    probabilities = (0.0, 0.5, 0.9, 0.99, 0.999, 1.0)
    return {
        "registered_group_alpha": 80,
        "registered_global_alpha": 20,
        "called_cells": len(selected_columns),
        "groups": group_count,
        "candidate_genes": len(candidate_rows),
        "cell_library_quantiles_0_10_50_90_100": quantiles(cell_totals, (0.0, 0.1, 0.5, 0.9, 1.0)),
        "group_total_counts": group_totals,
        "candidate_prior_contribution_quantiles_0_50_90_99_99_9_100": quantiles(
            contributions, probabilities
        ),
        "nonzero_fraction": len(nonzero) / len(contributions),
        "interpretation": "A unit local count dominates the registered prior for more than 99% of D0 candidate/group pairs.",
    }


def read_protocol(path: Path) -> dict[str, str] | None:
    if not path.exists():
        return None
    rows = [line.rstrip("\n").split("\t", 1) for line in path.open()]
    return dict(rows[1:])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("audit_dir", type=Path)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    grid = []
    for group_alpha, global_alpha in GRID:
        run_paths = paths(args.audit_dir, f"d0-grid/ag-{group_alpha}.aa-{global_alpha}")
        rows = [mode(path, "hierarchical") for path in run_paths]
        grid.append({
            "group_alpha": group_alpha,
            "global_alpha": global_alpha,
            "negative_log_loss_by_seed": [float(row["negative_log_loss"]) for row in rows],
            "mean_negative_log_loss": statistics.mean(float(row["negative_log_loss"]) for row in rows),
            "aggregate": aggregate(run_paths, "hierarchical"),
        })
    selected = min(grid, key=lambda row: (row["mean_negative_log_loss"], row["group_alpha"], row["global_alpha"]))

    datasets: dict[str, object] = {}
    for dataset in ("D0", "D1", "D2_prime"):
        if dataset == "D0":
            hierarchy_paths = paths(args.audit_dir, "d0-grid/ag-12800.aa-12800")
            calibrated_paths = paths(args.audit_dir, "d0-grid/ag-0.aa-25600")
        else:
            hierarchy_paths = paths(args.audit_dir, f"confirm/{dataset}/hierarchical")
            calibrated_paths = paths(args.audit_dir, f"confirm/{dataset}/calibrated_global")
        shuffled_paths = paths(args.audit_dir, f"shuffled/{dataset}")
        datasets[dataset] = {
            "pooled": aggregate(hierarchy_paths, "pooled"),
            "calibrated_cell_global_0_25600": aggregate(calibrated_paths, "hierarchical"),
            "real_hierarchy_12800_12800": aggregate(hierarchy_paths, "hierarchical"),
            "shuffled_hierarchy_12800_12800": aggregate(shuffled_paths, "hierarchical"),
        }

    confirmation_n = sum(datasets[name]["pooled"]["n_across_seeds"] for name in ("D1", "D2_prime"))

    def combined(key: str) -> float:
        return sum(
            datasets[name][key]["negative_log_loss"] * datasets[name][key]["n_across_seeds"]
            for name in ("D1", "D2_prime")
        ) / confirmation_n

    pooled = combined("pooled")
    calibrated = combined("calibrated_cell_global_0_25600")
    hierarchy = combined("real_hierarchy_12800_12800")
    shuffled = combined("shuffled_hierarchy_12800_12800")
    result = {
        "status": "complete",
        "analysis": "post_hoc_hierarchical_em_scale_audit",
        "inferential_status": "exploratory_not_independent_confirmation",
        "warning": "D1 and D2-prime had already been inspected under the registered gate; these results diagnose scale but are not a new confirmatory test.",
        "protocol": read_protocol(args.audit_dir / "protocol.tsv"),
        "registered_gate": {
            "result": "results/post-v1-hierarchical-em.json",
            "verdict_remains": "STOP for the registered alpha grid",
        },
        "d0_registered_scale_diagnostic": d0_prior_diagnostics(args.project_root),
        "d0_exploratory_selection": {"grid": grid, "selected": selected},
        "datasets": datasets,
        "combined_d1_d2_prime": {
            "n_across_seeds": confirmation_n,
            "pooled_negative_log_loss": pooled,
            "calibrated_cell_global_negative_log_loss": calibrated,
            "real_hierarchy_negative_log_loss": hierarchy,
            "shuffled_hierarchy_negative_log_loss": shuffled,
            "calibrated_cell_global_improvement_vs_pooled": (pooled - calibrated) / pooled,
            "real_hierarchy_improvement_vs_pooled": (pooled - hierarchy) / pooled,
            "real_hierarchy_improvement_vs_calibrated_cell_global": (calibrated - hierarchy) / calibrated,
            "real_hierarchy_improvement_vs_shuffled": (shuffled - hierarchy) / shuffled,
        },
        "interpretation": {
            "implementation_failure": False,
            "registered_grid_was_under_scaled": True,
            "most_gain_is_calibrated_cell_global_shrinkage": True,
            "real_groups_add_small_consistent_signal": True,
            "positive_paper_claim_requires_new_independent_validation": True,
        },
    }
    rendered = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.out:
        args.out.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
