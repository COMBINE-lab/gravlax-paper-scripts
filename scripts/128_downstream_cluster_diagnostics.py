#!/usr/bin/env python3
"""Post-gate diagnostics for Gene clustering failures; does not alter frozen verdicts."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import scipy.optimize
import scipy.sparse as sp
import yaml


def load_gate_module():
    path = Path(__file__).with_name("126_downstream_fidelity.py")
    spec = importlib.util.spec_from_file_location("downstream_gate", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


gate = load_gate_module()


def feature_symbols(path: Path) -> list[str]:
    symbols = []
    for line in path.open():
        fields = line.rstrip("\n").split("\t")
        symbols.append(fields[1] if len(fields) > 1 else fields[0])
    return symbols


def nearest_centroid_labels(
    oracle_pcs: np.ndarray, replay_pcs: np.ndarray, labels: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    categories = sorted(set(labels), key=lambda value: int(value))
    centroids = np.vstack([oracle_pcs[labels == value].mean(axis=0) for value in categories])

    def predict(pcs: np.ndarray) -> np.ndarray:
        squared = ((pcs[:, None, :] - centroids[None, :, :]) ** 2).sum(axis=2)
        return np.asarray([categories[i] for i in squared.argmin(axis=1)])

    return predict(oracle_pcs), predict(replay_pcs)


def analyze(config: dict[str, str], root: Path) -> dict[str, object]:
    dataset = config["id"]
    oracle_dir = root / config["oracle_gene"]
    replay_dir = root / config["replay_gene"]
    barcodes = gate.read_lines(oracle_dir / "barcodes.tsv")
    ids = gate.require_same_features(
        oracle_dir / "features.tsv", replay_dir / "features.tsv"
    )
    symbols = feature_symbols(oracle_dir / "features.tsv")
    columns = gate.replay_column_indices(replay_dir / "barcodes.tsv", barcodes)
    oracle_raw = gate.load_mtx_cells(oracle_dir / "matrix.mtx")
    replay_raw = gate.load_mtx_cells(replay_dir / "matrix.mtx", columns)
    detected = np.asarray(oracle_raw.getnnz(axis=0)).ravel() >= 10
    oracle_raw = oracle_raw[:, detected]
    replay_raw = replay_raw[:, detected]
    ids = np.asarray(ids, dtype=object)[detected]
    symbols = np.asarray(symbols, dtype=object)[detected]
    oracle_log = gate.normalize_log(oracle_raw)
    replay_log = gate.normalize_log(replay_raw)
    hvg = gate.oracle_hvgs(oracle_log, gate.N_HVG)
    oracle_pcs, replay_pcs, _, _, _, _ = gate.fit_shared_pca(
        oracle_log, replay_log, oracle_log, hvg
    )
    oracle_neighbors, _, oracle_adjacency = gate.neighbors(oracle_pcs)
    replay_neighbors, _, replay_adjacency = gate.neighbors(replay_pcs)
    oracle_labels = gate.leiden(oracle_adjacency)
    replay_labels = gate.leiden(replay_adjacency)
    oracle_categories = sorted(set(oracle_labels), key=lambda value: int(value))
    replay_categories = sorted(set(replay_labels), key=lambda value: int(value))
    oracle_index = {value: i for i, value in enumerate(oracle_categories)}
    replay_index = {value: i for i, value in enumerate(replay_categories)}
    contingency = np.zeros(
        (len(oracle_categories), len(replay_categories)), dtype=np.int64
    )
    for left, right in zip(oracle_labels, replay_labels):
        contingency[oracle_index[left], replay_index[right]] += 1
    row, column = scipy.optimize.linear_sum_assignment(contingency, maximize=True)
    matched = int(contingency[row, column].sum())
    replay_to_oracle = {
        replay_categories[j]: oracle_categories[i] for i, j in zip(row, column)
    }
    for j, label in enumerate(replay_categories):
        if label not in replay_to_oracle:
            replay_to_oracle[label] = oracle_categories[int(contingency[:, j].argmax())]
    remapped = np.asarray([replay_to_oracle[value] for value in replay_labels])
    label_stable = remapped == oracle_labels
    neighbor_jaccard = gate.neighbor_jaccards(oracle_neighbors, replay_neighbors)
    oracle_totals = np.asarray(oracle_raw.sum(axis=1)).ravel()
    replay_totals = np.asarray(replay_raw.sum(axis=1)).ravel()
    cell_total_relative_difference = np.divide(
        np.abs(replay_totals - oracle_totals),
        oracle_totals,
        out=np.zeros_like(oracle_totals),
        where=oracle_totals > 0,
    )
    cell_moved_mass = np.asarray(abs(oracle_raw - replay_raw).sum(axis=1)).ravel()
    cell_moved_mass = np.divide(
        cell_moved_mass,
        2 * oracle_totals,
        out=np.zeros_like(oracle_totals),
        where=oracle_totals > 0,
    )
    oracle_centroid, replay_centroid = nearest_centroid_labels(
        oracle_pcs, replay_pcs, oracle_labels
    )

    strata = []
    for label in oracle_categories:
        members = oracle_labels == label
        outside = ~members
        effect = (
            np.asarray(oracle_log[members].mean(axis=0)).ravel()
            - np.asarray(oracle_log[outside].mean(axis=0)).ravel()
        )
        top = np.argsort(effect, kind="stable")[-10:][::-1]
        row_values = contingency[oracle_index[label]]
        best_column = int(row_values.argmax())
        strata.append(
            {
                "oracle_cluster": label,
                "cells": int(members.sum()),
                "best_replay_cluster": replay_categories[best_column],
                "best_replay_overlap_fraction": float(
                    row_values[best_column] / members.sum()
                ),
                "mapped_label_stable_fraction": float(label_stable[members].mean()),
                "median_neighbor_jaccard": float(np.median(neighbor_jaccard[members])),
                "median_cell_total_relative_difference": float(
                    np.median(cell_total_relative_difference[members])
                ),
                "median_cell_moved_mass_fraction": float(
                    np.median(cell_moved_mass[members])
                ),
                "top_positive_markers": [
                    {"gene_id": str(ids[i]), "symbol": str(symbols[i])} for i in top
                ],
            }
        )

    return {
        "id": dataset,
        "role": config["role"],
        "oracle_cluster_count": len(oracle_categories),
        "replay_cluster_count": len(replay_categories),
        "hungarian_one_to_one_matched_fraction": matched / len(barcodes),
        "best_overlap_remapped_label_stable_fraction": float(label_stable.mean()),
        "oracle_vs_replay_fixed_centroid_label_agreement": float(
            np.mean(oracle_centroid == replay_centroid)
        ),
        "mismatched_cells": int((~label_stable).sum()),
        "mismatched_cell_median_total_relative_difference": float(
            np.median(cell_total_relative_difference[~label_stable])
            if np.any(~label_stable) else 0.0
        ),
        "stable_cell_median_total_relative_difference": float(
            np.median(cell_total_relative_difference[label_stable])
            if np.any(label_stable) else 0.0
        ),
        "mismatched_cell_median_moved_mass_fraction": float(
            np.median(cell_moved_mass[~label_stable]) if np.any(~label_stable) else 0.0
        ),
        "stable_cell_median_moved_mass_fraction": float(
            np.median(cell_moved_mass[label_stable]) if np.any(label_stable) else 0.0
        ),
        "contingency": {
            "oracle_labels": oracle_categories,
            "replay_labels": replay_categories,
            "counts": contingency.tolist(),
        },
        "oracle_cluster_strata": strata,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--primary-result", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest = yaml.safe_load(args.manifest.read_text())
    result = {
        "date": "2026-08-31",
        "status": "posthoc_failure_localization",
        "changes_frozen_verdicts": False,
        "primary_result_sha256": hashlib.sha256(
            args.primary_result.read_bytes()
        ).hexdigest(),
        "datasets": [],
    }
    for config in manifest["datasets"]:
        print(f"[{config['id']}] cluster diagnostics", flush=True)
        result["datasets"].append(analyze(config, args.project_root))
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
