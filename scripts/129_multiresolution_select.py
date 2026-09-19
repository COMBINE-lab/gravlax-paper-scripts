#!/usr/bin/env python3
"""Select stable downstream clustering settings from oracle matrices only."""
from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.util
import json
import platform
import sys
import warnings
from itertools import combinations
from pathlib import Path

import anndata
import numpy as np
import scanpy as sc
import sklearn
import yaml
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score


def load_gate() -> object:
    path = Path(__file__).with_name("126_downstream_fidelity.py")
    spec = importlib.util.spec_from_file_location("downstream_gate", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


gate = load_gate()


def leiden(adjacency, resolution: float, seed: int) -> np.ndarray:
    data = anndata.AnnData(X=np.zeros((adjacency.shape[0], 1), dtype=np.float32))
    data.obsp["connectivities"] = adjacency
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.leiden(
            data,
            adjacency=adjacency,
            resolution=resolution,
            random_state=seed,
            flavor="igraph",
            directed=False,
            n_iterations=2,
            key_added="cluster",
        )
    return data.obs["cluster"].cat.codes.to_numpy(dtype=np.int32)


def pairwise_ari(labels: list[np.ndarray]) -> list[float]:
    return [
        adjusted_rand_score(labels[left], labels[right])
        for left, right in combinations(range(len(labels)), 2)
    ]


def medoid(labels: list[np.ndarray]) -> tuple[np.ndarray, int, float]:
    scores = np.zeros(len(labels), dtype=np.float64)
    for left, right in combinations(range(len(labels)), 2):
        value = adjusted_rand_score(labels[left], labels[right])
        scores[left] += value
        scores[right] += value
    scores /= len(labels) - 1
    index = int(scores.argmax())
    return labels[index], index, float(scores[index])


def fit_oracle_pca(matrix, hvg: np.ndarray, components: int) -> np.ndarray:
    dense = matrix[:, hvg].toarray()
    mean = dense.mean(axis=0)
    std = dense.std(axis=0, ddof=1)
    std[std == 0] = 1.0
    maximum = min(components, dense.shape[0] - 1, dense.shape[1])
    model = PCA(n_components=maximum, svd_solver="arpack", random_state=gate.SEED)
    return model.fit_transform((dense - mean) / std)


def hash_order(dataset: str, replicate: int, barcodes: list[str]) -> np.ndarray:
    return np.asarray(
        sorted(
            range(len(barcodes)),
            key=lambda index: hashlib.sha256(
                f"multires|{dataset}|{replicate}|{barcodes[index]}".encode()
            ).digest(),
        ),
        dtype=np.int64,
    )


def marker_halves(dataset: str, barcodes: list[str]) -> np.ndarray:
    return np.asarray(
        [
            hashlib.sha256(f"multires-marker|{dataset}|{barcode}".encode()).digest()[0]
            & 1
            for barcode in barcodes
        ],
        dtype=np.int8,
    )


def marker_reproducibility(
    matrix, clusters: np.ndarray, halves: np.ndarray, top_k: int = 25
) -> dict[str, object]:
    per_cluster = []
    for cluster in sorted(set(map(int, clusters))):
        in_cluster = clusters == cluster
        left_inside = in_cluster & (halves == 0)
        left_outside = (~in_cluster) & (halves == 0)
        right_inside = in_cluster & (halves == 1)
        right_outside = (~in_cluster) & (halves == 1)
        counts = [
            int(left_inside.sum()),
            int(left_outside.sum()),
            int(right_inside.sum()),
            int(right_outside.sum()),
        ]
        if min(counts) < 10:
            per_cluster.append(
                {
                    "cluster": cluster,
                    "half_counts_inside_outside": counts,
                    "eligible": False,
                }
            )
            continue
        left_effect = (
            np.asarray(matrix[left_inside].mean(axis=0)).ravel()
            - np.asarray(matrix[left_outside].mean(axis=0)).ravel()
        )
        right_effect = (
            np.asarray(matrix[right_inside].mean(axis=0)).ravel()
            - np.asarray(matrix[right_outside].mean(axis=0)).ravel()
        )
        k = min(top_k, left_effect.size)
        left_top = set(np.argpartition(left_effect, -k)[-k:].tolist())
        right_top = set(np.argpartition(right_effect, -k)[-k:].tolist())
        per_cluster.append(
            {
                "cluster": cluster,
                "half_counts_inside_outside": counts,
                "eligible": True,
                "top25_positive_marker_jaccard": len(left_top & right_top)
                / len(left_top | right_top),
                "effect_pearson": gate.correlation(left_effect, right_effect),
            }
        )
    eligible = [item for item in per_cluster if item["eligible"]]
    all_groups = len(eligible) == len(per_cluster)
    return {
        "all_groups_evaluable": all_groups,
        "groups_evaluated": len(eligible),
        "groups_total": len(per_cluster),
        "median_top25_positive_marker_jaccard": (
            float(np.median([item["top25_positive_marker_jaccard"] for item in eligible]))
            if eligible
            else 0.0
        ),
        "median_effect_pearson": (
            float(np.median([item["effect_pearson"] for item in eligible]))
            if eligible
            else 0.0
        ),
        "per_cluster": per_cluster,
    }


def source_digests(directory: Path) -> list[dict[str, object]]:
    return [
        {
            "name": name,
            "bytes": (directory / name).stat().st_size,
            "sha256": gate.sha256(directory / name),
        }
        for name in ("barcodes.tsv", "features.tsv", "matrix.mtx")
    ]


def analyze_dataset(config: dict[str, object], manifest: dict[str, object], root: Path):
    dataset = str(config["id"])
    selection = manifest["selection"]
    dimensions = [int(value) for value in selection["dimensions"]]
    resolutions = [float(value) for value in selection["resolutions"]]
    full_seeds = [int(value) for value in selection["full_data_seeds"]]
    subsample_seeds = [int(value) for value in selection["subsample_seeds"]]
    replicate_count = int(selection["subsamples"])
    fraction = float(selection["subsample_fraction"])
    pca_components = int(selection["preprocessing"]["pca_fit_components"])
    thresholds = selection["eligibility"]
    oracle_dir = root / str(config["oracle_gene"])
    print(f"[{dataset}] loading oracle only", file=sys.stderr, flush=True)
    barcodes = gate.read_lines(oracle_dir / "barcodes.tsv")
    oracle_raw = gate.load_mtx_cells(oracle_dir / "matrix.mtx")
    detected = np.asarray(oracle_raw.getnnz(axis=0)).ravel() >= int(
        selection["preprocessing"]["minimum_oracle_cells_per_gene"]
    )
    oracle_log = gate.normalize_log(oracle_raw[:, detected])
    hvg = gate.oracle_hvgs(oracle_log, int(selection["preprocessing"]["hvg_count"]))
    full_pcs = fit_oracle_pca(oracle_log, hvg, pca_components)
    halves = marker_halves(dataset, barcodes)

    full_labels: dict[tuple[int, float], list[np.ndarray]] = {}
    full_medoids: dict[tuple[int, float], np.ndarray] = {}
    full_medoid_seed: dict[tuple[int, float], int] = {}
    full_medoid_centrality: dict[tuple[int, float], float] = {}
    for dimension in dimensions:
        print(f"[{dataset}] full d={dimension}", file=sys.stderr, flush=True)
        _, _, graph = gate.neighbors(full_pcs[:, :dimension])
        for resolution in resolutions:
            labels = [leiden(graph, resolution, seed) for seed in full_seeds]
            center, index, centrality = medoid(labels)
            key = (dimension, resolution)
            full_labels[key] = labels
            full_medoids[key] = center
            full_medoid_seed[key] = full_seeds[index]
            full_medoid_centrality[key] = centrality

    subsample_indices: list[np.ndarray] = []
    subsample_labels: list[dict[tuple[int, float], list[np.ndarray]]] = []
    count = int(np.floor(fraction * len(barcodes)))
    for replicate in range(replicate_count):
        indices = np.sort(hash_order(dataset, replicate, barcodes)[:count])
        subsample_indices.append(indices)
        pcs = fit_oracle_pca(oracle_log[indices], hvg, pca_components)
        labels_by_candidate: dict[tuple[int, float], list[np.ndarray]] = {}
        for dimension in dimensions:
            print(
                f"[{dataset}] subsample={replicate} d={dimension}",
                file=sys.stderr,
                flush=True,
            )
            _, _, graph = gate.neighbors(pcs[:, :dimension])
            for resolution in resolutions:
                labels_by_candidate[(dimension, resolution)] = [
                    leiden(graph, resolution, seed) for seed in subsample_seeds
                ]
        subsample_labels.append(labels_by_candidate)
        del pcs
        gc.collect()

    rows = []
    for dimension in dimensions:
        for resolution in resolutions:
            key = (dimension, resolution)
            optimizer_values = pairwise_ari(full_labels[key])
            subsample_medoids = []
            for replicate in range(replicate_count):
                labels = subsample_labels[replicate][key]
                optimizer_values.extend(pairwise_ari(labels))
                subsample_medoids.append(medoid(labels)[0])
            subsample_values = []
            for left, right in combinations(range(replicate_count), 2):
                _, left_positions, right_positions = np.intersect1d(
                    subsample_indices[left],
                    subsample_indices[right],
                    assume_unique=True,
                    return_indices=True,
                )
                subsample_values.append(
                    adjusted_rand_score(
                        subsample_medoids[left][left_positions],
                        subsample_medoids[right][right_positions],
                    )
                )
            center = full_medoids[key]
            cluster_sizes = np.bincount(center)
            adjacent_dimensions = []
            d_index = dimensions.index(dimension)
            for neighbor in (d_index - 1, d_index + 1):
                if 0 <= neighbor < len(dimensions):
                    adjacent_dimensions.append(
                        adjusted_rand_score(
                            center, full_medoids[(dimensions[neighbor], resolution)]
                        )
                    )
            adjacent_resolutions = []
            r_index = resolutions.index(resolution)
            for neighbor in (r_index - 1, r_index + 1):
                if 0 <= neighbor < len(resolutions):
                    adjacent_resolutions.append(
                        adjusted_rand_score(
                            center, full_medoids[(dimension, resolutions[neighbor])]
                        )
                    )
            markers = marker_reproducibility(oracle_log, center, halves)
            optimizer_median = float(np.median(optimizer_values))
            subsample_median = float(np.median(subsample_values))
            adjacent_dimension_best = float(max(adjacent_dimensions))
            adjacent_resolution_best = float(max(adjacent_resolutions))
            robustness = min(
                optimizer_median,
                subsample_median,
                adjacent_dimension_best,
                adjacent_resolution_best,
            )
            checks = {
                "cluster_count_min": len(cluster_sizes)
                >= int(thresholds["oracle_medoid_cluster_count_min"]),
                "cluster_count_max": len(cluster_sizes)
                <= int(thresholds["oracle_medoid_cluster_count_max"]),
                "minimum_cluster_cells": int(cluster_sizes.min())
                >= int(thresholds["minimum_cluster_cells"]),
                "maximum_cluster_fraction": float(cluster_sizes.max() / len(center))
                <= float(thresholds["maximum_cluster_fraction"]),
                "optimizer_seed_stability": optimizer_median
                >= float(thresholds["optimizer_seed_ari_median_min"]),
                "cell_subsample_stability": subsample_median
                >= float(thresholds["cell_subsample_ari_median_min"]),
                "adjacent_dimension_stability": adjacent_dimension_best
                >= float(thresholds["best_adjacent_dimension_medoid_ari_min"]),
                "adjacent_resolution_stability": adjacent_resolution_best
                >= float(thresholds["best_adjacent_resolution_medoid_ari_min"]),
                "all_marker_groups_evaluable": bool(markers["all_groups_evaluable"]),
                "marker_top25_reproducibility": float(
                    markers["median_top25_positive_marker_jaccard"]
                )
                >= float(thresholds["split_half_marker_top25_jaccard_median_min"]),
                "marker_effect_reproducibility": float(markers["median_effect_pearson"])
                >= float(thresholds["split_half_marker_effect_pearson_median_min"]),
            }
            rows.append(
                {
                    "dimension": dimension,
                    "resolution": resolution,
                    "oracle_medoid_seed": full_medoid_seed[key],
                    "oracle_medoid_centrality": full_medoid_centrality[key],
                    "oracle_medoid_cluster_count": len(cluster_sizes),
                    "minimum_cluster_cells": int(cluster_sizes.min()),
                    "maximum_cluster_fraction": float(cluster_sizes.max() / len(center)),
                    "optimizer_seed_ari_median": optimizer_median,
                    "cell_subsample_ari_median": subsample_median,
                    "best_adjacent_dimension_medoid_ari": adjacent_dimension_best,
                    "best_adjacent_resolution_medoid_ari": adjacent_resolution_best,
                    "oracle_robustness_score": robustness,
                    "split_half_markers": markers,
                    "eligibility_checks": checks,
                    "eligible": all(checks.values()),
                }
            )

    eligible = [row for row in rows if row["eligible"]]
    if eligible:
        maximum_clusters = max(row["oracle_medoid_cluster_count"] for row in eligible)
        finalists = [
            row
            for row in eligible
            if row["oracle_medoid_cluster_count"] >= maximum_clusters - 1
        ]
        finalists.sort(
            key=lambda row: (
                -row["oracle_robustness_score"],
                -row["split_half_markers"]["median_top25_positive_marker_jaccard"],
                row["dimension"],
                row["resolution"],
            )
        )
        chosen = finalists[0]
        locked = {
            "status": "stable_partition_locked",
            "dimension": chosen["dimension"],
            "resolution": chosen["resolution"],
            "oracle_medoid_cluster_count": chosen["oracle_medoid_cluster_count"],
            "oracle_robustness_score": chosen["oracle_robustness_score"],
            "maximum_eligible_cluster_count": maximum_clusters,
        }
    else:
        locked = {"status": "no_stable_partition"}

    result = {
        "id": dataset,
        "role": config["role"],
        "called_cells": len(barcodes),
        "retained_genes": int(detected.sum()),
        "hvg_count": int(hvg.sum()),
        "oracle_source_digests": source_digests(oracle_dir),
        "subsample_cell_counts": [len(indices) for indices in subsample_indices],
        "locked_selection": locked,
        "eligible_candidate_count": len(eligible),
        "candidates": rows,
    }
    del oracle_raw, oracle_log, full_pcs
    gc.collect()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest = yaml.safe_load(args.manifest.read_text())
    result = {
        "status": "complete",
        "analysis": "locked_postexploratory_oracle_only_multiresolution_selection",
        "manifest": str(args.manifest),
        "manifest_sha256": gate.sha256(args.manifest),
        "replay_files_opened": False,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scanpy": sc.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "datasets": [
            analyze_dataset(config, manifest, args.project_root)
            for config in manifest["datasets"]
        ],
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
