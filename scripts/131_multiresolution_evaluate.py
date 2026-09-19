#!/usr/bin/env python3
"""Evaluate replay only at oracle-locked multiresolution settings."""
from __future__ import annotations

import argparse
import gc
import importlib.util
import json
import platform
import sys
import warnings
from itertools import combinations, product
from pathlib import Path

import anndata
import numpy as np
import scanpy as sc
import scipy.sparse as sp
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


def json_default(value: object) -> object:
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


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


def distribution(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    return {
        "p10": float(np.quantile(array, 0.10)),
        "median": float(np.median(array)),
        "mean": float(array.mean()),
        "p90": float(np.quantile(array, 0.90)),
    }


def pairwise(labels: list[np.ndarray]) -> list[float]:
    return [
        adjusted_rand_score(labels[left], labels[right])
        for left, right in combinations(range(len(labels)), 2)
    ]


def cross(left: list[np.ndarray], right: list[np.ndarray]) -> list[float]:
    return [adjusted_rand_score(a, b) for a, b in product(left, right)]


def medoid(labels: list[np.ndarray], seeds: list[int]) -> tuple[np.ndarray, int, float]:
    scores = np.zeros(len(labels), dtype=np.float64)
    for left, right in combinations(range(len(labels)), 2):
        value = adjusted_rand_score(labels[left], labels[right])
        scores[left] += value
        scores[right] += value
    scores /= len(labels) - 1
    index = int(scores.argmax())
    return labels[index], seeds[index], float(scores[index])


def fit_shared(
    oracle_log: sp.csr_matrix,
    replay_log: sp.csr_matrix,
    hvg: np.ndarray,
    components: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, PCA]:
    oracle = oracle_log[:, hvg].toarray()
    replay = replay_log[:, hvg].toarray()
    mean = oracle.mean(axis=0)
    std = oracle.std(axis=0, ddof=1)
    std[std == 0] = 1.0
    maximum = min(components, oracle.shape[0] - 1, oracle.shape[1])
    model = PCA(n_components=maximum, svd_solver="arpack", random_state=gate.SEED)
    oracle_pcs = model.fit_transform((oracle - mean) / std)
    replay_pcs = model.transform((replay - mean) / std)
    return oracle_pcs, replay_pcs, mean, std, model


def project_control(
    control_log: sp.csr_matrix,
    hvg: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    model: PCA,
) -> np.ndarray:
    dense = control_log[:, hvg].toarray()
    return model.transform((dense - mean) / std)


def source_digests(directory: Path) -> list[dict[str, object]]:
    return [
        {
            "name": name,
            "bytes": (directory / name).stat().st_size,
            "sha256": gate.sha256(directory / name),
        }
        for name in ("barcodes.tsv", "features.tsv", "matrix.mtx")
    ]


def analyze_dataset(
    config: dict[str, object],
    locked: dict[str, object],
    manifest: dict[str, object],
    root: Path,
) -> dict[str, object]:
    dataset = str(config["id"])
    selection = locked["locked_selection"]
    if selection["status"] != "stable_partition_locked":
        return {
            "id": dataset,
            "role": config["role"],
            "status": "not_evaluated_no_stable_oracle_partition",
        }
    dimension = int(selection["dimension"])
    resolution = float(selection["resolution"])
    seeds = [int(value) for value in manifest["evaluation"]["ensemble_seeds"]]
    pca_components = int(manifest["selection"]["preprocessing"]["pca_fit_components"])
    print(
        f"[{dataset}] replay evaluation d={dimension} r={resolution}",
        file=sys.stderr,
        flush=True,
    )

    oracle_dir = root / str(config["oracle_gene"])
    replay_dir = root / str(config["replay_gene"])
    barcodes = gate.read_lines(oracle_dir / "barcodes.tsv")
    gate.require_same_features(oracle_dir / "features.tsv", replay_dir / "features.tsv")
    replay_columns = gate.replay_column_indices(
        replay_dir / "barcodes.tsv", barcodes
    )
    oracle_raw = gate.load_mtx_cells(oracle_dir / "matrix.mtx")
    replay_raw = gate.load_mtx_cells(replay_dir / "matrix.mtx", replay_columns)
    detected = np.asarray(oracle_raw.getnnz(axis=0)).ravel() >= int(
        manifest["selection"]["preprocessing"]["minimum_oracle_cells_per_gene"]
    )
    oracle_raw = oracle_raw[:, detected]
    replay_raw = replay_raw[:, detected]
    oracle_log = gate.normalize_log(oracle_raw)
    replay_log = gate.normalize_log(replay_raw)
    hvg = gate.oracle_hvgs(
        oracle_log, int(manifest["selection"]["preprocessing"]["hvg_count"])
    )
    oracle_pcs, replay_pcs, mean, std, model = fit_shared(
        oracle_log, replay_log, hvg, pca_components
    )
    oracle_knn, _, oracle_graph = gate.neighbors(oracle_pcs[:, :dimension])
    replay_knn, _, replay_graph = gate.neighbors(replay_pcs[:, :dimension])
    oracle_labels = [leiden(oracle_graph, resolution, seed) for seed in seeds]
    replay_labels = [leiden(replay_graph, resolution, seed) for seed in seeds]
    oracle_medoid, oracle_medoid_seed, oracle_centrality = medoid(
        oracle_labels, seeds
    )
    replay_medoid, replay_medoid_seed, replay_centrality = medoid(
        replay_labels, seeds
    )

    control_sources, control_selected = gate.choose_control_sources(
        dataset, barcodes, oracle_medoid
    )
    control_raw = oracle_raw[control_sources].tocsr()
    control_log = gate.normalize_log(control_raw)
    control_pcs = project_control(control_log, hvg, mean, std, model)
    control_knn, _, control_graph = gate.neighbors(control_pcs[:, :dimension])
    control_labels = [leiden(control_graph, resolution, seed) for seed in seeds]
    control_medoid, control_medoid_seed, control_centrality = medoid(
        control_labels, seeds
    )

    within_oracle = pairwise(oracle_labels)
    within_replay = pairwise(replay_labels)
    replay_cross = cross(oracle_labels, replay_labels)
    control_cross = cross(oracle_labels, control_labels)
    within_oracle_median = float(np.median(within_oracle))
    within_replay_median = float(np.median(within_replay))
    replay_cross_median = float(np.median(replay_cross))
    control_cross_median = float(np.median(control_cross))
    excess_loss = (
        (within_oracle_median + within_replay_median) / 2 - replay_cross_median
    )
    replay_metrics = gate.arm_metrics(
        oracle_medoid,
        replay_medoid,
        oracle_knn,
        replay_knn,
        oracle_log,
        replay_log,
        oracle_raw,
        replay_raw,
    )
    control_metrics = gate.arm_metrics(
        oracle_medoid,
        control_medoid,
        oracle_knn,
        control_knn,
        oracle_log,
        control_log,
        oracle_raw,
        control_raw,
    )
    partition_gates = manifest["evaluation"]["partition_gates"]
    fixed_gates = manifest["evaluation"]["fixed_group_gates"]
    sensitivity = manifest["evaluation"]["sensitivity_control"]
    checks = {
        "medoid_cross_ari": replay_metrics["ari"]
        >= float(partition_gates["medoid_cross_ari_min"]),
        "cross_excess_loss_vs_mean_within": excess_loss
        <= float(partition_gates["cross_excess_loss_vs_mean_within_max"]),
        "mean_marker_top25_jaccard": replay_metrics["markers"]["mean_top25_jaccard"]
        >= float(fixed_gates["mean_marker_top25_jaccard_min"]),
        "median_marker_effect_pearson": replay_metrics["markers"][
            "median_effect_pearson"
        ]
        >= float(fixed_gates["median_marker_effect_pearson_min"]),
        "pseudobulk_pearson": replay_metrics["pseudobulk_pearson"]
        >= float(fixed_gates["pseudobulk_pearson_min"]),
        "pseudobulk_spearman": replay_metrics["pseudobulk_spearman"]
        >= float(fixed_gates["pseudobulk_spearman_min"]),
        "replay_minus_control_cross_ari": replay_cross_median - control_cross_median
        >= float(sensitivity["replay_minus_control_cross_ari_min"]),
        "replay_minus_control_neighbor_jaccard": replay_metrics[
            "neighbor_jaccard_median"
        ]
        - control_metrics["neighbor_jaccard_median"]
        >= float(sensitivity["replay_minus_control_neighbor_jaccard_min"]),
    }
    component_verdicts = {
        "partition_retention": (
            "PASS"
            if checks["medoid_cross_ari"]
            and checks["cross_excess_loss_vs_mean_within"]
            else "FAIL"
        ),
        "marker_retention": (
            "PASS"
            if checks["mean_marker_top25_jaccard"]
            and checks["median_marker_effect_pearson"]
            else "FAIL"
        ),
        "pseudobulk_retention": (
            "PASS"
            if checks["pseudobulk_pearson"] and checks["pseudobulk_spearman"]
            else "FAIL"
        ),
        "control_sensitivity": (
            "PASS"
            if checks["replay_minus_control_cross_ari"]
            and checks["replay_minus_control_neighbor_jaccard"]
            else "FAIL"
        ),
        "composite": "PASS" if all(checks.values()) else "FAIL",
    }
    result = {
        "id": dataset,
        "role": config["role"],
        "status": "complete",
        "locked_dimension": dimension,
        "locked_resolution": resolution,
        "locked_oracle_cluster_count": selection["oracle_medoid_cluster_count"],
        "called_cells": len(barcodes),
        "retained_genes": int(detected.sum()),
        "hvg_count": int(hvg.sum()),
        "oracle_source_digests": source_digests(oracle_dir),
        "replay_source_digests": source_digests(replay_dir),
        "ensemble": {
            "seeds": seeds,
            "within_oracle_ari": distribution(within_oracle),
            "within_replay_ari": distribution(within_replay),
            "replay_cross_ari": distribution(replay_cross),
            "control_cross_ari": distribution(control_cross),
            "cross_excess_loss_vs_mean_within": excess_loss,
            "replay_minus_control_cross_ari": replay_cross_median
            - control_cross_median,
            "oracle_medoid_seed": oracle_medoid_seed,
            "replay_medoid_seed": replay_medoid_seed,
            "control_medoid_seed": control_medoid_seed,
            "oracle_medoid_centrality": oracle_centrality,
            "replay_medoid_centrality": replay_centrality,
            "control_medoid_centrality": control_centrality,
        },
        "replay": replay_metrics,
        "control": control_metrics,
        "control_selected_cells": int(control_selected.sum()),
        "replay_minus_control_neighbor_jaccard": replay_metrics[
            "neighbor_jaccard_median"
        ]
        - control_metrics["neighbor_jaccard_median"],
        "gate": {
            "checks": checks,
            "component_verdicts": component_verdicts,
            "verdict": component_verdicts["composite"],
        },
    }
    del oracle_raw, replay_raw, control_raw, oracle_log, replay_log, control_log
    gc.collect()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest = yaml.safe_load(args.manifest.read_text())
    selection = json.loads(args.selection.read_text())
    manifest_sha256 = gate.sha256(args.manifest)
    if selection["status"] != "complete":
        raise SystemExit("oracle selection is not complete")
    if selection["manifest_sha256"] != manifest_sha256:
        raise SystemExit("oracle selection manifest digest does not match evaluation manifest")
    locked_by_id = {item["id"]: item for item in selection["datasets"]}
    datasets = [
        analyze_dataset(config, locked_by_id[config["id"]], manifest, args.project_root)
        for config in manifest["datasets"]
    ]
    human = [
        item
        for item in datasets
        if item["role"] == "confirmatory_human" and item["status"] == "complete"
    ]
    mouse = [
        item
        for item in datasets
        if item["role"] == "cross_chemistry_stress_test"
        and item["status"] == "complete"
    ]
    result = {
        "status": "complete",
        "analysis": "locked_postexploratory_multiresolution_replay_evaluation",
        "manifest": str(args.manifest),
        "manifest_sha256": manifest_sha256,
        "oracle_selection": str(args.selection),
        "oracle_selection_sha256": gate.sha256(args.selection),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scanpy": sc.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "datasets": datasets,
        "summary": {
            "human_partition_retention_verdict": (
                "PASS"
                if len(human) == 3
                and all(
                    item["gate"]["component_verdicts"]["partition_retention"]
                    == "PASS"
                    for item in human
                )
                else "FAIL"
            ),
            "human_marker_retention_verdict": (
                "PASS"
                if len(human) == 3
                and all(
                    item["gate"]["component_verdicts"]["marker_retention"]
                    == "PASS"
                    for item in human
                )
                else "FAIL"
            ),
            "human_pseudobulk_retention_verdict": (
                "PASS"
                if len(human) == 3
                and all(
                    item["gate"]["component_verdicts"]["pseudobulk_retention"]
                    == "PASS"
                    for item in human
                )
                else "FAIL"
            ),
            "human_control_sensitivity_verdict": (
                "PASS"
                if len(human) == 3
                and all(
                    item["gate"]["component_verdicts"]["control_sensitivity"]
                    == "PASS"
                    for item in human
                )
                else "FAIL"
            ),
            "human_composite_verdict": (
                "PASS"
                if len(human) == 3
                and all(item["gate"]["verdict"] == "PASS" for item in human)
                else "FAIL"
            ),
            "mouse_5_prime_component_verdicts": (
                mouse[0]["gate"]["component_verdicts"]
                if len(mouse) == 1
                else "NOT_EVALUATED"
            ),
            "frozen_resolution_0_8_verdict_changed": False,
        },
    }
    args.out.write_text(json.dumps(result, indent=2, default=json_default) + "\n")


if __name__ == "__main__":
    main()
