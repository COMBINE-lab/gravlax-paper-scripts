#!/usr/bin/env python3
"""Run the preregistered Gene and velocity downstream-fidelity gate."""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import sys
import time
from pathlib import Path

import anndata
import numpy as np
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import scipy.stats
import sklearn
import yaml
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.neighbors import NearestNeighbors


SEED = 1729
TARGET_SUM = 10_000.0
N_HVG = 2_000
N_PCS = 30
N_NEIGHBORS = 15
LEIDEN_RESOLUTION = 0.8


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_lines(path: Path) -> list[str]:
    return [line.rstrip("\n") for line in path.open()]


def feature_ids(path: Path) -> list[str]:
    return [line.split("\t", 1)[0] for line in path.open()]


def require_same_features(left: Path, right: Path) -> list[str]:
    a = feature_ids(left)
    b = feature_ids(right)
    if a != b:
        raise SystemExit(f"feature order differs: {left} versus {right}")
    return a


def replay_column_indices(raw_barcodes: Path, selected: list[str]) -> np.ndarray:
    wanted = {barcode: i for i, barcode in enumerate(selected)}
    result = np.full(len(selected), -1, dtype=np.int64)
    for raw_index, line in enumerate(raw_barcodes.open()):
        local = wanted.get(line.rstrip("\n"))
        if local is not None:
            result[local] = raw_index
    if np.any(result < 0):
        missing = [selected[i] for i in np.flatnonzero(result < 0)[:5]]
        raise SystemExit(f"selected barcodes absent from {raw_barcodes}: {missing}")
    return result


def load_mtx_cells(path: Path, columns: np.ndarray | None = None) -> sp.csr_matrix:
    matrix = scipy.io.mmread(str(path))
    if not sp.issparse(matrix):
        matrix = sp.coo_matrix(matrix)
    matrix = matrix.tocsc()
    if columns is not None:
        matrix = matrix[:, columns]
    result = matrix.T.tocsr().astype(np.float64)
    result.eliminate_zeros()
    return result


def normalize_log(matrix: sp.csr_matrix) -> sp.csr_matrix:
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    scale = np.divide(TARGET_SUM, totals, out=np.zeros_like(totals), where=totals > 0)
    result = sp.diags(scale) @ matrix
    result = result.tocsr()
    result.data = np.log1p(result.data)
    return result


def oracle_hvgs(log_matrix: sp.csr_matrix, n_top: int) -> np.ndarray:
    data = anndata.AnnData(X=log_matrix.copy())
    sc.pp.highly_variable_genes(
        data,
        flavor="seurat",
        n_top_genes=min(n_top, data.n_vars),
        inplace=True,
    )
    mask = np.asarray(data.var["highly_variable"], dtype=bool)
    if mask.sum() < min(50, data.n_vars):
        raise SystemExit(f"too few HVGs selected: {mask.sum()}")
    return mask


def fit_shared_pca(
    oracle: sp.csr_matrix,
    replay: sp.csr_matrix,
    control: sp.csr_matrix,
    hvg: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, PCA]:
    oracle_dense = oracle[:, hvg].toarray()
    replay_dense = replay[:, hvg].toarray()
    control_dense = control[:, hvg].toarray()
    mean = oracle_dense.mean(axis=0)
    std = oracle_dense.std(axis=0, ddof=1)
    std[std == 0] = 1.0
    oracle_scaled = (oracle_dense - mean) / std
    replay_scaled = (replay_dense - mean) / std
    control_scaled = (control_dense - mean) / std
    n_components = min(N_PCS, oracle_scaled.shape[0] - 1, oracle_scaled.shape[1])
    pca = PCA(n_components=n_components, svd_solver="arpack", random_state=SEED)
    oracle_pcs = pca.fit_transform(oracle_scaled)
    replay_pcs = pca.transform(replay_scaled)
    control_pcs = pca.transform(control_scaled)
    return oracle_pcs, replay_pcs, control_pcs, mean, std, pca


def neighbors(pcs: np.ndarray) -> tuple[np.ndarray, np.ndarray, sp.csr_matrix]:
    k = min(N_NEIGHBORS, pcs.shape[0] - 1)
    model = NearestNeighbors(n_neighbors=k + 1, metric="euclidean", n_jobs=1)
    model.fit(pcs)
    distances, indices = model.kneighbors(pcs)
    indices = indices[:, 1:]
    distances = distances[:, 1:]
    rows = np.repeat(np.arange(pcs.shape[0]), k)
    weights = 1.0 / (1.0 + distances.ravel())
    adjacency = sp.csr_matrix(
        (weights, (rows, indices.ravel())), shape=(pcs.shape[0], pcs.shape[0])
    )
    adjacency = adjacency.maximum(adjacency.T)
    return indices, distances, adjacency


def leiden(adjacency: sp.csr_matrix) -> np.ndarray:
    data = anndata.AnnData(X=sp.csr_matrix((adjacency.shape[0], 1)))
    data.obsp["connectivities"] = adjacency
    sc.tl.leiden(
        data,
        adjacency=adjacency,
        resolution=LEIDEN_RESOLUTION,
        random_state=SEED,
        flavor="igraph",
        directed=False,
        n_iterations=2,
        key_added="cluster",
    )
    return np.asarray(data.obs["cluster"].astype(str))


def neighbor_jaccards(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    values = np.zeros(left.shape[0], dtype=np.float64)
    for i in range(left.shape[0]):
        a = set(map(int, left[i]))
        b = set(map(int, right[i]))
        values[i] = len(a & b) / len(a | b)
    return values


def correlation(left: np.ndarray, right: np.ndarray, method: str = "pearson") -> float:
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    if np.all(left == left[0]) or np.all(right == right[0]):
        return 1.0 if np.array_equal(left, right) else 0.0
    if method == "spearman":
        return float(scipy.stats.spearmanr(left, right).statistic)
    return float(np.corrcoef(left, right)[0, 1])


def marker_metrics(
    oracle: sp.csr_matrix,
    other: sp.csr_matrix,
    clusters: np.ndarray,
    top_k: int = 25,
) -> dict[str, object]:
    per_cluster = []
    for cluster in sorted(set(clusters), key=lambda value: int(value)):
        inside = clusters == cluster
        outside = ~inside
        effect_oracle = (
            np.asarray(oracle[inside].mean(axis=0)).ravel()
            - np.asarray(oracle[outside].mean(axis=0)).ravel()
        )
        effect_other = (
            np.asarray(other[inside].mean(axis=0)).ravel()
            - np.asarray(other[outside].mean(axis=0)).ravel()
        )
        k = min(top_k, effect_oracle.size)
        top_oracle = set(np.argpartition(np.abs(effect_oracle), -k)[-k:].tolist())
        top_other = set(np.argpartition(np.abs(effect_other), -k)[-k:].tolist())
        jaccard = len(top_oracle & top_other) / len(top_oracle | top_other)
        per_cluster.append(
            {
                "cluster": cluster,
                "cells": int(inside.sum()),
                "top25_jaccard": jaccard,
                "effect_pearson": correlation(effect_oracle, effect_other),
            }
        )
    return {
        "per_cluster": per_cluster,
        "mean_top25_jaccard": float(
            np.mean([item["top25_jaccard"] for item in per_cluster])
        ),
        "median_effect_pearson": float(
            np.median([item["effect_pearson"] for item in per_cluster])
        ),
    }


def choose_control_sources(
    dataset: str, barcodes: list[str], clusters: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    n_select = int(math.floor(0.20 * len(barcodes)))
    ordered = sorted(
        range(len(barcodes)),
        key=lambda i: hashlib.sha256(f"{dataset}|{barcodes[i]}".encode()).digest(),
    )
    selected = np.asarray(ordered[:n_select], dtype=np.int64)
    labels = sorted(set(clusters), key=lambda value: int(value))
    members = {label: np.flatnonzero(clusters == label) for label in labels}
    next_label = {label: labels[(i + 1) % len(labels)] for i, label in enumerate(labels)}
    sources = np.arange(len(barcodes), dtype=np.int64)
    used = {label: 0 for label in labels}
    for cell in selected:
        target = next_label[clusters[cell]]
        donors = members[target]
        sources[cell] = donors[used[target] % len(donors)]
        used[target] += 1
    return sources, selected


def arm_metrics(
    oracle_clusters: np.ndarray,
    other_clusters: np.ndarray,
    oracle_neighbors: np.ndarray,
    other_neighbors: np.ndarray,
    oracle_log: sp.csr_matrix,
    other_log: sp.csr_matrix,
    oracle_raw: sp.csr_matrix,
    other_raw: sp.csr_matrix,
) -> dict[str, object]:
    jaccard = neighbor_jaccards(oracle_neighbors, other_neighbors)
    marker = marker_metrics(oracle_log, other_log, oracle_clusters)
    oracle_bulk = np.asarray(oracle_raw.sum(axis=0)).ravel()
    other_bulk = np.asarray(other_raw.sum(axis=0)).ravel()
    return {
        "ari": float(adjusted_rand_score(oracle_clusters, other_clusters)),
        "nmi": float(normalized_mutual_info_score(oracle_clusters, other_clusters)),
        "neighbor_jaccard_mean": float(jaccard.mean()),
        "neighbor_jaccard_median": float(np.median(jaccard)),
        "neighbor_jaccard_p10": float(np.quantile(jaccard, 0.10)),
        "markers": marker,
        "pseudobulk_pearson": correlation(oracle_bulk, other_bulk),
        "pseudobulk_spearman": correlation(oracle_bulk, other_bulk, "spearman"),
    }


def gene_gate(replay: dict[str, object], control: dict[str, object]) -> dict[str, object]:
    checks = {
        "ari_at_least_0_95": replay["ari"] >= 0.95,
        "nmi_at_least_0_95": replay["nmi"] >= 0.95,
        "median_neighbor_jaccard_at_least_0_80": replay["neighbor_jaccard_median"] >= 0.80,
        "mean_marker_top25_jaccard_at_least_0_80": replay["markers"]["mean_top25_jaccard"] >= 0.80,
        "median_marker_effect_pearson_at_least_0_98": replay["markers"]["median_effect_pearson"] >= 0.98,
        "pseudobulk_pearson_at_least_0_999": replay["pseudobulk_pearson"] >= 0.999,
        "pseudobulk_spearman_at_least_0_999": replay["pseudobulk_spearman"] >= 0.999,
        "ari_exceeds_control_by_0_05": replay["ari"] - control["ari"] >= 0.05,
        "neighbor_jaccard_exceeds_control_by_0_05": (
            replay["neighbor_jaccard_median"] - control["neighbor_jaccard_median"] >= 0.05
        ),
    }
    return {"checks": checks, "verdict": "PASS" if all(checks.values()) else "FAIL"}


def smoothing_matrix(neighbor_indices: np.ndarray) -> sp.csr_matrix:
    n, k = neighbor_indices.shape
    rows = np.repeat(np.arange(n), k + 1)
    columns = np.concatenate(
        [np.column_stack([np.arange(n), neighbor_indices]).ravel()]
    )
    values = np.full(rows.size, 1.0 / (k + 1), dtype=np.float64)
    return sp.csr_matrix((values, (rows, columns)), shape=(n, n))


def velocity_cosines(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    numerator = np.einsum("ij,ij->i", left, right)
    denominator = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)


def transition_top_sets(
    pcs: np.ndarray, velocity: np.ndarray, neighbor_indices: np.ndarray, top_k: int = 3
) -> list[set[int]]:
    result = []
    for i, candidates in enumerate(neighbor_indices):
        delta = pcs[candidates] - pcs[i]
        denominator = np.linalg.norm(delta, axis=1) * np.linalg.norm(velocity[i])
        scores = np.divide(
            delta @ velocity[i], denominator, out=np.full(len(candidates), -np.inf), where=denominator > 0
        )
        positive = np.flatnonzero(scores > 0)
        if positive.size:
            order = positive[np.argsort(scores[positive], kind="stable")[-top_k:]]
            result.append(set(map(int, candidates[order])))
        else:
            result.append(set())
    return result


def velocity_metrics(
    oracle_spliced: sp.csr_matrix,
    oracle_unspliced: sp.csr_matrix,
    replay_spliced: sp.csr_matrix,
    replay_unspliced: sp.csr_matrix,
    hvg: np.ndarray,
    gene_std: np.ndarray,
    pca: PCA,
    oracle_pcs: np.ndarray,
    oracle_neighbors: np.ndarray,
) -> dict[str, object]:
    eligible = (
        np.asarray(oracle_spliced.getnnz(axis=0)).ravel() >= 20
    ) & (
        np.asarray(oracle_unspliced.getnnz(axis=0)).ravel() >= 10
    ) & hvg
    hvg_positions = np.flatnonzero(hvg)
    eligible_global = np.flatnonzero(eligible)
    position_of_global = {value: i for i, value in enumerate(hvg_positions)}
    component_positions = np.asarray(
        [position_of_global[value] for value in eligible_global], dtype=np.int64
    )
    if eligible_global.size < 20:
        raise SystemExit(f"too few velocity genes: {eligible_global.size}")

    def normalize_layers(
        spliced: sp.csr_matrix, unspliced: sp.csr_matrix
    ) -> tuple[sp.csr_matrix, sp.csr_matrix]:
        totals = np.asarray((spliced + unspliced).sum(axis=1)).ravel()
        scale = np.divide(TARGET_SUM, totals, out=np.zeros_like(totals), where=totals > 0)
        diagonal = sp.diags(scale)
        return (diagonal @ spliced).tocsr(), (diagonal @ unspliced).tocsr()

    oracle_s, oracle_u = normalize_layers(oracle_spliced, oracle_unspliced)
    replay_s, replay_u = normalize_layers(replay_spliced, replay_unspliced)
    smoother = smoothing_matrix(oracle_neighbors)
    oracle_s = smoother @ oracle_s[:, eligible_global]
    oracle_u = smoother @ oracle_u[:, eligible_global]
    replay_s = smoother @ replay_s[:, eligible_global]
    replay_u = smoother @ replay_u[:, eligible_global]
    sum_s = np.asarray(oracle_s.sum(axis=0)).ravel()
    gamma = np.divide(
        np.asarray(oracle_u.sum(axis=0)).ravel(),
        sum_s,
        out=np.zeros_like(sum_s),
        where=sum_s > 0,
    )
    oracle_residual = oracle_u.toarray() - oracle_s.toarray() * gamma
    replay_residual = replay_u.toarray() - replay_s.toarray() * gamma
    loadings = pca.components_[:, component_positions]
    scale = gene_std[component_positions]
    oracle_velocity = (oracle_residual / scale) @ loadings.T
    replay_velocity = (replay_residual / scale) @ loadings.T
    cosines = velocity_cosines(oracle_velocity, replay_velocity)
    oracle_transitions = transition_top_sets(
        oracle_pcs, oracle_velocity, oracle_neighbors, top_k=3
    )
    replay_transitions = transition_top_sets(
        oracle_pcs, replay_velocity, oracle_neighbors, top_k=3
    )
    transition_jaccard = []
    both_empty = 0
    for left, right in zip(oracle_transitions, replay_transitions):
        if not left and not right:
            both_empty += 1
            transition_jaccard.append(1.0)
        else:
            transition_jaccard.append(len(left & right) / len(left | right))
    transition_jaccard = np.asarray(transition_jaccard)
    checks = {
        "median_cosine_at_least_0_90": float(np.median(cosines)) >= 0.90,
        "positive_cosine_fraction_at_least_0_95": float(np.mean(cosines > 0)) >= 0.95,
        "median_top3_transition_jaccard_at_least_0_80": (
            float(np.median(transition_jaccard)) >= 0.80
        ),
    }
    median_cosine = float(np.median(cosines))
    verdict = "PASS" if all(checks.values()) else ("MARGINAL" if median_cosine >= 0.75 else "STOP")
    return {
        "eligible_genes": int(eligible_global.size),
        "median_cosine": median_cosine,
        "mean_cosine": float(cosines.mean()),
        "cosine_p10": float(np.quantile(cosines, 0.10)),
        "positive_cosine_fraction": float(np.mean(cosines > 0)),
        "median_top3_transition_jaccard": float(np.median(transition_jaccard)),
        "mean_top3_transition_jaccard": float(transition_jaccard.mean()),
        "both_transition_sets_empty_fraction": both_empty / len(transition_jaccard),
        "checks": checks,
        "verdict": verdict,
    }


def input_record(paths: list[Path], root: Path) -> list[dict[str, object]]:
    records = []
    for path in paths:
        records.append(
            {
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return records


def analyze_dataset(config: dict[str, str], root: Path) -> dict[str, object]:
    started = time.monotonic()
    dataset = config["id"]
    oracle_gene = root / config["oracle_gene"]
    replay_gene = root / config["replay_gene"]
    oracle_velocity = root / config["oracle_velocity"]
    replay_velocity = root / config["replay_velocity"]
    print(f"[{dataset}] loading Gene matrices", file=sys.stderr, flush=True)

    barcodes = read_lines(oracle_gene / "barcodes.tsv")
    if len(set(barcodes)) != len(barcodes):
        raise SystemExit(f"duplicate oracle barcodes in {dataset}")
    gene_ids = require_same_features(
        oracle_gene / "features.tsv", replay_gene / "features.tsv"
    )
    replay_columns = replay_column_indices(replay_gene / "barcodes.tsv", barcodes)
    oracle_raw = load_mtx_cells(oracle_gene / "matrix.mtx")
    replay_raw = load_mtx_cells(replay_gene / "matrix.mtx", replay_columns)
    if oracle_raw.shape != replay_raw.shape or oracle_raw.shape[1] != len(gene_ids):
        raise SystemExit(f"Gene matrix dimensions differ in {dataset}")

    detected = np.asarray(oracle_raw.getnnz(axis=0)).ravel() >= 10
    oracle_raw = oracle_raw[:, detected]
    replay_raw = replay_raw[:, detected]
    retained_gene_ids = np.asarray(gene_ids, dtype=object)[detected]
    oracle_log = normalize_log(oracle_raw)
    replay_log = normalize_log(replay_raw)
    hvg = oracle_hvgs(oracle_log, N_HVG)

    # Establish oracle groups before constructing the preregistered cross-cluster control.
    temporary = oracle_log[:, hvg].toarray()
    temporary_mean = temporary.mean(axis=0)
    temporary_std = temporary.std(axis=0, ddof=1)
    temporary_std[temporary_std == 0] = 1.0
    temporary_pca = PCA(
        n_components=min(N_PCS, temporary.shape[0] - 1, temporary.shape[1]),
        svd_solver="arpack",
        random_state=SEED,
    )
    temporary_pcs = temporary_pca.fit_transform(
        (temporary - temporary_mean) / temporary_std
    )
    temporary_neighbors, _, temporary_adjacency = neighbors(temporary_pcs)
    oracle_clusters_for_control = leiden(temporary_adjacency)
    control_sources, control_selected = choose_control_sources(
        dataset, barcodes, oracle_clusters_for_control
    )
    control_raw = oracle_raw[control_sources].tocsr()
    control_log = normalize_log(control_raw)

    oracle_pcs, replay_pcs, control_pcs, _, gene_std, pca = fit_shared_pca(
        oracle_log, replay_log, control_log, hvg
    )
    oracle_neighbors, _, oracle_adjacency = neighbors(oracle_pcs)
    replay_neighbors, _, replay_adjacency = neighbors(replay_pcs)
    control_neighbors, _, control_adjacency = neighbors(control_pcs)
    oracle_clusters = leiden(oracle_adjacency)
    replay_clusters = leiden(replay_adjacency)
    control_clusters = leiden(control_adjacency)
    if not np.array_equal(oracle_clusters, oracle_clusters_for_control):
        raise SystemExit("oracle clustering changed between control construction and final fit")

    replay_metrics = arm_metrics(
        oracle_clusters,
        replay_clusters,
        oracle_neighbors,
        replay_neighbors,
        oracle_log,
        replay_log,
        oracle_raw,
        replay_raw,
    )
    control_metrics = arm_metrics(
        oracle_clusters,
        control_clusters,
        oracle_neighbors,
        control_neighbors,
        oracle_log,
        control_log,
        oracle_raw,
        control_raw,
    )
    gene_result = {
        "called_cells": len(barcodes),
        "retained_genes": int(detected.sum()),
        "highly_variable_genes": int(hvg.sum()),
        "oracle_clusters": len(set(oracle_clusters)),
        "control_selected_cells": int(control_selected.size),
        "replay": replay_metrics,
        "positive_control": control_metrics,
        "gate": gene_gate(replay_metrics, control_metrics),
    }

    print(f"[{dataset}] loading velocity matrices", file=sys.stderr, flush=True)
    require_same_features(
        oracle_velocity / "features.tsv", replay_velocity / "features.tsv"
    )
    velocity_oracle_barcodes = read_lines(oracle_velocity / "barcodes.tsv")
    if velocity_oracle_barcodes != barcodes:
        raise SystemExit(f"Gene/velocity oracle barcode order differs in {dataset}")
    velocity_columns = replay_column_indices(replay_velocity / "barcodes.tsv", barcodes)
    oracle_s = load_mtx_cells(oracle_velocity / "spliced.mtx")[:, detected]
    oracle_u = load_mtx_cells(oracle_velocity / "unspliced.mtx")[:, detected]
    replay_s = load_mtx_cells(replay_velocity / "spliced.mtx", velocity_columns)[:, detected]
    replay_u = load_mtx_cells(replay_velocity / "unspliced.mtx", velocity_columns)[:, detected]
    velocity_result = velocity_metrics(
        oracle_s,
        oracle_u,
        replay_s,
        replay_u,
        hvg,
        gene_std,
        pca,
        oracle_pcs,
        oracle_neighbors,
    )
    del oracle_raw, replay_raw, control_raw, oracle_log, replay_log, control_log
    del oracle_s, oracle_u, replay_s, replay_u
    gc.collect()

    paths = [
        oracle_gene / "matrix.mtx",
        oracle_gene / "features.tsv",
        oracle_gene / "barcodes.tsv",
        replay_gene / "matrix.mtx",
        replay_gene / "features.tsv",
        replay_gene / "barcodes.tsv",
        oracle_velocity / "spliced.mtx",
        oracle_velocity / "unspliced.mtx",
        oracle_velocity / "features.tsv",
        oracle_velocity / "barcodes.tsv",
        replay_velocity / "spliced.mtx",
        replay_velocity / "unspliced.mtx",
        replay_velocity / "features.tsv",
        replay_velocity / "barcodes.tsv",
    ]
    return {
        "id": dataset,
        "role": config["role"],
        "organism": config["organism"],
        "chemistry": config["chemistry"],
        "gene": gene_result,
        "velocity": velocity_result,
        "retained_gene_id_sha256": hashlib.sha256(
            ("\n".join(map(str, retained_gene_ids)) + "\n").encode()
        ).hexdigest(),
        "inputs": input_record(paths, root),
        "wall_seconds": time.monotonic() - started,
    }


def package_versions() -> dict[str, str]:
    names = [
        "scanpy",
        "anndata",
        "numpy",
        "scipy",
        "scikit-learn",
        "igraph",
        "leidenalg",
        "pyyaml",
    ]
    return {name: importlib.metadata.version(name) for name in names}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--dataset", action="append", help="run only these dataset ids")
    args = parser.parse_args()

    manifest = yaml.safe_load(args.manifest.read_text())
    selected = set(args.dataset or [])
    configs = [
        item for item in manifest["datasets"] if not selected or item["id"] in selected
    ]
    if selected != {item["id"] for item in configs} and selected:
        raise SystemExit(f"unknown datasets: {sorted(selected - {item['id'] for item in configs})}")
    start = time.monotonic()
    result = {
        "date": "2026-08-31",
        "preregistered_spec_commit": manifest["preregistered_spec_commit"],
        "claim_scope": (
            "downstream stability against fresh STARsolo; no cell-type truth and no claim of "
            "exact fresh annotation-aware processing"
        ),
        "software": {
            "python": platform.python_version(),
            "packages": package_versions(),
            "seed": SEED,
            "hostname": platform.node(),
            "platform": platform.platform(),
            "threads": int(os.environ.get("OMP_NUM_THREADS", "1")),
        },
        "parameters": {
            "target_sum": TARGET_SUM,
            "hvg_flavor": "seurat",
            "hvg_count": N_HVG,
            "principal_components": N_PCS,
            "neighbors": N_NEIGHBORS,
            "leiden_resolution": LEIDEN_RESOLUTION,
            "control_cell_fraction": 0.20,
        },
        "datasets": [],
    }
    for config in configs:
        result["datasets"].append(analyze_dataset(config, args.project_root))
        args.out.parent.mkdir(parents=True, exist_ok=True)
        partial = dict(result)
        partial["status"] = "partial"
        args.out.write_text(json.dumps(partial, indent=2) + "\n")

    confirmatory = [item for item in result["datasets"] if item["role"] == "confirmatory"]
    stress = [item for item in result["datasets"] if item["role"] != "confirmatory"]
    result["summary"] = {
        "human_gene_verdict": (
            "PASS" if confirmatory and all(
                item["gene"]["gate"]["verdict"] == "PASS" for item in confirmatory
            ) else "FAIL"
        ),
        "human_velocity_verdict": (
            "PASS" if confirmatory and all(
                item["velocity"]["verdict"] == "PASS" for item in confirmatory
            ) else (
                "STOP" if any(item["velocity"]["verdict"] == "STOP" for item in confirmatory)
                else "MARGINAL"
            )
        ),
        "stress_test_gene_verdicts": {
            item["id"]: item["gene"]["gate"]["verdict"] for item in stress
        },
        "stress_test_velocity_verdicts": {
            item["id"]: item["velocity"]["verdict"] for item in stress
        },
        "wall_seconds": time.monotonic() - start,
    }
    result["status"] = "complete"
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
