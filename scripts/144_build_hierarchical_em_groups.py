#!/usr/bin/env python3
"""Build leakage-controlled replay-only macro groups for hierarchical EM."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import warnings
from itertools import combinations
from pathlib import Path

import anndata
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import sklearn
import yaml
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score


def load_downstream_helpers() -> object:
    path = Path(__file__).with_name("126_downstream_fidelity.py")
    spec = importlib.util.spec_from_file_location("downstream_gate", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


gate = load_downstream_helpers()


def leiden(adjacency: sp.csr_matrix, resolution: float, seed: int) -> np.ndarray:
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


def medoid(labels: list[np.ndarray], seeds: list[int]) -> tuple[np.ndarray, int, float]:
    scores = np.zeros(len(labels), dtype=np.float64)
    for left, right in combinations(range(len(labels)), 2):
        value = adjusted_rand_score(labels[left], labels[right])
        scores[left] += value
        scores[right] += value
    scores /= len(labels) - 1
    index = int(scores.argmax())
    return labels[index], seeds[index], float(scores[index])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_groups(path: Path, barcodes: list[str], labels: np.ndarray) -> None:
    with path.open("w") as handle:
        for barcode, label in zip(barcodes, labels, strict=True):
            handle.write(f"{barcode}\tg{int(label)}\n")


def shuffled_labels(dataset: str, barcodes: list[str], labels: np.ndarray) -> np.ndarray:
    order = sorted(
        range(len(barcodes)),
        key=lambda i: hashlib.sha256(
            f"hierarchical-em-shuffle-v1\0{dataset}\0{barcodes[i]}".encode()
        ).digest(),
    )
    source = np.sort(labels)
    shuffled = np.empty_like(labels)
    shuffled[np.asarray(order, dtype=np.int64)] = source
    return shuffled


def candidate_union(candidate_dir: Path, dataset: str, seeds: list[int]) -> set[str]:
    genes: set[str] = set()
    for seed in seeds:
        path = candidate_dir / f"{dataset}.seed-{seed}.genes.txt"
        if not path.is_file():
            raise SystemExit(f"missing candidate-gene input: {path}")
        genes.update(line.strip() for line in path.open() if line.strip())
    return genes


def build_dataset(
    config: dict[str, object],
    manifest: dict[str, object],
    project_root: Path,
    candidate_dir: Path,
    out_dir: Path,
) -> dict[str, object]:
    dataset = str(config["id"])
    replay_dir = project_root / str(config["replay_gene"])
    called_path = project_root / str(config["oracle_called_barcodes"])
    barcodes = gate.read_lines(called_path)
    columns = gate.replay_column_indices(replay_dir / "barcodes.tsv", barcodes)
    raw = gate.load_mtx_cells(replay_dir / "matrix.mtx", columns)
    features = gate.feature_ids(replay_dir / "features.tsv")
    if raw.shape[1] != len(features):
        raise SystemExit(f"{dataset}: matrix/feature dimension mismatch")

    seeds = [int(value) for value in manifest["targets"]["seeds"]]
    excluded = candidate_union(candidate_dir, dataset, seeds)
    detected = np.asarray(raw.getnnz(axis=0)).ravel() >= 10
    not_candidate = np.asarray([gene not in excluded for gene in features], dtype=bool)
    eligible = detected & not_candidate
    minimum = int(manifest["grouping"]["minimum_eligible_genes"])
    if int(eligible.sum()) < minimum:
        raise SystemExit(
            f"{dataset}: only {int(eligible.sum())} eligible genes remain after candidate exclusion"
        )

    raw = raw[:, eligible]
    log = gate.normalize_log(raw)
    hvg = gate.oracle_hvgs(log, int(manifest["grouping"]["hvg_count_max"]))
    dense = log[:, hvg].toarray()
    mean = dense.mean(axis=0)
    std = dense.std(axis=0, ddof=1)
    std[std == 0] = 1.0
    components = min(40, dense.shape[0] - 1, dense.shape[1] - 1)
    if components < int(config["dimension"]):
        raise SystemExit(f"{dataset}: insufficient PCA rank ({components})")
    model = PCA(n_components=components, svd_solver="arpack", random_state=1729)
    pcs = model.fit_transform((dense - mean) / std)
    _, _, graph = gate.neighbors(pcs[:, : int(config["dimension"])])
    leiden_seeds = [int(value) for value in manifest["grouping"]["leiden_seeds"]]
    labels = [leiden(graph, float(config["resolution"]), seed) for seed in leiden_seeds]
    selected, selected_seed, centrality = medoid(labels, leiden_seeds)
    counts = np.bincount(selected)
    if int(counts.min()) < 20:
        raise SystemExit(f"{dataset}: smallest selected group has only {int(counts.min())} cells")

    shuffled = shuffled_labels(dataset, barcodes, selected)
    real_path = out_dir / f"{dataset}.real.tsv"
    shuffled_path = out_dir / f"{dataset}.shuffled.tsv"
    write_groups(real_path, barcodes, selected)
    write_groups(shuffled_path, barcodes, shuffled)
    return {
        "id": dataset,
        "called_cells": len(barcodes),
        "candidate_genes_excluded": len(excluded),
        "detected_genes_before_exclusion": int(detected.sum()),
        "eligible_genes_after_exclusion": int(eligible.sum()),
        "hvg_count": int(hvg.sum()),
        "dimension": int(config["dimension"]),
        "resolution": float(config["resolution"]),
        "leiden_seeds": leiden_seeds,
        "medoid_seed": selected_seed,
        "medoid_mean_within_arm_ari": centrality,
        "group_sizes": counts.tolist(),
        "real_groups": {"path": str(real_path), "sha256": sha256(real_path)},
        "shuffled_groups": {
            "path": str(shuffled_path),
            "sha256": sha256(shuffled_path),
            "sizes_preserved": bool(np.array_equal(np.sort(np.bincount(shuffled)), np.sort(counts))),
        },
        "called_barcodes_sha256": sha256(called_path),
        "replay_matrix_sha256": sha256(replay_dir / "matrix.mtx"),
        "replay_features_sha256": sha256(replay_dir / "features.tsv"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest = yaml.safe_load(args.manifest.read_text())
    args.out_dir.mkdir(parents=True, exist_ok=False)
    datasets = [
        build_dataset(config, manifest, args.project_root, args.candidate_dir, args.out_dir)
        for config in manifest["datasets"]
    ]
    result = {
        "status": "complete",
        "analysis": "replay_only_candidate_excluded_hierarchical_em_groups",
        "manifest": str(args.manifest),
        "manifest_sha256": sha256(args.manifest),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scanpy": sc.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "datasets": datasets,
    }
    (args.out_dir / "groups.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
