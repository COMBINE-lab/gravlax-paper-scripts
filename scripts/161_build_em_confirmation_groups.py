#!/usr/bin/env python3
"""Select and write leakage-controlled groups for untouched EM confirmation."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import statistics
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


def load_helpers() -> object:
    path = Path(__file__).with_name("126_downstream_fidelity.py")
    spec = importlib.util.spec_from_file_location("downstream_gate", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


gate = load_helpers()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_candidates(candidate_dir: Path, dataset: str, seeds: list[int]) -> set[str]:
    genes: set[str] = set()
    for seed in seeds:
        path = candidate_dir / f"{dataset}.seed-{seed}.genes.txt"
        if not path.is_file():
            raise SystemExit(f"missing candidate file: {path}")
        genes.update(line.strip() for line in path.open() if line.strip())
    return genes


def replay_column_indices(path: Path, barcodes: list[str]) -> list[int]:
    """Find a small called-cell set in a potentially multi-million-barcode replay axis."""
    positions = {barcode: offset for offset, barcode in enumerate(barcodes)}
    if len(positions) != len(barcodes):
        raise SystemExit("called-cell barcode list contains duplicates")
    columns = [-1] * len(barcodes)
    with path.open() as handle:
        for replay_index, line in enumerate(handle):
            barcode = line.rstrip("\n")
            called_index = positions.get(barcode)
            if called_index is not None:
                columns[called_index] = replay_index
    missing = [barcodes[index] for index, column in enumerate(columns) if column < 0]
    if missing:
        raise SystemExit(
            f"{len(missing)} called-cell barcodes are absent from replay axis; first={missing[0]}"
        )
    return columns


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


def partition_summary(labels: list[np.ndarray], seeds: list[int]) -> dict[str, object]:
    pairwise: list[float] = []
    centrality = np.zeros(len(labels), dtype=np.float64)
    for left, right in combinations(range(len(labels)), 2):
        value = float(adjusted_rand_score(labels[left], labels[right]))
        pairwise.append(value)
        centrality[left] += value
        centrality[right] += value
    centrality /= len(labels) - 1
    medoid_index = int(centrality.argmax())
    selected = labels[medoid_index]
    counts = np.bincount(selected)
    return {
        "median_pairwise_seed_ari": statistics.median(pairwise),
        "minimum_pairwise_seed_ari": min(pairwise),
        "mean_pairwise_seed_ari": statistics.fmean(pairwise),
        "medoid_seed": seeds[medoid_index],
        "medoid_mean_ari": float(centrality[medoid_index]),
        "medoid_labels": selected,
        "cluster_count": int(len(counts)),
        "minimum_group_cells": int(counts.min()),
        "group_sizes": counts.tolist(),
    }


def shuffled_labels(dataset: str, barcodes: list[str], labels: np.ndarray) -> np.ndarray:
    order = sorted(
        range(len(barcodes)),
        key=lambda index: hashlib.sha256(
            f"em-confirmation-shuffle-v1\0{dataset}\0{barcodes[index]}".encode()
        ).digest(),
    )
    shuffled = np.empty_like(labels)
    shuffled[np.asarray(order, dtype=np.int64)] = np.sort(labels)
    return shuffled


def write_groups(path: Path, barcodes: list[str], labels: np.ndarray) -> None:
    with path.open("w") as handle:
        for barcode, label in zip(barcodes, labels, strict=True):
            handle.write(f"{barcode}\tg{int(label)}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", choices=("D4", "D3"), required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--replay-dir", type=Path, required=True)
    parser.add_argument("--called-barcodes", type=Path, required=True)
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest = yaml.safe_load(args.manifest.read_text())
    grouping = manifest["groups"]
    seeds = [int(value) for value in manifest["targets"]["seeds"]]
    leiden_seeds = [int(value) for value in grouping["leiden_seeds"]]
    barcodes = gate.read_lines(args.called_barcodes)
    columns = replay_column_indices(args.replay_dir / "barcodes.tsv", barcodes)
    raw = gate.load_mtx_cells(args.replay_dir / "matrix.mtx", columns)
    features = gate.feature_ids(args.replay_dir / "features.tsv")
    if raw.shape != (len(barcodes), len(features)):
        raise SystemExit(f"matrix dimensions {raw.shape} do not match barcodes/features")

    excluded = read_candidates(args.candidate_dir, args.dataset, seeds)
    detected = np.asarray(raw.getnnz(axis=0)).ravel() >= int(
        grouping["detected_gene_minimum_cells"]
    )
    not_candidate = np.asarray([gene not in excluded for gene in features], dtype=bool)
    eligible_genes = detected & not_candidate
    if int(eligible_genes.sum()) < int(grouping["minimum_eligible_genes"]):
        raise SystemExit("too few eligible genes remain after candidate exclusion")

    log = gate.normalize_log(raw[:, eligible_genes])
    hvg = gate.oracle_hvgs(log, int(grouping["hvg_count_max"]))
    dense = log[:, hvg].toarray()
    mean = dense.mean(axis=0)
    std = dense.std(axis=0, ddof=1)
    std[std == 0] = 1.0
    maximum_dimension = max(int(value) for value in grouping["dimensions"])
    components = min(maximum_dimension, dense.shape[0] - 1, dense.shape[1] - 1)
    if components < maximum_dimension:
        raise SystemExit(
            f"insufficient PCA rank {components} for registered dimension {maximum_dimension}"
        )
    model = PCA(n_components=components, svd_solver="arpack", random_state=1729)
    pcs = model.fit_transform((dense - mean) / std)

    cluster_min, cluster_max = (int(value) for value in grouping["eligible_cluster_count"])
    min_group_cells = int(grouping["minimum_group_cells"])
    min_ari = float(grouping["minimum_median_pairwise_seed_ARI"])
    records: list[dict[str, object]] = []
    stored_labels: dict[tuple[int, float], np.ndarray] = {}
    for dimension_value in grouping["dimensions"]:
        dimension = int(dimension_value)
        _, _, graph = gate.neighbors(pcs[:, :dimension])
        for resolution_value in grouping["resolutions"]:
            resolution = float(resolution_value)
            labels = [leiden(graph, resolution, seed) for seed in leiden_seeds]
            summary = partition_summary(labels, leiden_seeds)
            selected = summary.pop("medoid_labels")
            eligible = (
                cluster_min <= int(summary["cluster_count"]) <= cluster_max
                and int(summary["minimum_group_cells"]) >= min_group_cells
                and float(summary["median_pairwise_seed_ari"]) >= min_ari
            )
            records.append(
                {
                    "dimension": dimension,
                    "resolution": resolution,
                    **summary,
                    "eligible": eligible,
                }
            )
            stored_labels[(dimension, resolution)] = selected

    eligible_records = [record for record in records if bool(record["eligible"])]
    if not eligible_records:
        args.out_dir.mkdir(parents=True, exist_ok=False)
        result = {
            "status": "STOP_no_eligible_group_configuration",
            "dataset": args.dataset,
            "configurations": records,
        }
        (args.out_dir / "groups.json").write_text(json.dumps(result, indent=2) + "\n")
        raise SystemExit("no registered grouping configuration is eligible")
    selected_record = min(
        eligible_records,
        key=lambda record: (
            -float(record["median_pairwise_seed_ari"]),
            int(record["dimension"]),
            abs(float(record["resolution"]) - 0.50),
            float(record["resolution"]),
        ),
    )
    key = (int(selected_record["dimension"]), float(selected_record["resolution"]))
    selected_labels = stored_labels[key]
    shuffled = shuffled_labels(args.dataset, barcodes, selected_labels)
    args.out_dir.mkdir(parents=True, exist_ok=False)
    real_path = args.out_dir / f"{args.dataset}.real.tsv"
    shuffled_path = args.out_dir / f"{args.dataset}.shuffled.tsv"
    write_groups(real_path, barcodes, selected_labels)
    write_groups(shuffled_path, barcodes, shuffled)

    result = {
        "status": "complete",
        "analysis": "untouched_candidate_excluded_replay_only_group_selection",
        "dataset": args.dataset,
        "manifest": str(args.manifest),
        "manifest_sha256": sha256(args.manifest),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scanpy": sc.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "called_cells": len(barcodes),
        "candidate_genes_excluded": len(excluded),
        "detected_genes_before_exclusion": int(detected.sum()),
        "eligible_genes_after_exclusion": int(eligible_genes.sum()),
        "hvg_count": int(hvg.sum()),
        "selection_rule": grouping["selection"],
        "tie_break": grouping["tie_break"],
        "selected": selected_record,
        "eligible_configuration_count": len(eligible_records),
        "configuration_count": len(records),
        "configurations": records,
        "real_groups": {"path": str(real_path), "sha256": sha256(real_path)},
        "shuffled_groups": {
            "path": str(shuffled_path),
            "sha256": sha256(shuffled_path),
            "sizes_preserved": bool(
                np.array_equal(
                    np.sort(np.bincount(shuffled)), np.sort(np.bincount(selected_labels))
                )
            ),
        },
        "called_barcodes_sha256": sha256(args.called_barcodes),
        "replay_barcodes_sha256": sha256(args.replay_dir / "barcodes.tsv"),
        "replay_matrix_sha256": sha256(args.replay_dir / "matrix.mtx"),
        "replay_features_sha256": sha256(args.replay_dir / "features.tsv"),
        "candidate_files": {
            f"seed_{seed}": sha256(args.candidate_dir / f"{args.dataset}.seed-{seed}.genes.txt")
            for seed in seeds
        },
    }
    (args.out_dir / "groups.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
