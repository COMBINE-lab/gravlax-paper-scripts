#!/usr/bin/env python3
"""Build the pre-locked D5 broad-cell group map without using FNBP1 evidence."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def decode(values: np.ndarray) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def zscore(values: np.ndarray) -> np.ndarray:
    std = float(values.std())
    if std == 0.0:
        return np.zeros_like(values, dtype=np.float64)
    return (values - float(values.mean())) / std


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    manifest = json.loads(args.manifest.read_text())
    expected = manifest["dataset"]["filtered_feature_bc_matrix_h5_sha256"]
    observed = sha256(args.h5)
    if observed != expected:
        raise SystemExit(f"H5 digest mismatch: expected {expected}, observed {observed}")
    if manifest["status"] != "LABEL_RULE_LOCKED_BEFORE_ANY_CELL_TYPE_EVENT_QUERY":
        raise SystemExit("manifest does not carry the pre-query lock status")

    modules: dict[str, dict[str, list[str]]] = manifest["marker_modules"]
    marker_genes = sorted({gene for lineage in modules.values() for genes in lineage.values() for gene in genes})
    with h5py.File(args.h5) as handle:
        matrix = handle["matrix"]
        barcodes = decode(matrix["barcodes"][:])
        names = decode(matrix["features"]["name"][:])
        feature_types = decode(matrix["features"]["feature_type"][:])
        shape = tuple(int(value) for value in matrix["shape"][:])
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]

    if shape != (len(names), len(barcodes)):
        raise SystemExit(f"unexpected matrix shape {shape} for {len(names)} features x {len(barcodes)} cells")
    if len(barcodes) != manifest["dataset"]["called_nuclei"]:
        raise SystemExit(f"expected {manifest['dataset']['called_nuclei']} called nuclei, found {len(barcodes)}")
    if len(set(barcodes)) != len(barcodes):
        raise SystemExit("duplicate H5 barcodes")

    rows_by_name: dict[str, list[int]] = {}
    for row, (name, feature_type) in enumerate(zip(names, feature_types)):
        if feature_type == "Gene Expression":
            rows_by_name.setdefault(name, []).append(row)
    missing = sorted(set(marker_genes) - rows_by_name.keys())
    if missing:
        raise SystemExit(f"missing locked marker genes: {missing}")

    marker_slot = {gene: offset for offset, gene in enumerate(marker_genes)}
    row_to_slot = np.full(len(names), -1, dtype=np.int16)
    for gene, slot in marker_slot.items():
        for row in rows_by_name[gene]:
            row_to_slot[row] = slot
    is_gex = np.asarray([feature_type == "Gene Expression" for feature_type in feature_types])
    is_fnbp1 = np.asarray([name == "FNBP1" and feature_type == "Gene Expression" for name, feature_type in zip(names, feature_types)])
    counts = np.zeros((len(marker_genes), len(barcodes)), dtype=np.float64)
    library = np.zeros(len(barcodes), dtype=np.float64)
    for cell in range(len(barcodes)):
        start, end = int(indptr[cell]), int(indptr[cell + 1])
        cell_rows = indices[start:end]
        cell_values = data[start:end].astype(np.float64, copy=False)
        keep_gex = is_gex[cell_rows]
        keep_denominator = keep_gex & ~is_fnbp1[cell_rows]
        library[cell] = float(cell_values[keep_denominator].sum())
        slots = row_to_slot[cell_rows]
        keep_markers = slots >= 0
        np.add.at(counts[:, cell], slots[keep_markers], cell_values[keep_markers])
    if np.any(library <= 0):
        raise SystemExit(f"{int((library <= 0).sum())} called nuclei have no non-FNBP1 GEX UMIs")

    transformed = np.log1p(counts * (10000.0 / library)[None, :])
    marker_z = np.vstack([np.clip(zscore(row), -3.0, 3.0) for row in transformed])
    lineage_scores: dict[str, np.ndarray] = {}
    lineage_detected: dict[str, np.ndarray] = {}
    selected_submodule: dict[str, np.ndarray] = {}
    for lineage, submodules in modules.items():
        sub_names = list(submodules)
        raw_scores = []
        detections = []
        for sub_name in sub_names:
            slots = [marker_slot[gene] for gene in submodules[sub_name]]
            raw_scores.append(zscore(marker_z[slots, :].mean(axis=0)))
            detections.append((counts[slots, :] > 0).sum(axis=0))
        score_matrix = np.vstack(raw_scores)
        detection_matrix = np.vstack(detections)
        selected = np.argmax(score_matrix, axis=0)
        columns = np.arange(len(barcodes))
        lineage_scores[lineage] = score_matrix[selected, columns]
        lineage_detected[lineage] = detection_matrix[selected, columns]
        selected_submodule[lineage] = np.asarray([sub_names[index] for index in selected])

    lineage_names = list(modules)
    score_matrix = np.vstack([lineage_scores[name] for name in lineage_names])
    detected_matrix = np.vstack([lineage_detected[name] for name in lineage_names])
    order = np.argsort(score_matrix, axis=0)
    top_index, runner_index = order[-1, :], order[-2, :]
    columns = np.arange(len(barcodes))
    top_score = score_matrix[top_index, columns]
    runner_score = score_matrix[runner_index, columns]
    top_detected = detected_matrix[top_index, columns]
    runner_detected = detected_matrix[runner_index, columns]
    margin = top_score - runner_score

    thresholds = manifest["label_rule"]["thresholds"]
    minimum_detected = int(thresholds["minimum_detected_markers_in_selected_submodule"])
    minimum_score = float(thresholds["minimum_lineage_score"])
    minimum_margin = float(thresholds["minimum_top_vs_runner_margin"])
    high_umi = float(np.quantile(library, float(thresholds["high_GEX_UMI_quantile"])))
    labels: list[str] = []
    reasons: list[str] = []
    for cell in range(len(barcodes)):
        top = lineage_names[int(top_index[cell])]
        top_eligible = top_detected[cell] >= minimum_detected
        runner_eligible = runner_detected[cell] >= minimum_detected
        conflict = (
            top_eligible
            and runner_eligible
            and top_score[cell] >= minimum_score
            and runner_score[cell] >= minimum_score
            and margin[cell] < minimum_margin
        )
        high_count_conflict = (
            library[cell] > high_umi
            and top_eligible
            and runner_eligible
            and runner_score[cell] >= 0.0
        )
        if conflict or high_count_conflict:
            labels.append("ambiguous_or_doublet")
            reasons.append("score_conflict" if conflict else "high_count_multilineage")
        elif top_eligible and top_score[cell] >= minimum_score and margin[cell] >= minimum_margin:
            labels.append(top)
            reasons.append("assigned")
        else:
            labels.append("unresolved")
            reasons.append("insufficient_marker_score_or_margin")

    archive_barcodes = [barcode.removesuffix("-1") for barcode in barcodes]
    if len(set(archive_barcodes)) != len(archive_barcodes):
        raise SystemExit("removing the 10x GEM-group suffix creates duplicate archive barcodes")
    group_path = args.out / "broad-cell-groups.tsv"
    detail_path = args.out / "broad-cell-label-details.tsv"
    with group_path.open("w") as handle:
        for barcode, label in zip(archive_barcodes, labels):
            handle.write(f"{barcode}\t{label}\n")
    with detail_path.open("w") as handle:
        handle.write("barcode\tgroup\treason\ttop_lineage\trunner_lineage\ttop_score\trunner_score\tmargin\tGEX_UMIs_excluding_FNBP1\n")
        for cell, barcode in enumerate(barcodes):
            handle.write(
                f"{barcode}\t{labels[cell]}\t{reasons[cell]}\t{lineage_names[int(top_index[cell])]}\t"
                f"{lineage_names[int(runner_index[cell])]}\t{top_score[cell]:.8g}\t"
                f"{runner_score[cell]:.8g}\t{margin[cell]:.8g}\t{int(library[cell])}\n"
            )

    group_counts = {label: labels.count(label) for label in sorted(set(labels))}
    summary = {
        "schema": "gravlax.fnbp1-broad-cell-groups.v1",
        "status": "LABELS_FROZEN_BEFORE_EVENT_QUERY",
        "manifest": str(args.manifest),
        "manifest_sha256": sha256(args.manifest),
        "input": {"path": str(args.h5), "sha256": observed, "called_nuclei": len(barcodes)},
        "groups": group_counts,
        "high_GEX_UMI_threshold_excluding_FNBP1": high_umi,
        "marker_genes_present": marker_genes,
        "group_map": group_path.name,
        "group_map_sha256": sha256(group_path),
        "detail_table": detail_path.name,
        "detail_table_sha256": sha256(detail_path),
        "software": {"h5py": h5py.__version__, "numpy": np.__version__},
    }
    summary_path = args.out / "groups-summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
