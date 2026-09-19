#!/usr/bin/env python3
"""Build called-cell and frozen marker-rule group maps for the federation screen.

Annotation-aware count matrices are used only to define cell scopes.  Event discovery and
reduction use the annotation-free archives.  The T/monocyte rule is the rule used by the earlier
APA analysis: T = CD3E > 0 and LYZ == 0; M = LYZ > 3 and CD3E == 0.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.io import mmread


DATASETS = {
    "D0": "runs/oracle/full/v49/Solo.out/Gene/filtered",
    "D1": "runs/oracle/d1/v49/Solo.out/Gene/filtered",
    "D2": "runs/oracle/d2/v49/Solo.out/Gene/filtered",
    "D2p": "runs/oracle/d2p/v49/Solo.out/Gene/filtered",
    "D3": "runs/post-v1/em-promotion-d3-prep-r1/starsolo/Solo.out/Gene/filtered",
    "D4": "runs/post-v1/em-confirmation-d4-prep-r1/starsolo/Solo.out/Gene/filtered",
}
TM_DATASETS = {"D0", "D1", "D3", "D4"}


def read_features(path: Path) -> list[str]:
    return [line.rstrip("\n").split("\t")[1] for line in path.open()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    summary: dict[str, object] = {
        "schema": "gravlax.federated-biology.groups.v1",
        "rule": {"T": "CD3E > 0 and LYZ == 0", "M": "LYZ > 3 and CD3E == 0"},
        "datasets": {},
    }
    for dataset, relative in DATASETS.items():
        source = args.root / relative
        barcodes = [line.strip() for line in (source / "barcodes.tsv").open()]
        called_path = args.out / f"{dataset}.called.tsv"
        with called_path.open("w") as handle:
            for barcode in barcodes:
                handle.write(f"{barcode}\tcalled\n")
        record: dict[str, object] = {
            "source": relative,
            "called_cells": len(barcodes),
            "called_map": called_path.name,
        }
        if dataset in TM_DATASETS:
            features = read_features(source / "features.tsv")
            gene_index = {gene: index for index, gene in enumerate(features)}
            missing = {"CD3E", "LYZ"} - gene_index.keys()
            if missing:
                raise SystemExit(f"{dataset} lacks marker genes: {sorted(missing)}")
            matrix = mmread(source / "matrix.mtx").tocsr()
            if matrix.shape != (len(features), len(barcodes)):
                raise SystemExit(
                    f"{dataset} matrix shape {matrix.shape} does not match "
                    f"{len(features)} features x {len(barcodes)} barcodes"
                )
            cd3e = np.asarray(matrix[gene_index["CD3E"], :].todense()).ravel()
            lyz = np.asarray(matrix[gene_index["LYZ"], :].todense()).ravel()
            groups = [
                (barcode, "T")
                for barcode, cd3e_count, lyz_count in zip(barcodes, cd3e, lyz)
                if cd3e_count > 0 and lyz_count == 0
            ]
            groups.extend(
                (barcode, "M")
                for barcode, cd3e_count, lyz_count in zip(barcodes, cd3e, lyz)
                if lyz_count > 3 and cd3e_count == 0
            )
            tm_path = args.out / f"{dataset}.TM.tsv"
            with tm_path.open("w") as handle:
                for barcode, group in groups:
                    handle.write(f"{barcode}\t{group}\n")
            record.update(
                {
                    "T_cells": sum(group == "T" for _, group in groups),
                    "M_cells": sum(group == "M" for _, group in groups),
                    "TM_map": tm_path.name,
                }
            )
        summary["datasets"][dataset] = record

    (args.out / "groups.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
