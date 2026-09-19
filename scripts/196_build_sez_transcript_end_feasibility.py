#!/usr/bin/env python3
"""Build the event-blind GSE234790 transcript-end design and frozen cell maps.

This program reads only published cell metadata.  It never reads an alignment,
archive, endpoint, poly(A) catalogue, or gene-level molecular count.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED_METADATA_SHA256 = "5a26c7e02c02a06ded7eec9e6c24fa62c6a9718ab37a3c212b472ed331a6c593"
BARCODE = re.compile(r"([ACGT]{16})-1$")
DONORS = tuple("ABCDEFGH")
YOUNG = frozenset("ABCD")

# Published cluster identities from Puvogel et al. (eNeuro 2024).  Only these
# immutable published labels define the biological groups; no endpoint-derived
# or Gravlax-derived clustering is used.
CLUSTER_GROUP = {
    "8": "ependymal",
    "0": "astro_nsc",
    "7": "astro_nsc",
    "16": "astro_nsc",
    "9": "neuroblast",
    "11": "immature_neuron",
    "13": "immature_neuron",
    "14": "immature_neuron",
    "2": "mature_neuron",
    "3": "mature_neuron",
    "10": "mature_neuron",
    "12": "mature_neuron",
    "19": "mature_neuron",
    "4": "opc",
    "1": "oligodendrocyte",
    "5": "oligodendrocyte",
    "6": "microglia",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    observed_sha256 = sha256(args.metadata)
    if observed_sha256 != EXPECTED_METADATA_SHA256:
        raise SystemExit(
            f"metadata digest mismatch: {observed_sha256} != {EXPECTED_METADATA_SHA256}"
        )

    rows: dict[str, list[tuple[str, str]]] = defaultdict(list)
    cluster_counts: Counter[tuple[str, str]] = Counter()
    group_counts: Counter[tuple[str, str]] = Counter()
    ages: dict[str, float] = {}
    seen: set[tuple[str, str]] = set()
    other_clusters: Counter[str] = Counter()
    with gzip.open(args.metadata, "rt", newline="") as handle:
        reader = csv.DictReader(handle)
        if "sample" not in (reader.fieldnames or ()) or "integrated_snn_res.0.4" not in (
            reader.fieldnames or ()
        ):
            raise SystemExit("metadata lacks required sample or published-cluster columns")
        cell_column = (reader.fieldnames or [None])[0]
        if cell_column is None:
            raise SystemExit("metadata lacks a cell identifier column")
        for row in reader:
            donor = row["sample"].strip()
            if donor not in DONORS:
                raise SystemExit(f"unexpected donor {donor!r}")
            match = BARCODE.search(row[cell_column])
            if match is None:
                raise SystemExit(f"cannot extract 10x barcode from {row[cell_column]!r}")
            barcode = match.group(1)
            key = (donor, barcode)
            if key in seen:
                raise SystemExit(f"duplicate donor/barcode {donor}/{barcode}")
            seen.add(key)
            cluster = row["integrated_snn_res.0.4"].strip()
            cluster_counts[(donor, cluster)] += 1
            age = float(row["age"])
            if donor in ages and ages[donor] != age:
                raise SystemExit(f"donor {donor} has inconsistent ages")
            ages[donor] = age
            group = CLUSTER_GROUP.get(cluster, "other")
            if group == "other":
                other_clusters[cluster] += 1
            rows[donor].append((barcode, group))
            group_counts[(donor, group)] += 1

    if set(rows) != set(DONORS):
        raise SystemExit(f"metadata donors {sorted(rows)} do not match {list(DONORS)}")
    args.out_root.mkdir(parents=True, exist_ok=True)
    group_paths: dict[str, Path] = {}
    for donor in DONORS:
        path = args.out_root / f"donor-{donor}.tsv"
        with path.open("w", newline="") as handle:
            for barcode, group in sorted(rows[donor]):
                handle.write(f"{barcode}\t{group}\n")
        group_paths[donor] = path

    args.design.parent.mkdir(parents=True, exist_ok=True)
    with args.design.open("w", newline="") as handle:
        handle.write("sample\tcondition\tarchive\tgroups\n")
        for donor in DONORS:
            condition = "young" if donor in YOUNG else "middle"
            archive = (
                f"../../../runs/post-v1/sez-transcript-end-atlas-r1/"
                f"donor-{donor}/sez.called.aie"
            )
            groups = (
                f"../../results/post-v1-sez-transcript-end-feasibility/"
                f"groups/donor-{donor}.tsv"
            )
            handle.write(f"{donor}\t{condition}\t{archive}\t{groups}\n")

    tracked_groups = sorted(set(CLUSTER_GROUP.values()) | {"other"})
    summary = {
        "schema": "gravlax.sez-transcript-end-feasibility.v1",
        "status": "PASS_EVENT_BLIND",
        "endpoint_evidence_inspected": False,
        "inputs": {
            "metadata": {
                "logical_path": "data/multidonor-fnbp1/gse234790/GSE234790_meta_data_counts.csv.gz",
                "sha256": observed_sha256,
                "bytes": args.metadata.stat().st_size,
            }
        },
        "published_cluster_mapping": dict(sorted(CLUSTER_GROUP.items(), key=lambda item: int(item[0]))),
        "primary_groups": ["astro_nsc", "mature_neuron"],
        "donors": [],
        "other_published_clusters": dict(sorted(other_clusters.items(), key=lambda item: int(item[0]))),
    }
    for donor in DONORS:
        counts = {group: group_counts[(donor, group)] for group in tracked_groups}
        summary["donors"].append(
            {
                "sample": donor,
                "condition": "young" if donor in YOUNG else "middle",
                "age": ages[donor],
                "tracked_nuclei": len(rows[donor]),
                "groups": counts,
                "primary_minimum_nuclei": min(counts["astro_nsc"], counts["mature_neuron"]),
                "groups_sha256": sha256(group_paths[donor]),
            }
        )
    if any(row["primary_minimum_nuclei"] < 150 for row in summary["donors"]):
        summary["status"] = "STOP_EVENT_BLIND_PRIMARY_GROUP_DEPTH"
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
