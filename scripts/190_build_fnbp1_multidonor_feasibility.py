#!/usr/bin/env python3
"""Event-blind feasibility and frozen broad-cell maps for GSE234790.

The FNBP1 row is recognized only so it can be skipped wholesale.  Its values are
never tokenized, counted, or retained.  No sequence or splice-event input is read.
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

import numpy as np


BARCODE_RE = re.compile(r"([ACGT]{16})-1$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def zscore(values: np.ndarray) -> np.ndarray:
    std = float(values.std())
    if std == 0.0:
        return np.zeros_like(values, dtype=np.float64)
    return (values - float(values.mean())) / std


def verify_input(path: Path, spec: dict[str, object]) -> None:
    observed_size = path.stat().st_size
    observed_sha = sha256(path)
    if observed_size != spec.get("bytes", observed_size):
        raise SystemExit(f"size mismatch for {path}: {observed_size} != {spec['bytes']}")
    if observed_sha != spec["sha256"]:
        raise SystemExit(f"digest mismatch for {path}: {observed_sha} != {spec['sha256']}")


def read_metadata(path: Path) -> tuple[list[str], dict[str, dict[str, object]]]:
    records: dict[str, dict[str, object]] = {}
    order: list[str] = []
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle)
        first = reader.fieldnames[0] if reader.fieldnames else None
        if first is None:
            raise SystemExit("metadata has no row-name column")
        for row in reader:
            cell = row[first]
            donor = row["sample"].strip()
            if donor not in set("ABCDEFGH"):
                raise SystemExit(f"unexpected donor {donor!r}")
            if cell in records:
                raise SystemExit(f"duplicate metadata cell {cell}")
            records[cell] = {
                "donor": donor,
                "age": float(row["age"]),
                "sex": row["sex"],
                "ethnicity": row["ethnicity"],
                "PMI_hours": float(row["PMI"]),
            }
            order.append(cell)
    return order, records


def parse_counts(
    path: Path, marker_genes: list[str]
) -> tuple[list[str], dict[str, np.ndarray], np.ndarray, int]:
    marker_set = set(marker_genes)
    marker_counts: dict[str, np.ndarray] = {}
    skipped_fnbp1_rows = 0
    with gzip.open(path, "rb") as handle:
        header = handle.readline().decode().strip().split()
        cells = [token.strip('"') for token in header]
        library_excluding_fnbp1 = np.zeros(len(cells), dtype=np.int64)
        for line_number, line in enumerate(handle, 2):
            if not line.startswith(b'"'):
                raise SystemExit(f"counts line {line_number} does not start with a quoted gene")
            closing = line.find(b'" ', 1)
            if closing < 0:
                raise SystemExit(f"counts line {line_number} has no gene delimiter")
            gene = line[1:closing].decode()
            if gene == "FNBP1":
                skipped_fnbp1_rows += 1
                continue
            values = np.fromstring(line[closing + 2 :], sep=" ", dtype=np.int64)
            if values.size != len(cells):
                raise SystemExit(
                    f"counts line {line_number} ({gene}) has {values.size} values; expected {len(cells)}"
                )
            library_excluding_fnbp1 += values
            if gene in marker_set:
                if gene in marker_counts:
                    marker_counts[gene] += values
                else:
                    marker_counts[gene] = values.copy()
    return cells, marker_counts, library_excluding_fnbp1, skipped_fnbp1_rows


def label_donor(
    cells: list[str],
    columns: np.ndarray,
    marker_counts: dict[str, np.ndarray],
    library: np.ndarray,
    modules: dict[str, dict[str, list[str]]],
    thresholds: dict[str, object],
) -> tuple[list[str], dict[str, object]]:
    donor_library = library[columns].astype(np.float64)
    if np.any(donor_library <= 0):
        raise SystemExit(f"{int((donor_library <= 0).sum())} donor nuclei have no non-FNBP1 counts")
    marker_genes = sorted(marker_counts)
    marker_slot = {gene: index for index, gene in enumerate(marker_genes)}
    counts = np.vstack([marker_counts[gene][columns] for gene in marker_genes]).astype(np.float64)
    transformed = np.log1p(counts * (10000.0 / donor_library)[None, :])
    marker_z = np.vstack([np.clip(zscore(row), -3.0, 3.0) for row in transformed])

    lineage_names = list(modules)
    lineage_scores: list[np.ndarray] = []
    lineage_detected: list[np.ndarray] = []
    for lineage in lineage_names:
        sub_scores: list[np.ndarray] = []
        sub_detected: list[np.ndarray] = []
        for genes in modules[lineage].values():
            slots = [marker_slot[gene] for gene in genes]
            sub_scores.append(zscore(marker_z[slots, :].mean(axis=0)))
            sub_detected.append((counts[slots, :] > 0).sum(axis=0))
        score_matrix = np.vstack(sub_scores)
        detected_matrix = np.vstack(sub_detected)
        selected = np.argmax(score_matrix, axis=0)
        donor_columns = np.arange(len(columns))
        lineage_scores.append(score_matrix[selected, donor_columns])
        lineage_detected.append(detected_matrix[selected, donor_columns])

    score_matrix = np.vstack(lineage_scores)
    detected_matrix = np.vstack(lineage_detected)
    order = np.argsort(score_matrix, axis=0)
    top_index, runner_index = order[-1, :], order[-2, :]
    donor_columns = np.arange(len(columns))
    top_score = score_matrix[top_index, donor_columns]
    runner_score = score_matrix[runner_index, donor_columns]
    top_detected = detected_matrix[top_index, donor_columns]
    runner_detected = detected_matrix[runner_index, donor_columns]
    margin = top_score - runner_score

    minimum_detected = int(thresholds["minimum_detected_markers_in_selected_submodule"])
    minimum_score = float(thresholds["minimum_lineage_score"])
    minimum_margin = float(thresholds["minimum_top_vs_runner_margin"])
    high_umi = float(np.quantile(donor_library, float(thresholds["high_GEX_UMI_quantile"])))
    labels: list[str] = []
    for cell in range(len(columns)):
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
            donor_library[cell] > high_umi
            and top_eligible
            and runner_eligible
            and runner_score[cell] >= 0.0
        )
        if conflict or high_count_conflict:
            labels.append("ambiguous_or_doublet")
        elif top_eligible and top_score[cell] >= minimum_score and margin[cell] >= minimum_margin:
            labels.append(top)
        else:
            labels.append("unresolved")

    detected_by_marker = {
        gene: int((marker_counts[gene][columns] > 0).sum()) for gene in marker_genes
    }
    summary: dict[str, object] = {
        "nuclei": len(columns),
        "groups": dict(sorted(Counter(labels).items())),
        "high_GEX_UMI_threshold_excluding_FNBP1": high_umi,
        "minimum_marker_positive_nuclei": min(detected_by_marker.values()),
        "marker_positive_nuclei": detected_by_marker,
    }
    return labels, summary


def read_ena(path: Path) -> dict[str, object]:
    per_donor: dict[str, dict[str, int]] = defaultdict(lambda: {"runs": 0, "fastq_bytes": 0, "read_pairs": 0})
    run_accessions: set[str] = set()
    experiments: set[str] = set()
    md5_complete = True
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            run_accessions.add(row["run_accession"])
            experiments.add(row["experiment_accession"])
            match = re.search(r", ([A-H]), snRNA-seq$", row["sample_title"])
            if not match:
                raise SystemExit(f"cannot resolve donor from ENA title {row['sample_title']!r}")
            donor = match.group(1)
            byte_fields = row["fastq_bytes"].split(";")
            md5_fields = row["fastq_md5"].split(";")
            ftp_fields = row["fastq_ftp"].split(";")
            if len(byte_fields) != 2 or len(md5_fields) != 2 or len(ftp_fields) != 2:
                md5_complete = False
            if not all(re.fullmatch(r"[0-9a-f]{32}", value) for value in md5_fields):
                md5_complete = False
            per_donor[donor]["runs"] += 1
            per_donor[donor]["fastq_bytes"] += sum(int(value) for value in byte_fields)
            per_donor[donor]["read_pairs"] += int(row["read_count"])
    return {
        "runs": len(run_accessions),
        "experiments": len(experiments),
        "paired_fastq_bytes": sum(row["fastq_bytes"] for row in per_donor.values()),
        "read_pairs": sum(row["read_pairs"] for row in per_donor.values()),
        "complete_paired_fastq_md5": md5_complete,
        "per_donor": dict(sorted(per_donor.items())),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--ena", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    if manifest["status"] != "PROSPECTIVE_GATE_FROZEN_BEFORE_ANY_FNBP1_SEQUENCE_OR_EVENT_QUERY":
        raise SystemExit("manifest is not a frozen pre-event gate")
    modules = manifest["marker_modules"]
    marker_genes = sorted(
        {gene for lineage in modules.values() for genes in lineage.values() for gene in genes}
    )
    if "FNBP1" in marker_genes:
        raise SystemExit("FNBP1 occurs in locked marker modules")
    processed = manifest["dataset"]["processed_inputs"]
    verify_input(args.counts, processed["counts"])
    verify_input(args.metadata, processed["metadata"])
    if sha256(args.ena) != processed["ena_run_manifest"]["sha256"]:
        raise SystemExit("ENA manifest digest mismatch")

    metadata_order, metadata = read_metadata(args.metadata)
    cells, marker_counts, library, skipped_fnbp1_rows = parse_counts(args.counts, marker_genes)
    if cells != metadata_order:
        raise SystemExit("counts header and metadata cell order differ")
    missing = sorted(set(marker_genes) - set(marker_counts))
    if missing:
        raise SystemExit(f"missing locked markers: {missing}")
    if skipped_fnbp1_rows != 1:
        raise SystemExit(f"expected to skip one opaque FNBP1 row; skipped {skipped_fnbp1_rows}")

    args.out.mkdir(parents=True, exist_ok=True)
    groups_dir = args.out / "groups"
    groups_dir.mkdir(exist_ok=True)
    donor_summaries: dict[str, dict[str, object]] = {}
    donor_ids = [entry["id"] for entry in manifest["dataset"]["donors"]]
    for donor in donor_ids:
        columns = np.asarray([index for index, cell in enumerate(cells) if metadata[cell]["donor"] == donor])
        labels, donor_summary = label_donor(
            cells, columns, marker_counts, library, modules, manifest["label_rule"]["thresholds"]
        )
        group_path = groups_dir / f"donor-{donor}.tsv"
        seen_barcodes: set[str] = set()
        with group_path.open("w") as handle:
            for column, label in zip(columns, labels):
                match = BARCODE_RE.search(cells[int(column)])
                if not match:
                    raise SystemExit(f"cannot extract 10x barcode from {cells[int(column)]}")
                barcode = match.group(1)
                if barcode in seen_barcodes:
                    raise SystemExit(f"duplicate donor-{donor} barcode {barcode}")
                seen_barcodes.add(barcode)
                handle.write(f"{barcode}\t{label}\n")
        donor_summary["group_map"] = str(group_path)
        donor_summary["group_map_sha256"] = sha256(group_path)
        donor_summaries[donor] = donor_summary

    ena = read_ena(args.ena)
    checks = {
        "eight_distinct_metadata_donors": sorted(donor_summaries) == donor_ids,
        "eight_distinct_sequence_experiments": ena["experiments"] == 8,
        "at_least_four_donors_with_1000_nuclei": sum(
            int(row["nuclei"]) >= 1000 for row in donor_summaries.values()
        ) >= 4,
        "every_locked_marker_present": not missing,
        "at_least_four_donors_with_50_microglia_immune_labels": sum(
            int(row["groups"].get("microglia_immune", 0)) >= 50 for row in donor_summaries.values()
        ) >= 4,
        "paired_fastq_md5_complete": ena["complete_paired_fastq_md5"],
    }
    result = {
        "schema": "gravlax.fnbp1-multidonor-sez-feasibility.v1",
        "status": "PASS" if all(checks.values()) else "STOP_UNDERPOWERED",
        "event_evidence_inspected": False,
        "fnbp1_counts_parsed": False,
        "opaque_fnbp1_rows_skipped": skipped_fnbp1_rows,
        "manifest": str(args.manifest),
        "manifest_sha256": sha256(args.manifest),
        "inputs": {
            "counts_sha256": sha256(args.counts),
            "metadata_sha256": sha256(args.metadata),
            "ena_manifest_sha256": sha256(args.ena),
        },
        "nuclei": len(cells),
        "donors": donor_summaries,
        "sequence": ena,
        "checks": checks,
        "next_action": (
            "raw sequence acquisition and fixed-event evaluation are authorized by the frozen gate"
            if all(checks.values())
            else "STOP: do not acquire or query raw sequence evidence"
        ),
        "software": {"numpy": np.__version__},
    }
    result_path = args.out / "feasibility.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
