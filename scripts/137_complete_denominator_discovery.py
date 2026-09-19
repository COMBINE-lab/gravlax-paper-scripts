#!/usr/bin/env python3
"""Evaluate prospective v32->v49 discovery with an explicit complete denominator."""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence


GENE_ID = re.compile(r'(?:^|;\s*)gene_id "([^"]+)"')
GENE_NAME = re.compile(r'(?:^|;\s*)gene_name "([^"]+)"')
GENE_TYPE = re.compile(r'(?:^|;\s*)gene_type "([^"]+)"')


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int
    strand: str
    name: str
    umis: int = 0
    cells: int = 0
    gene_name: str = ""
    gene_type: str = ""


class IntervalIndex:
    """Small immutable interval index using start order plus a prefix maximum end."""

    def __init__(self, rows: Iterable[Interval], stranded: bool = True):
        grouped: dict[tuple[str, str], list[Interval]] = defaultdict(list)
        for row in rows:
            grouped[(row.chrom, row.strand if stranded else "*")].append(row)
        self.groups: dict[tuple[str, str], tuple[list[int], list[int], list[Interval]]] = {}
        for key, vals in grouped.items():
            vals.sort(key=lambda x: (x.start, x.end, x.name))
            starts: list[int] = []
            prefix_max: list[int] = []
            high = -1
            for val in vals:
                starts.append(val.start)
                high = max(high, val.end)
                prefix_max.append(high)
            self.groups[key] = (starts, prefix_max, vals)
        self.stranded = stranded

    def overlaps(self, chrom: str, start: int, end: int, strand: str = "*") -> list[Interval]:
        key = (chrom, strand if self.stranded else "*")
        group = self.groups.get(key)
        if group is None:
            return []
        starts, prefix_max, vals = group
        i = bisect.bisect_left(starts, end) - 1
        found: list[Interval] = []
        while i >= 0 and prefix_max[i] > start:
            val = vals[i]
            if val.end > start:
                found.append(val)
            i -= 1
        found.sort(key=lambda x: (x.start, x.end, x.name))
        return found


def stable_gene_id(value: str) -> str:
    # Preserve GENCODE's pseudoautosomal ``_PAR_Y`` locus suffix while removing
    # only the numeric version (ENSG....7_PAR_Y -> ENSG..._PAR_Y).
    return re.sub(r"\.\d+(?=_PAR_Y$|$)", "", value)


def normalize_chrom(value: str) -> str:
    if value.startswith("chr"):
        return value
    if value == "MT":
        return "chrM"
    if value.isdigit() or value in {"X", "Y", "M"}:
        return "chr" + value
    return value


def parse_gtf_genes(path: Path) -> dict[str, Interval]:
    genes: dict[str, Interval] = {}
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or "\tgene\t" not in line:
                continue
            fields = line.rstrip("\n").split("\t", 8)
            if len(fields) != 9 or fields[2] != "gene":
                continue
            match = GENE_ID.search(fields[8])
            if match is None:
                raise ValueError(f"gene record without gene_id in {path}")
            gene_id = stable_gene_id(match.group(1))
            name = GENE_NAME.search(fields[8])
            kind = GENE_TYPE.search(fields[8])
            row = Interval(
                normalize_chrom(fields[0]),
                int(fields[3]) - 1,
                int(fields[4]),
                fields[6],
                gene_id,
                gene_name=name.group(1) if name else "",
                gene_type=kind.group(1) if kind else "",
            )
            previous = genes.get(gene_id)
            if previous is not None and previous != row:
                raise ValueError(f"duplicate stable gene_id with unequal records: {gene_id}")
            genes[gene_id] = row
    return genes


def read_called_barcodes(path: Path) -> set[str]:
    with path.open() as handle:
        return {line.rstrip("\n") for line in handle if line.strip()}


def selected_columns(raw_barcodes: Path, called: set[str]) -> set[int]:
    selected: set[int] = set()
    with raw_barcodes.open() as handle:
        for column, line in enumerate(handle, start=1):
            if line.rstrip("\n") in called:
                selected.add(column)
    if len(selected) != len(called):
        raise ValueError(
            f"called/raw barcode mismatch for {raw_barcodes}: {len(selected)} of {len(called)}"
        )
    return selected


def matrix_gene_sums(raw_dir: Path, filtered_dir: Path) -> dict[str, int]:
    called = read_called_barcodes(filtered_dir / "barcodes.tsv")
    columns = selected_columns(raw_dir / "barcodes.tsv", called)
    features: list[str] = []
    with (raw_dir / "features.tsv").open() as handle:
        for line in handle:
            features.append(stable_gene_id(line.split("\t", 1)[0]))
    sums = [0] * len(features)
    dimensions_seen = False
    with (raw_dir / "matrix.mtx").open() as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.split()
            if not dimensions_seen:
                if len(fields) != 3 or int(fields[0]) != len(features):
                    raise ValueError(f"unexpected MatrixMarket dimensions in {raw_dir}")
                dimensions_seen = True
                continue
            row, column, value = int(fields[0]), int(fields[1]), float(fields[2])
            if column in columns:
                if not value.is_integer():
                    raise ValueError(f"non-integral Gene count in {raw_dir}: {value}")
                sums[row - 1] += int(value)
    if not dimensions_seen:
        raise ValueError(f"missing MatrixMarket dimensions in {raw_dir}")
    result: dict[str, int] = defaultdict(int)
    for gene, value in zip(features, sums):
        result[gene] += value
    return dict(result)


def parse_discovery(path: Path, prefix: str) -> list[Interval]:
    rows: list[Interval] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"{path}:{line_number}: expected six columns")
            chrom, start, end, strand, umis, cells = fields
            rows.append(
                Interval(
                    normalize_chrom(chrom), int(start), int(end), strand,
                    f"{prefix}_{line_number:06d}", int(umis), int(cells)
                )
            )
    return rows


def parse_polyasite(path: Path) -> dict[tuple[str, str], list[int]]:
    sites: dict[tuple[str, str], list[int]] = defaultdict(list)
    import gzip
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[5] not in {"+", "-"}:
                continue
            point = (int(fields[1]) + int(fields[2])) // 2
            sites[(normalize_chrom(fields[0]), fields[5])].append(point)
    for values in sites.values():
        values.sort()
    return sites


def polyasite_support(sites: dict[tuple[str, str], list[int]], row: Interval) -> bool:
    values = sites.get((row.chrom, row.strand), [])
    if row.strand == "+":
        start, end = row.start, row.end + 500
    else:
        start, end = max(0, row.start - 500), row.end
    i = bisect.bisect_left(values, start)
    return i < len(values) and values[i] < end


def parse_star_junctions(path: Path) -> IntervalIndex:
    rows: list[Interval] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) < 9:
                raise ValueError(f"{path}:{line_number}: malformed STAR junction")
            strand = {1: "+", 2: "-"}.get(int(fields[3]))
            unique_reads = int(fields[6])
            if strand is None or unique_reads < 3:
                continue
            rows.append(
                Interval(
                    normalize_chrom(fields[0]), int(fields[1]) - 1, int(fields[2]),
                    strand, f"SJ_{line_number}", unique_reads
                )
            )
    return IntervalIndex(rows)


def star_splice_support(index: IntervalIndex, row: Interval) -> bool:
    return any(
        junction.start >= row.start and junction.end <= row.end
        for junction in index.overlaps(row.chrom, row.start, row.end, row.strand)
    )


def matching(
    truth: Sequence[Interval], candidates: Sequence[Interval]
) -> tuple[dict[int, int], dict[int, list[int]]]:
    candidate_index = IntervalIndex(candidates)
    by_name = {candidate.name: i for i, candidate in enumerate(candidates)}
    adjacency: dict[int, list[int]] = {}
    for gene_index, gene in enumerate(truth):
        overlaps = candidate_index.overlaps(gene.chrom, gene.start, gene.end, gene.strand)
        overlaps.sort(
            key=lambda candidate: (
                -min(gene.end, candidate.end) + max(gene.start, candidate.start),
                -candidate.umis,
                candidate.chrom,
                candidate.start,
                candidate.end,
                candidate.name,
            )
        )
        adjacency[gene_index] = [by_name[candidate.name] for candidate in overlaps]

    candidate_to_gene: dict[int, int] = {}

    def augment(gene_index: int, seen: set[int]) -> bool:
        for candidate in adjacency[gene_index]:
            if candidate in seen:
                continue
            seen.add(candidate)
            previous = candidate_to_gene.get(candidate)
            if previous is None or augment(previous, seen):
                candidate_to_gene[candidate] = gene_index
                return True
        return False

    order = sorted(
        range(len(truth)),
        key=lambda i: (len(adjacency[i]), truth[i].chrom, truth[i].start, truth[i].end, truth[i].name),
    )
    sys.setrecursionlimit(max(2000, len(truth) * 4))
    for gene_index in order:
        augment(gene_index, set())
    gene_to_candidate = {gene: candidate for candidate, gene in candidate_to_gene.items()}
    return gene_to_candidate, adjacency


def verdict(value: float, pass_at: float, marginal_at: float) -> str:
    if value >= pass_at:
        return "pass"
    if value >= marginal_at:
        return "marginal"
    return "fail"


def ratio(numerator: int | float, denominator: int | float) -> float:
    return float(numerator) / float(denominator) if denominator else 0.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(root: Path, path: Path) -> dict[str, int | str]:
    try:
        logical = str(path.relative_to(root))
    except ValueError:
        logical = path.name
    return {"path": logical, "bytes": path.stat().st_size, "sha256": sha256(path)}


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize_threshold(
    threshold: int,
    new_genes: dict[str, Interval],
    d0_counts: dict[str, int],
    strata: dict[str, str],
    candidates: list[Interval],
) -> tuple[dict[str, object], dict[int, int], dict[int, list[int]], list[Interval]]:
    truth = sorted(
        (gene for gene_id, gene in new_genes.items() if d0_counts.get(gene_id, 0) >= threshold),
        key=lambda x: (x.chrom, x.start, x.end, x.strand, x.name),
    )
    matches, adjacency = matching(truth, candidates)
    by_stratum: dict[str, dict[str, object]] = {}
    for stratum in ("clean", "antisense_only", "same_strand"):
        indices = [i for i, gene in enumerate(truth) if strata[gene.name] == stratum]
        hit = [i for i in indices if i in matches]
        total_mass = sum(d0_counts[truth[i].name] for i in indices)
        hit_mass = sum(d0_counts[truth[i].name] for i in hit)
        by_stratum[stratum] = {
            "truth_genes": len(indices),
            "matched_genes": len(hit),
            "recall": ratio(len(hit), len(indices)),
            "truth_umi_mass": total_mass,
            "matched_umi_mass": hit_mass,
            "umi_weighted_recall": ratio(hit_mass, total_mass),
            "permissive_any_overlap_genes": sum(bool(adjacency[i]) for i in indices),
            "permissive_any_overlap_recall": ratio(sum(bool(adjacency[i]) for i in indices), len(indices)),
        }
    total_mass = sum(d0_counts[gene.name] for gene in truth)
    hit_mass = sum(d0_counts[truth[i].name] for i in matches)
    result = {
        "minimum_d0_gene_umis": threshold,
        "truth_genes": len(truth),
        "matched_genes": len(matches),
        "unmatched_genes": len(truth) - len(matches),
        "macro_recall": ratio(len(matches), len(truth)),
        "truth_umi_mass": total_mass,
        "matched_umi_mass": hit_mass,
        "umi_weighted_recall": ratio(hit_mass, total_mass),
        "permissive_any_overlap_genes": sum(bool(adjacency[i]) for i in range(len(truth))),
        "permissive_any_overlap_recall": ratio(sum(bool(adjacency[i]) for i in range(len(truth))), len(truth)),
        "strata": by_stratum,
    }
    return result, matches, adjacency, truth


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--d0-discovery", type=Path, required=True)
    parser.add_argument("--d1-discovery", type=Path, required=True)
    parser.add_argument("--aie-bin", type=Path, required=True)
    parser.add_argument("--code-commit", required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--truth-tsv", type=Path, required=True)
    parser.add_argument("--candidate-tsv", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()

    paths = {
        "gencode_v32": root / "annotations/gencode.v32.annotation.gtf",
        "gencode_v49": root / "annotations/gencode.v49.annotation.gtf",
        "polyasite": root / "runs/apa2/polyasite.bed.gz",
        "d0_star_junctions": root / "runs/oracle/full/v49/SJ.out.tab",
        "d0_discovery": args.d0_discovery.resolve(),
        "d1_discovery": args.d1_discovery.resolve(),
        "d0_gene_matrix": root / "runs/oracle/full/v49/Solo.out/Gene/raw/matrix.mtx",
        "d0_gene_features": root / "runs/oracle/full/v49/Solo.out/Gene/raw/features.tsv",
        "d0_raw_barcodes": root / "runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv",
        "d0_filtered_barcodes": root / "runs/oracle/full/v49/Solo.out/Gene/filtered/barcodes.tsv",
        "d1_gene_matrix": root / "runs/oracle/d1/v49/Solo.out/Gene/raw/matrix.mtx",
        "d1_gene_features": root / "runs/oracle/d1/v49/Solo.out/Gene/raw/features.tsv",
        "d1_raw_barcodes": root / "runs/oracle/d1/v49/Solo.out/Gene/raw/barcodes.tsv",
        "d1_filtered_barcodes": root / "runs/oracle/d1/v49/Solo.out/Gene/filtered/barcodes.tsv",
        "aie_binary": args.aie_bin.resolve(),
    }

    v32 = parse_gtf_genes(paths["gencode_v32"])
    v49 = parse_gtf_genes(paths["gencode_v49"])
    new_genes = {gene_id: gene for gene_id, gene in v49.items() if gene_id not in v32}
    v32_unstranded = IntervalIndex(v32.values(), stranded=False)
    strata: dict[str, str] = {}
    for gene_id, gene in new_genes.items():
        overlaps = v32_unstranded.overlaps(gene.chrom, gene.start, gene.end)
        if not overlaps:
            strata[gene_id] = "clean"
        elif any(other.strand == gene.strand for other in overlaps):
            strata[gene_id] = "same_strand"
        else:
            strata[gene_id] = "antisense_only"

    d0_raw = root / "runs/oracle/full/v49/Solo.out/Gene/raw"
    d0_filtered = root / "runs/oracle/full/v49/Solo.out/Gene/filtered"
    d1_raw = root / "runs/oracle/d1/v49/Solo.out/Gene/raw"
    d1_filtered = root / "runs/oracle/d1/v49/Solo.out/Gene/filtered"
    d0_counts = matrix_gene_sums(d0_raw, d0_filtered)
    d1_counts = matrix_gene_sums(d1_raw, d1_filtered)
    d0_candidates = parse_discovery(paths["d0_discovery"], "D0")
    d1_candidates = parse_discovery(paths["d1_discovery"], "D1")
    d1_index = IntervalIndex(d1_candidates)
    polyasite = parse_polyasite(paths["polyasite"])
    star = parse_star_junctions(paths["d0_star_junctions"])

    threshold_results: dict[str, dict[str, object]] = {}
    primary_matches: dict[int, int] = {}
    primary_adjacency: dict[int, list[int]] = {}
    primary_truth: list[Interval] = []
    for threshold in (10, 50, 100):
        result, matches, adjacency, truth = summarize_threshold(
            threshold, new_genes, d0_counts, strata, d0_candidates
        )
        threshold_results[str(threshold)] = result
        if threshold == 50:
            primary_matches, primary_adjacency, primary_truth = matches, adjacency, truth

    candidate_to_truth = {candidate: gene for gene, candidate in primary_matches.items()}
    primary_truth_names = {gene.name for gene in primary_truth}
    all_new_index = IntervalIndex(new_genes.values())
    v49_preexisting_index = IntervalIndex(gene for gene_id, gene in v49.items() if gene_id in v32)

    truth_rows: list[dict[str, object]] = []
    for gene_index, gene in enumerate(primary_truth):
        candidate_index = primary_matches.get(gene_index)
        candidate = d0_candidates[candidate_index] if candidate_index is not None else None
        d1_hits = d1_index.overlaps(gene.chrom, gene.start, gene.end, gene.strand)
        truth_rows.append({
            "gene_id": gene.name,
            "gene_name": gene.gene_name,
            "gene_type": gene.gene_type,
            "chrom": gene.chrom,
            "start": gene.start,
            "end": gene.end,
            "strand": gene.strand,
            "d0_gene_umis": d0_counts.get(gene.name, 0),
            "d1_gene_umis": d1_counts.get(gene.name, 0),
            "stratum": strata[gene.name],
            "matched": int(candidate is not None),
            "matched_candidate": candidate.name if candidate else "",
            "candidate_umis": candidate.umis if candidate else 0,
            "permissive_any_overlap": int(bool(primary_adjacency[gene_index])),
            "d1_expression": int(d1_counts.get(gene.name, 0) >= 10),
            "d1_discovery": int(bool(d1_hits)),
            "polyasite": int(polyasite_support(polyasite, gene)),
            "star_splice": int(star_splice_support(star, gene)),
        })

    candidate_rows: list[dict[str, object]] = []
    for candidate_index, candidate in enumerate(d0_candidates):
        matched_gene_index = candidate_to_truth.get(candidate_index)
        matched_gene = primary_truth[matched_gene_index] if matched_gene_index is not None else None
        new_overlaps = all_new_index.overlaps(candidate.chrom, candidate.start, candidate.end, candidate.strand)
        old_overlaps = v49_preexisting_index.overlaps(
            candidate.chrom, candidate.start, candidate.end, candidate.strand
        )
        if matched_gene is not None:
            classification = "primary_match"
        elif new_overlaps:
            classification = "other_v49_added_gene"
        elif old_overlaps:
            classification = "preexisting_v49_gene_span"
        else:
            classification = "unresolved"
        d1_hit = bool(d1_index.overlaps(candidate.chrom, candidate.start, candidate.end, candidate.strand))
        pas_hit = polyasite_support(polyasite, candidate)
        sj_hit = star_splice_support(star, candidate)
        candidate_rows.append({
            "candidate_id": candidate.name,
            "chrom": candidate.chrom,
            "start": candidate.start,
            "end": candidate.end,
            "strand": candidate.strand,
            "d0_umis": candidate.umis,
            "d0_cells": candidate.cells,
            "annotation_history_class": classification,
            "matched_primary_gene": matched_gene.name if matched_gene else "",
            "overlapping_v49_added_genes": ",".join(gene.name for gene in new_overlaps),
            "d1_discovery": int(d1_hit),
            "polyasite": int(pas_hit),
            "star_splice": int(sj_hit),
            "support_count": int(d1_hit) + int(pas_hit) + int(sj_hit),
        })

    class_counts = Counter(row["annotation_history_class"] for row in candidate_rows)
    class_mass: Counter[str] = Counter()
    support_combinations: Counter[str] = Counter()
    support_mass: Counter[str] = Counter()
    for row in candidate_rows:
        cls = str(row["annotation_history_class"])
        class_mass[cls] += int(row["d0_umis"])
        combination = "".join(
            code for code, field in (("D", "d1_discovery"), ("P", "polyasite"), ("S", "star_splice"))
            if row[field]
        ) or "none"
        support_combinations[combination] += 1
        support_mass[combination] += int(row["d0_umis"])

    ranked = sorted(
        candidate_rows,
        key=lambda row: (
            -int(row["d0_umis"]), str(row["chrom"]), int(row["start"]),
            int(row["end"]), str(row["strand"]), str(row["candidate_id"])
        ),
    )
    labels = [row["annotation_history_class"] == "primary_match" for row in ranked]
    missing = len(primary_truth) - sum(labels)
    labels.extend([True] * missing)
    positives_seen = 0
    ap_sum = 0.0
    for rank, label in enumerate(labels, start=1):
        if label:
            positives_seen += 1
            ap_sum += positives_seen / rank
    average_precision = ap_sum / len(primary_truth) if primary_truth else 0.0
    precision_at: dict[str, dict[str, float | int]] = {}
    for k in (100, 500, 1000):
        observed = ranked[:k]
        positives = sum(row["annotation_history_class"] == "primary_match" for row in observed)
        precision_at[str(k)] = {
            "rows": len(observed), "positives": positives, "precision": ratio(positives, len(observed))
        }

    matched_candidates = len(primary_matches)
    total_candidate_mass = sum(candidate.umis for candidate in d0_candidates)
    recurrent = [row for row in candidate_rows if row["d1_discovery"]]
    recurrent_mass = sum(int(row["d0_umis"]) for row in recurrent)
    primary = threshold_results["50"]
    clean_recall = float(primary["strata"]["clean"]["recall"])
    macro_recall = float(primary["macro_recall"])

    if sum(int(value["truth_genes"]) for value in primary["strata"].values()) != len(primary_truth):
        raise AssertionError("truth stratum accounting failed")
    if int(primary["matched_genes"]) + int(primary["unmatched_genes"]) != len(primary_truth):
        raise AssertionError("matched/unmatched truth accounting failed")
    if sum(class_counts.values()) != len(d0_candidates):
        raise AssertionError("candidate class accounting failed")
    if len(set(candidate_to_truth)) != len(primary_matches):
        raise AssertionError("matching is not one-to-one")

    def flag_summary(rows: list[dict[str, object]], fields: tuple[str, ...]) -> dict[str, object]:
        return {
            field: {
                "rows": sum(bool(row[field]) for row in rows),
                "fraction": ratio(sum(bool(row[field]) for row in rows), len(rows)),
            }
            for field in fields
        }

    truth_support = {
        "all": flag_summary(truth_rows, ("d1_expression", "d1_discovery", "polyasite", "star_splice")),
        "by_match": {
            label: flag_summary(
                [row for row in truth_rows if bool(row["matched"]) == matched],
                ("d1_expression", "d1_discovery", "polyasite", "star_splice"),
            )
            for label, matched in (("matched", True), ("unmatched", False))
        },
        "by_stratum": {
            stratum: flag_summary(
                [row for row in truth_rows if row["stratum"] == stratum],
                ("d1_expression", "d1_discovery", "polyasite", "star_splice"),
            )
            for stratum in ("clean", "antisense_only", "same_strand")
        },
    }
    candidate_recurrence_by_class = {}
    for classification in sorted(class_counts):
        rows = [row for row in candidate_rows if row["annotation_history_class"] == classification]
        recurrence_rows = [row for row in rows if row["d1_discovery"]]
        mass = sum(int(row["d0_umis"]) for row in rows)
        class_recurrence_mass = sum(int(row["d0_umis"]) for row in recurrence_rows)
        candidate_recurrence_by_class[classification] = {
            "loci": len(rows),
            "recurrent_loci": len(recurrence_rows),
            "recurrence_fraction": ratio(len(recurrence_rows), len(rows)),
            "d0_umi_mass": mass,
            "recurrent_d0_umi_mass": class_recurrence_mass,
            "recurrence_umi_fraction": ratio(class_recurrence_mass, mass),
        }

    result = {
        "schema": "gravlax.complete-denominator-discovery.v1",
        "protocol_lock_commit": "522bad43daf8d22c5141d7c5b4eb3f8fe9ca6563",
        "gravlax_code_commit": args.code_commit,
        "configuration": {
            "primary_min_d0_gene_umis": 50,
            "sensitivity_thresholds": [10, 100],
            "discovery_min_umis": 10,
            "discovery_merge_gap_bp": 1000,
            "matching": "same-chromosome/strand interval overlap; deterministic maximum-cardinality one-to-one",
            "d1_expression_min_umis": 10,
            "polyasite_downstream_bp": 500,
            "star_min_unique_reads": 3,
        },
        "inputs": {name: identity(root, path) for name, path in sorted(paths.items())},
        "annotation_history": {
            "v32_genes": len(v32),
            "v49_genes": len(v49),
            "v49_added_stable_gene_ids": len(new_genes),
        },
        "thresholds": threshold_results,
        "primary_gate": {
            "clean_recall_verdict": verdict(clean_recall, 0.90, 0.75),
            "all_strata_macro_recall_verdict": verdict(macro_recall, 0.60, 0.35),
            "accounting": "pass",
        },
        "annotation_history_confirmation": {
            "emitted_d0_candidates": len(d0_candidates),
            "matched_primary_truth": matched_candidates,
            "precision_at_emission_floor": ratio(matched_candidates, len(d0_candidates)),
            "recall_at_emission_floor": ratio(matched_candidates, len(primary_truth)),
            "ranked_event_average_precision": average_precision,
            "precision_at": precision_at,
            "candidate_classes": {
                key: {"loci": class_counts[key], "d0_umi_mass": class_mass[key]}
                for key in sorted(class_counts)
            },
        },
        "independent_support": {
            "d1_discovery_candidates": len(d1_candidates),
            "d0_candidate_d1_recurrence_loci": len(recurrent),
            "d0_candidate_d1_recurrence_fraction": ratio(len(recurrent), len(d0_candidates)),
            "d0_candidate_d1_recurrence_umi_mass": recurrent_mass,
            "d0_candidate_d1_recurrence_umi_fraction": ratio(recurrent_mass, total_candidate_mass),
            "d0_candidate_d1_recurrence_by_annotation_history_class": candidate_recurrence_by_class,
            "primary_truth_support": truth_support,
            "support_combinations": {
                key: {"loci": support_combinations[key], "d0_umi_mass": support_mass[key]}
                for key in sorted(support_combinations)
            },
        },
        "accounting": {
            "truth_rows": len(truth_rows),
            "candidate_rows": len(candidate_rows),
            "truth_strata_sum": sum(int(value["truth_genes"]) for value in primary["strata"].values()),
            "truth_matched_plus_unmatched": int(primary["matched_genes"]) + int(primary["unmatched_genes"]),
            "candidate_class_sum": sum(class_counts.values()),
        },
        "interpretation_guard": (
            "Precision is annotation-history confirmation only; unsupported loci are unresolved, "
            "not established biological false positives."
        ),
    }

    write_tsv(args.truth_tsv, list(truth_rows[0]) if truth_rows else [], truth_rows)
    write_tsv(args.candidate_tsv, list(candidate_rows[0]) if candidate_rows else [], candidate_rows)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(
        f"primary truth {len(primary_truth)}; matched {matched_candidates}; "
        f"clean recall {clean_recall:.1%}; all-strata recall {macro_recall:.1%}; "
        f"D1 recurrence {ratio(len(recurrent), len(d0_candidates)):.1%} loci"
    )


if __name__ == "__main__":
    main()
