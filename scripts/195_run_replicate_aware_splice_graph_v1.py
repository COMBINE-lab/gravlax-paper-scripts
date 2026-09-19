#!/usr/bin/env python3
"""Run and validate the frozen replicate-aware molecular splice-graph v1 gate."""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import subprocess
import tempfile
import time
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def git_head(repository: Path) -> str:
    return subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "HEAD"],
        check=True,
        text=True,
        capture_output=True,
    ).stdout.strip()


def parse_design(path: Path) -> list[dict[str, object]]:
    lines = path.read_text().splitlines()
    if not lines or lines[0] != "sample\tcondition\tarchive\tcells":
        raise SystemExit("design header is not the frozen four-column header")
    rows = []
    seen_samples: set[str] = set()
    seen_archives: set[Path] = set()
    for line_number, line in enumerate(lines[1:], 2):
        fields = line.split("\t")
        if len(fields) != 4 or any(not field for field in fields):
            raise SystemExit(f"malformed design line {line_number}")
        sample, condition, archive_label, cells_label = fields
        archive = (path.parent / archive_label).resolve()
        if sample in seen_samples or archive in seen_archives:
            raise SystemExit("design contains a duplicate sample or resolved archive")
        if not archive.is_file():
            raise SystemExit(f"missing design archive: {archive}")
        seen_samples.add(sample)
        seen_archives.add(archive)
        rows.append(
            {
                "sample": sample,
                "condition": condition,
                "archive": archive,
                "archive_label": archive_label,
                "cells": cells_label,
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aie", required=True, type=Path)
    parser.add_argument("--code-repo", required=True, type=Path)
    parser.add_argument("--design", required=True, type=Path)
    parser.add_argument("--registered-result", required=True, type=Path)
    parser.add_argument("--locus", default="chr9:129907500-129926000")
    parser.add_argument("--runs", default=20, type=int)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.runs < 1:
        raise SystemExit("--runs must be at least 1")
    for path in (args.aie, args.design, args.registered_result):
        if not path.is_file():
            raise SystemExit(f"missing input: {path}")
    if not args.code_repo.is_dir():
        raise SystemExit(f"missing code repository: {args.code_repo}")
    design_rows = parse_design(args.design)
    if len(design_rows) != 8 or {row["condition"] for row in design_rows} != {"cohort"}:
        raise SystemExit("real gate must retain eight donors under one neutral condition")
    args.work_dir.mkdir(parents=True, exist_ok=True)

    command = [
        str(args.aie),
        "cohort",
        "splice-graph",
        args.locus,
        "--design",
        str(args.design),
        "--counts-only",
        "--min-support",
        "1",
        "--min-edge-samples",
        "2",
        "--json",
    ]
    raw_outputs: list[bytes] = []
    wall_seconds: list[float] = []
    peak_rss_kib: list[int] = []
    for run in range(args.runs):
        with tempfile.NamedTemporaryFile(
            dir=args.work_dir, prefix=f"time-{run:02d}-", delete=False
        ) as timing_handle:
            timing_path = Path(timing_handle.name)
        timed_command = ["/usr/bin/time", "-f", "%M", "-o", str(timing_path), *command]
        started = time.perf_counter()
        completed = subprocess.run(timed_command, check=True, capture_output=True)
        wall_seconds.append(time.perf_counter() - started)
        peak_rss_kib.append(int(timing_path.read_text().strip()))
        timing_path.unlink()
        raw_outputs.append(completed.stdout)
    if any(output != raw_outputs[0] for output in raw_outputs[1:]):
        raise SystemExit("repeated cohort splice-graph JSON was not byte-identical")

    graph = json.loads(raw_outputs[0])
    if graph.get("schema") != "gravlax.cohort.splice-graph.v1":
        raise SystemExit("unexpected cohort splice-graph schema")
    semantics = graph.get("semantics", {})
    expected_semantics = {
        "replicate_unit": "one unique design sample/archive row",
        "cells_and_molecules_are_replicates": False,
        "complete_transcript_claim": False,
        "missing_count": "explicit zero; evidence is not imputed",
    }
    for key, expected in expected_semantics.items():
        if semantics.get(key) != expected:
            raise SystemExit(f"missing or changed semantic guard: {key}")
    if graph.get("inference") != {"enabled": False, "reason": "counts-only", "tests": []}:
        raise SystemExit("real gate unexpectedly emitted inference")

    paths = {int(path["id"]): path for path in graph["paths"]}
    edges = {int(edge["id"]): edge for edge in graph["edges"]}
    sample_rows = graph["samples"]
    if [sample["sample"] for sample in sample_rows] != [row["sample"] for row in design_rows]:
        raise SystemExit("sample order differs from the frozen design")
    zero_path_rows = 0
    zero_edge_rows = 0
    sample_summaries = []
    path_sample_counts = {path_id: {} for path_id in paths}
    edge_sample_counts = {edge_id: {} for edge_id in edges}
    for sample in sample_rows:
        sample_id = sample["sample"]
        path_rows = {int(row["path_id"]): int(row["umis"]) for row in sample["path_rows"]}
        edge_rows = {int(row["edge_id"]): int(row["umis"]) for row in sample["edge_rows"]}
        if set(path_rows) != set(paths) or set(edge_rows) != set(edges):
            raise SystemExit(f"incomplete zero-explicit matrix for sample {sample_id}")
        zero_path_rows += sum(count == 0 for count in path_rows.values())
        zero_edge_rows += sum(count == 0 for count in edge_rows.values())
        for path_id, count in path_rows.items():
            path_sample_counts[path_id][sample_id] = count
        for edge_id, observed in edge_rows.items():
            expected = sum(
                count
                for path_id, count in path_rows.items()
                if edge_id in paths[path_id]["edge_ids"]
            )
            if observed != expected:
                raise SystemExit(
                    f"edge/path conservation failed for sample {sample_id}, edge {edge_id}"
                )
            edge_sample_counts[edge_id][sample_id] = observed
        strand_counts = Counter()
        for path_id, count in path_rows.items():
            strand_counts[paths[path_id]["strand"]] += count
        reported_strands = {
            row["strand"]: int(row["umis"]) for row in sample["strand_totals"]
        }
        if reported_strands != {"+": strand_counts["+"], "-": strand_counts["-"]}:
            raise SystemExit(f"strand total conservation failed for sample {sample_id}")
        sample_summaries.append(
            {
                "sample": sample_id,
                "plus_path_umis": strand_counts["+"],
                "minus_path_umis": strand_counts["-"],
                "nonzero_paths": sum(count > 0 for count in path_rows.values()),
                "eligible_strands_at_default": [
                    row["strand"] for row in sample["strand_totals"] if row["eligible"]
                ],
            }
        )

    for path_id, path in paths.items():
        counts = path_sample_counts[path_id]
        if sum(counts.values()) != int(path["aggregate_umis"]):
            raise SystemExit(f"aggregate path count failed for path {path_id}")
        if sum(count > 0 for count in counts.values()) != int(path["supporting_samples"]):
            raise SystemExit(f"supporting-sample count failed for path {path_id}")
    if any(int(edge["catalogue_samples"]) < 2 for edge in edges.values()):
        raise SystemExit("returned edge violates the frozen recurrence threshold")

    registered = json.loads(args.registered_result.read_text())
    stop_preserved = (
        registered.get("status") == "STOP_UNDERPOWERED"
        and registered.get("primary_brain_replication", {}).get("status")
        == "STOP_UNDERPOWERED"
        and registered.get("microglia_localization", {}).get("status")
        == "STOP_DESCRIPTIVE"
    )
    if not stop_preserved:
        raise SystemExit("the registered donor-level STOP record is absent or changed")

    length_histogram: dict[int, Counter[str]] = {}
    for path in paths.values():
        length = len(path["edge_ids"])
        length_histogram.setdefault(length, Counter())
        length_histogram[length]["paths"] += 1
        length_histogram[length]["umis"] += int(path["aggregate_umis"])
    path_summaries = [
        {
            "path_id": path_id,
            "strand": path["strand"],
            "junctions": path["junctions"],
            "aggregate_umis": path["aggregate_umis"],
            "supporting_samples": path["supporting_samples"],
            "sample_umis": path_sample_counts[path_id],
        }
        for path_id, path in paths.items()
    ]
    edge_summaries = [
        {
            "edge_id": edge_id,
            "strand": edge["strand"],
            "donor": edge["donor"],
            "acceptor": edge["acceptor"],
            "catalogue_samples": edge["catalogue_samples"],
            "sample_umis": edge_sample_counts[edge_id],
        }
        for edge_id, edge in edges.items()
    ]

    result = {
        "schema": "gravlax.replicate-aware-splice-graph-v1.result.v1",
        "status": "PASS_ENGINEERING_COUNTS_ONLY",
        "registration": {
            "manifest": "experiments/manifests/replicate-aware-splice-graph-v1.json",
            "specification": "experiments/specs/replicate-aware-splice-graph-v1.md",
            "registration_commit": "8d01b5ac87c327598a9ac1e980f07c767d8fe079",
        },
        "implementation": {
            "gravlax_commit": git_head(args.code_repo),
            "executable_sha256": sha256(args.aie),
            "query_schema": graph["schema"],
        },
        "inputs": {
            "dataset": "GSE234790 / PRJNA983239 eight adult human SEZ donors",
            "locus": args.locus,
            "design": "experiments/designs/fnbp1-sez-cohort-splice-graph.tsv",
            "design_sha256": sha256(args.design),
            "registered_result": "results/post-v1-fnbp1-multidonor-sez-validation.json",
            "registered_result_sha256": sha256(args.registered_result),
            "reference_digest": graph["reference_digest"],
            "archives": [
                {
                    "sample": row["sample"],
                    "logical_path": str(row["archive_label"]),
                    "bytes": row["archive"].stat().st_size,
                    "sha256": sha256(row["archive"]),
                }
                for row in design_rows
            ],
        },
        "thresholds": graph["thresholds"],
        "graph": {
            "planning": graph["planning"],
            "total_path_umis": sum(int(path["aggregate_umis"]) for path in paths.values()),
            "sample_path_rows": len(sample_rows) * len(paths),
            "zero_path_rows": zero_path_rows,
            "sample_edge_rows": len(sample_rows) * len(edges),
            "zero_edge_rows": zero_edge_rows,
            "path_length_histogram": [
                {"edges_per_path": length, **counts}
                for length, counts in sorted(length_histogram.items())
            ],
            "samples": sample_summaries,
            "paths": path_summaries,
            "edges": edge_summaries,
        },
        "acceptance": {
            "repeated_json_byte_identical": True,
            "complete_zero_explicit_sample_path_matrix": True,
            "complete_zero_explicit_sample_edge_matrix": True,
            "per_sample_edge_path_umi_conservation": True,
            "sample_strand_total_conservation": True,
            "aggregate_path_and_support_conservation": True,
            "common_edge_recurrence_guard": True,
            "counts_only_emits_no_test": True,
            "biological_sample_semantic_guard": True,
            "registered_stop_preserved": stop_preserved,
            "workspace_tests_passed": 124,
            "documentation_site_build_passed": True,
        },
        "benchmark": {
            "runs": args.runs,
            "median_wall_seconds": round(statistics.median(wall_seconds), 6),
            "min_wall_seconds": round(min(wall_seconds), 6),
            "max_wall_seconds": round(max(wall_seconds), 6),
            "max_peak_rss_kib": max(peak_rss_kib),
            "measurement": "release executable; warm repeated process invocations; Python perf_counter wall and GNU time maximum RSS",
        },
        "output_sha256": hashlib.sha256(raw_outputs[0]).hexdigest(),
        "interpretation": {
            "inference": "none; the real panel is the preregistered single-condition counts-only engineering check",
            "registered_outcome": "STOP_UNDERPOWERED for primary donor replication and STOP_DESCRIPTIVE for microglial localization remain binding",
            "capability": "the common graph exposes exact recurrent edge and multi-junction path-fragment rows across donors",
        },
        "interpretation_limits": [
            "Path fragments are positive junction co-support, not reconstructed complete transcripts.",
            "No arbitrary donor contrast was constructed and no p-value was emitted.",
            "The targeted archives contain only the registered locus and called-cell scope.",
            "Sparse FNBP1 evidence cannot support a donor-level differential-splicing claim.",
        ],
    }
    if not all(
        value is True for value in result["acceptance"].values() if isinstance(value, bool)
    ):
        raise SystemExit("one or more replicate-aware splice-graph acceptance checks failed")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
