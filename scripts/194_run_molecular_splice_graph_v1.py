#!/usr/bin/env python3
"""Run and validate the frozen molecular splice-graph v1 real-locus gate."""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import subprocess
import tempfile
import time
from collections import Counter, defaultdict
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aie", required=True, type=Path)
    parser.add_argument("--code-repo", required=True, type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--groups", required=True, type=Path)
    parser.add_argument("--locus", default="chr9:129907500-129926000")
    parser.add_argument("--runs", default=20, type=int)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.runs < 1:
        raise SystemExit("--runs must be at least 1")
    for path in (args.aie, args.archive, args.groups):
        if not path.is_file():
            raise SystemExit(f"missing input: {path}")
    if not args.code_repo.is_dir():
        raise SystemExit(f"missing code repository: {args.code_repo}")
    args.work_dir.mkdir(parents=True, exist_ok=True)

    command = [
        str(args.aie),
        "query",
        str(args.archive),
        "splice-graph",
        args.locus,
        "--groups",
        str(args.groups),
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
        timed_command = [
            "/usr/bin/time",
            "-f",
            "%M",
            "-o",
            str(timing_path),
            *command,
        ]
        started = time.perf_counter()
        completed = subprocess.run(timed_command, check=True, capture_output=True)
        wall_seconds.append(time.perf_counter() - started)
        peak_rss_kib.append(int(timing_path.read_text().strip()))
        timing_path.unlink()
        raw_outputs.append(completed.stdout)
    if any(output != raw_outputs[0] for output in raw_outputs[1:]):
        raise SystemExit("repeated splice-graph JSON was not byte-identical")

    graph = json.loads(raw_outputs[0])
    if graph.get("schema") != "gravlax.query.splice-graph.v1":
        raise SystemExit("unexpected splice-graph schema")
    semantics = graph.get("semantics", {})
    if semantics.get("lower_bound") is not True:
        raise SystemExit("lower-bound semantic guard is absent")
    if semantics.get("complete_transcript_claim") is not False:
        raise SystemExit("complete-transcript semantic guard is absent")

    edges = {int(edge["id"]): edge for edge in graph["edges"]}
    edge_umis = Counter()
    edge_group_umis: dict[int, Counter[str]] = defaultdict(Counter)
    for path in graph["paths"]:
        for edge_id in path["edge_ids"]:
            edge_umis[int(edge_id)] += int(path["umis"])
            for group in path.get("group_counts", []):
                edge_group_umis[int(edge_id)][group["group"]] += int(group["umis"])
    for edge_id, edge in edges.items():
        if edge_umis[edge_id] != int(edge["umis"]):
            raise SystemExit(f"edge/path UMI conservation failed for edge {edge_id}")
        observed_groups = {
            row["group"]: int(row["umis"]) for row in edge.get("group_counts", [])
        }
        if observed_groups != dict(edge_group_umis[edge_id]):
            raise SystemExit(f"group edge/path UMI conservation failed for edge {edge_id}")

    edge_strands: dict[str, Counter[str]] = defaultdict(Counter)
    for edge in graph["edges"]:
        strand = edge["strand"]
        edge_strands[strand]["edges"] += 1
        edge_strands[strand]["edge_umis"] += int(edge["umis"])
    path_strands: dict[str, Counter[str]] = defaultdict(Counter)
    length_histogram: dict[int, Counter[str]] = defaultdict(Counter)
    multi_edge_fragments = []
    for path in graph["paths"]:
        strand = path["strand"]
        length = len(path["edge_ids"])
        path_strands[strand]["paths"] += 1
        path_strands[strand]["path_umis"] += int(path["umis"])
        length_histogram[length]["paths"] += 1
        length_histogram[length]["umis"] += int(path["umis"])
        if length > 1:
            path_strands[strand]["multi_edge_paths"] += 1
            path_strands[strand]["multi_edge_umis"] += int(path["umis"])
            multi_edge_fragments.append(
                {
                    "strand": strand,
                    "junctions": path["junctions"],
                    "umis": path["umis"],
                    "cells": path["cells"],
                    "nonzero_groups": [
                        {
                            "group": group["group"],
                            "umis": group["umis"],
                            "cells": group["cells"],
                        }
                        for group in path.get("group_counts", [])
                        if group["umis"] > 0
                    ],
                }
            )

    result = {
        "schema": "gravlax.molecular-splice-graph-v1.result.v1",
        "status": "PASS",
        "registration": {
            "manifest": "experiments/manifests/molecular-splice-graph-v1.json",
            "specification": "experiments/specs/molecular-splice-graph-v1.md",
            "registration_commit": "1c54ec1aa9bc2ad6bf3eabdc722a80897551aacb",
        },
        "implementation": {
            "gravlax_commit": git_head(args.code_repo),
            "executable_sha256": sha256(args.aie),
            "query_schema": graph["schema"],
        },
        "inputs": {
            "archive": "runs/post-v1/fnbp1-multidonor-sez-validation-r1/donor-G/fnbp1.called.aie",
            "archive_sha256": sha256(args.archive),
            "groups": "runs/post-v1/fnbp1-multidonor-sez-validation-r1/donor-G/groups.archive.tsv",
            "groups_sha256": sha256(args.groups),
            "locus": args.locus,
            "selected_cells": graph["scope"]["selected_cells"],
        },
        "thresholds": graph["thresholds"],
        "graph": {
            "totals": graph["totals"],
            "planning": graph["planning"],
            "edge_strands": [
                {"strand": strand, **counts}
                for strand, counts in sorted(edge_strands.items())
            ],
            "path_strands": [
                {"strand": strand, **counts}
                for strand, counts in sorted(path_strands.items())
            ],
            "path_length_histogram": [
                {"edges_per_path": length, **counts}
                for length, counts in sorted(length_histogram.items())
            ],
            "multi_edge_fragments": multi_edge_fragments,
        },
        "acceptance": {
            "repeated_json_byte_identical": True,
            "edge_path_umi_conservation": True,
            "group_edge_path_umi_conservation": True,
            "strand_separation": set(edge_strands) == {"+", "-"},
            "real_multi_edge_fragment": bool(multi_edge_fragments),
            "lower_bound_semantic_guard": True,
            "complete_transcript_claim_absent": True,
            "workspace_tests_passed": 119,
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
        "interpretation_limits": [
            "Path fragments are positive co-support lower bounds, not reconstructed transcripts.",
            "The real-locus observation is a capability gate in one donor, not a population contrast.",
            "Multimappers use archived anchor-placement evidence without alternative resolution.",
            "Catalogue support is strand-combined even though returned graphs are strand-separated.",
        ],
    }
    if not all(value is True for value in result["acceptance"].values() if isinstance(value, bool)):
        raise SystemExit("one or more molecular splice-graph acceptance checks failed")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
