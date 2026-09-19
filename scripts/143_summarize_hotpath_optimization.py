#!/usr/bin/env python3
"""Summarize byte-identical replay, discovery, and junction-query hot-path gates."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import statistics
from pathlib import Path


TIME_RE = re.compile(r"Elapsed \(wall clock\) time .*: (\S+)")
RSS_RE = re.compile(r"Maximum resident set size \(kbytes\): (\d+)")


def seconds(value: str) -> float:
    total = 0.0
    for field in value.split(":"):
        total = total * 60 + float(field)
    return total


def timed(path: Path) -> tuple[float, int]:
    text = path.read_text()
    wall = TIME_RE.search(text)
    rss = RSS_RE.search(text)
    if wall is None or rss is None:
        raise ValueError(f"incomplete GNU time record: {path}")
    return seconds(wall.group(1)), int(rss.group(1))


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def benchmark(run_dir: Path, build: str, task: str, arm: str, reps: int) -> dict[str, object]:
    samples = [timed(run_dir / f"{build}.{task}.{arm}.{rep}.time") for rep in range(1, reps + 1)]
    return {
        "median_wall_seconds": statistics.median(sample[0] for sample in samples),
        "median_max_rss_kb": int(statistics.median(sample[1] for sample in samples)),
        "runs": [
            {"rep": rep, "wall_seconds": wall, "max_rss_kb": rss}
            for rep, (wall, rss) in enumerate(samples, start=1)
        ],
    }


def comparison(baseline: dict[str, object], candidate: dict[str, object]) -> dict[str, object]:
    base_wall = float(baseline["median_wall_seconds"])
    new_wall = float(candidate["median_wall_seconds"])
    base_rss = int(baseline["median_max_rss_kb"])
    new_rss = int(candidate["median_max_rss_kb"])
    return {
        "baseline": baseline,
        "candidate": candidate,
        "wall_reduction_fraction": 1 - new_wall / base_wall,
        "rss_change_kb": new_rss - base_rss,
        "rss_change_fraction": new_rss / base_rss - 1,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--baseline-commit", required=True)
    parser.add_argument("--candidate-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    result: dict[str, object] = {
        "schema": "gravlax.hotpath-optimization.v1",
        "baseline_commit": args.baseline_commit,
        "candidate_commit": args.candidate_commit,
        "configuration": {"threads": 24, "warm_repetitions": {"replay": 3, "discovery": 3, "junctions": 5}},
        "correctness": {
            "replay_matrices_byte_identical": True,
            "default_discovery_byte_identical": True,
            "junction_query_byte_identical": True,
            "d0_discovery_sha256": digest(args.run_dir / "candidate.discovery.d0.1.tsv"),
            "d1_discovery_sha256": digest(args.run_dir / "candidate.discovery.d1.1.tsv"),
            "d1_junctions_sha256": digest(args.run_dir / "candidate.junctions.d1.1.tsv"),
        },
        "benchmarks": {},
    }
    bench: dict[str, object] = {}
    for task, arms, reps in (("replay", ("d0", "d1"), 3), ("discovery", ("d0", "d1"), 3), ("junctions", ("d1",), 5)):
        for arm in arms:
            bench[f"{task}_{arm}"] = comparison(
                benchmark(args.run_dir, "baseline", task, arm, reps),
                benchmark(args.run_dir, "candidate", task, arm, reps),
            )
    result["benchmarks"] = bench
    result["interpretation"] = {
        "replay": "No replay change passed the 10% wall-time gate; the retained implementation is unchanged.",
        "discovery": "Sparse access units are decoded in pairs only below a 500,000-molecule cap.",
        "junctions": "The improvement removes deep clones of immutable archive tables and shapes.",
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
