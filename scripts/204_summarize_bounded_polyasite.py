#!/usr/bin/env python3
"""Fail-closed summary for the bounded-memory PolyASite hot path."""

import argparse
import hashlib
import json
import re
from pathlib import Path


TABLES = ("sites.tsv", "genes.tsv", "fragment-kernel.tsv")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_time(path: Path) -> dict:
    text = path.read_text()
    wall_match = re.search(r"Elapsed \(wall clock\) time.*: ([0-9:.]+)", text)
    rss_match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    status_match = re.search(r"Exit status: (\d+)", text)
    if not wall_match or not rss_match or not status_match:
        raise ValueError(f"incomplete GNU time record: {path}")
    fields = [float(value) for value in wall_match.group(1).split(":")]
    wall = fields[-1]
    if len(fields) >= 2:
        wall += 60 * fields[-2]
    if len(fields) == 3:
        wall += 3600 * fields[-3]
    return {
        "wall_seconds": wall,
        "max_rss_kib": int(rss_match.group(1)),
        "exit_status": int(status_match.group(1)),
        "record_sha256": sha256(path),
    }


def adjacent_time(directory: Path) -> Path:
    return Path(str(directory) + ".time.txt")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frozen", type=Path, required=True)
    parser.add_argument("--prior-optimized", type=Path, required=True)
    parser.add_argument("--bounded", type=Path, required=True)
    parser.add_argument("--gravlax-commit", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.gravlax_commit):
        raise SystemExit("--gravlax-commit must be a full lowercase 40-hex commit")

    frozen_digests = {table: sha256(args.frozen / table) for table in TABLES}
    prior_digests = {table: sha256(args.prior_optimized / table) for table in TABLES}
    bounded_digests = {table: sha256(args.bounded / table) for table in TABLES}
    exact = frozen_digests == prior_digests == bounded_digests
    prior = read_time(adjacent_time(args.prior_optimized))
    bounded = read_time(adjacent_time(args.bounded))
    rss_gate_kib = int(2.25 * 1024 * 1024)
    gates = {
        "scientific_tables_byte_identical": exact,
        "exit_status_zero": bounded["exit_status"] == 0,
        "wall_seconds_le_30": bounded["wall_seconds"] <= 30.0,
        "max_rss_kib_le_2_25_gib": bounded["max_rss_kib"] <= rss_gate_kib,
    }
    if not all(gates.values()):
        raise SystemExit(f"bounded PolyASite gate failed: {gates}")

    result = {
        "schema": "gravlax.polyasite-bounded-memory.v1",
        "status": "PASS",
        "gates": gates,
        "thresholds": {"wall_seconds": 30.0, "max_rss_kib": rss_gate_kib},
        "prior_optimized": prior,
        "bounded": bounded,
        "comparison": {
            "rss_reduction_fraction": 1.0
            - bounded["max_rss_kib"] / prior["max_rss_kib"],
            "wall_ratio": bounded["wall_seconds"] / prior["wall_seconds"],
        },
        "implementation": {
            "sample_frontier": "at most two adjacent archives and 256 MiB aggregate compressed archive bytes",
            "chunk_frontier": 8,
            "endpoint_frontier": "one chromosome per active archive",
            "counts": "flat sorted packed coordinate rows; downstream gene ranges are borrowed",
            "classification": "read-only prefetched COC blocks with parallel chunk classification",
        },
        "input_digests": {
            "gravlax_commit": args.gravlax_commit,
            "tables": frozen_digests,
            "frozen_summary_sha256": sha256(args.frozen / "summary.json"),
            "bounded_summary_sha256": sha256(args.bounded / "summary.json"),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
