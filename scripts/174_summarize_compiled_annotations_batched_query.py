#!/usr/bin/env python3
"""Apply the preregistered compiled-annotation and batch-query gates."""
from __future__ import annotations

import argparse
import hashlib
import json
import statistics
from pathlib import Path


def timings(path: Path) -> tuple[list[float], list[int]]:
    walls, rss = [], []
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        walls.append(float(fields[-2]))
        rss.append(int(fields[-1]))
    return walls, rss


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = args.root

    aic = root / "annotation/gencode.v49.annotation.aic"
    gtf_walls, gtf_rss = timings(root / "annotation-query/gtf-times.tsv")
    aic_walls, aic_rss = timings(root / "annotation-query/aic-times.tsv")
    batch_walls, batch_rss = timings(root / "batch/times.tsv")
    independent_walls, _ = timings(root / "independent/times.tsv")
    gtf_median = statistics.median(gtf_walls)
    aic_median = statistics.median(aic_walls)
    batch_median = statistics.median(batch_walls)
    independent_median = statistics.median(independent_walls)
    planning = json.loads((root / "batch/run-1.json").read_text())["planning"]
    correctness = json.loads((root / "batch/correctness.json").read_text())
    aic_fraction = aic.stat().st_size / args.gtf.stat().st_size
    query_fraction = aic_median / gtf_median
    batch_fraction = batch_median / independent_median

    acceptance = {
        "compiled_builds_byte_identical": (
            digest(aic) == digest(root / "annotation/gencode.v49.annotation.repeat.aic")
        ),
        "compiled_artifact_at_most_5pct_of_gtf": aic_fraction <= 0.05,
        "compiled_query_at_most_25pct_of_gtf_wall": query_fraction <= 0.25,
        "compiled_query_peak_rss_at_most_1gib": max(aic_rss) <= 1048576,
        "gtf_aic_query_byte_identical": (
            digest(root / "annotation-query/gtf-1.json")
            == digest(root / "annotation-query/aic-1.json")
        ),
        "gtf_aic_replay_byte_identical": all(
            digest(root / f"replay-gtf/output/{name}")
            == digest(root / f"replay-aic/output/{name}")
            for name in ("matrix.mtx", "features.tsv", "barcodes.tsv")
        ),
        "gtf_aic_discovery_byte_identical": (
            digest(root / "discovery/gtf.tsv") == digest(root / "discovery/aic.tsv")
        ),
        "all_32_batch_answers_exact": correctness["queries_checked"] == 32,
        "batch_chunk_decode_reduction_at_least_75pct": (
            planning["chunk_decode_reduction_fraction"] >= 0.75
        ),
        "batch_wall_at_most_35pct_of_independent": batch_fraction <= 0.35,
        "batch_peak_rss_at_most_1gib": max(batch_rss) <= 1048576,
    }
    if not all(acceptance.values()):
        failed = [name for name, passed in acceptance.items() if not passed]
        raise SystemExit(f"gate failure: {', '.join(failed)}")

    compile_walls, compile_rss = timings(root / "annotation/compile-time.tsv")
    replay_gtf_walls, replay_gtf_rss = timings(root / "replay-gtf/time.tsv")
    replay_aic_walls, replay_aic_rss = timings(root / "replay-aic/time.tsv")
    result = {
        "schema": "gravlax.compiled-annotations-batched-query.v1",
        "status": "PASS",
        "gravlax_commit": args.commit,
        "configuration": {"dataset": "D0", "threads": 24},
        "compiled_annotation": {
            "source_gtf_bytes": args.gtf.stat().st_size,
            "artifact_bytes": aic.stat().st_size,
            "artifact_fraction_of_gtf": aic_fraction,
            "artifact_sha256": digest(aic),
            "compile_wall_seconds": compile_walls[0],
            "compile_peak_rss_kib": compile_rss[0],
            "annotation_query_gtf_median_wall_seconds": gtf_median,
            "annotation_query_aic_median_wall_seconds": aic_median,
            "annotation_query_wall_fraction": query_fraction,
            "annotation_query_aic_max_rss_kib": max(aic_rss),
            "replay_gtf_wall_seconds": replay_gtf_walls[0],
            "replay_aic_wall_seconds": replay_aic_walls[0],
            "replay_gtf_peak_rss_kib": replay_gtf_rss[0],
            "replay_aic_peak_rss_kib": replay_aic_rss[0],
        },
        "batched_query": {
            "queries": correctness["queries_checked"],
            "independent_chunk_decodes": planning["independent_chunk_decodes"],
            "unique_chunk_decodes": planning["unique_chunk_decodes"],
            "chunk_decode_reduction_fraction": planning["chunk_decode_reduction_fraction"],
            "batch_median_wall_seconds": batch_median,
            "independent_median_wall_seconds": independent_median,
            "batch_wall_fraction": batch_fraction,
            "batch_max_rss_kib": max(batch_rss),
            "batch_output_sha256": digest(root / "batch/run-1.json"),
        },
        "acceptance": acceptance,
        "decision": {
            "compiled_annotation": "publish_aic_as_public_v1_artifact",
            "batched_query": "publish_region_and_exact_junction_plan_v1",
        },
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
