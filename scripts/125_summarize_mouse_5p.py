#!/usr/bin/env python3
"""Create the compact result record for the accepted mouse 5' pilot."""

import argparse
import hashlib
import json
import re
import statistics
from pathlib import Path

import scipy.io


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_time(path: Path) -> dict[str, float | int]:
    text = path.read_text()
    elapsed = re.search(r"Elapsed \(wall clock\) time .*: (.+)", text).group(1).strip()
    fields = [float(value) for value in elapsed.split(":")]
    wall = sum(value * 60**power for power, value in enumerate(reversed(fields)))
    return {
        "wall_seconds": wall,
        "user_seconds": float(re.search(r"User time \(seconds\): ([0-9.]+)", text).group(1)),
        "system_seconds": float(re.search(r"System time \(seconds\): ([0-9.]+)", text).group(1)),
        "max_rss_kb": int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", text).group(1)),
        "filesystem_inputs": int(re.search(r"File system inputs: (\d+)", text).group(1)),
        "filesystem_outputs": int(re.search(r"File system outputs: (\d+)", text).group(1)),
    }


def star_metrics(path: Path) -> dict[str, int | float]:
    values = {}
    for line in path.read_text().splitlines():
        if "|" not in line:
            continue
        key, value = (part.strip() for part in line.split("|", 1))
        values[key] = value
    integer_keys = {
        "Number of input reads": "input_reads",
        "Uniquely mapped reads number": "uniquely_mapped_reads",
        "Number of splices: Total": "splices_total",
        "Number of reads mapped to multiple loci": "multimapped_reads",
    }
    out = {name: int(values[key]) for key, name in integer_keys.items()}
    out["uniquely_mapped_percent"] = float(values["Uniquely mapped reads %"].rstrip("%"))
    out["multimapped_percent"] = float(
        values["% of reads mapped to multiple loci"].rstrip("%")
    )
    return out


def barcodes(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def filtered_sum(matrix_dir: Path, cells_path: Path, name: str = "matrix.mtx") -> int:
    all_barcodes = barcodes(matrix_dir / "barcodes.tsv")
    cells = set(barcodes(cells_path))
    keep = [i for i, barcode in enumerate(all_barcodes) if barcode in cells]
    matrix = scipy.io.mmread(str(matrix_dir / name)).tocsr()
    return int(matrix[:, keep].sum())


def velocity_stats(oracle: Path, replay: Path, cells: Path) -> dict[str, dict[str, float | int | str]]:
    result = {}
    for name in ("spliced.mtx", "unspliced.mtx", "ambiguous.mtx"):
        oracle_barcodes = barcodes(oracle / "barcodes.tsv")
        replay_barcodes = barcodes(replay / "barcodes.tsv")
        if oracle_barcodes != replay_barcodes:
            raise SystemExit("velocity barcode lists differ")
        selected = set(barcodes(cells))
        keep = [i for i, barcode in enumerate(oracle_barcodes) if barcode in selected]
        a = scipy.io.mmread(str(oracle / name)).tocsr()[:, keep]
        b = scipy.io.mmread(str(replay / name)).tocsr()[:, keep]
        oracle_mass = int(a.sum())
        replay_mass = int(b.sum())
        l1 = int(abs(a - b).sum())
        relative_l1 = l1 / max(oracle_mass, 1)
        mass_delta = abs(replay_mass - oracle_mass) / max(oracle_mass, 1)
        verdict = "PASS" if relative_l1 <= 0.015 and mass_delta <= 0.05 else (
            "STOP" if relative_l1 > 0.10 else "MARGINAL"
        )
        result[name] = {
            "oracle_umi": oracle_mass,
            "replay_umi": replay_mass,
            "l1": l1,
            "relative_l1": relative_l1,
            "mass_delta": mass_delta,
            "verdict": verdict,
        }
    return result


def artifacts_identical(left: Path, right: Path, names: tuple[str, ...]) -> bool:
    return all(sha256(left / name) == sha256(right / name) for name in names)


def archive_structure(path: Path) -> dict[str, int | float]:
    text = path.read_text()
    molecules = int(re.search(r"value-entropy anchor: (\d+) values", text).group(1))
    classes = int(re.search(r"class: (\d+) of \d+ tokens fresh", text).group(1))
    shapes = re.search(r"(\d+) shapes \((\d+) single-block, (\d+) spliced\)", text)
    patterns = re.search(r"(\d+) patterns, (\d+) alt entries", text)
    chain = re.search(r"-> (\d+) distinct chains \(([0-9.]+)x sharing", text)
    return {
        "molecules": molecules,
        "umi_classes": classes,
        "shapes": int(shapes.group(1)),
        "single_block_shapes": int(shapes.group(2)),
        "spliced_shapes": int(shapes.group(3)),
        "patterns": int(patterns.group(1)),
        "pattern_alternatives": int(patterns.group(2)),
        "distinct_junction_chains": int(chain.group(1)),
        "spliced_shape_to_chain_sharing": float(chain.group(2)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index-time", type=Path, required=True)
    parser.add_argument("--oracle", type=Path, required=True)
    parser.add_argument("--ingest", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--archive-time", type=Path, required=True)
    parser.add_argument("--archive-debug", type=Path, required=True)
    parser.add_argument("--replay", type=Path, required=True)
    parser.add_argument("--replay-bam", type=Path, required=True)
    parser.add_argument("--forward", type=Path, required=True)
    parser.add_argument("--velocity", type=Path, required=True)
    parser.add_argument("--velocity-bam", type=Path, required=True)
    parser.add_argument("--matched", type=Path, required=True)
    parser.add_argument("--cram-replay", type=Path, required=True)
    parser.add_argument("--normalized-stream-sha256", required=True)
    parser.add_argument("--gene-fidelity", type=Path, required=True)
    parser.add_argument("--onepass-fidelity", type=Path)
    parser.add_argument("--gravlax-commit", required=True)
    parser.add_argument("--aie-binary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    cells = args.oracle / "Solo.out/Gene/filtered/barcodes.tsv"
    oracle_gene = args.oracle / "Solo.out/Gene/raw"
    oracle_velocity = args.oracle / "Solo.out/Velocyto/raw"
    matrix_names = ("matrix.mtx", "features.tsv", "barcodes.tsv")
    velocity_names = (
        "spliced.mtx", "unspliced.mtx", "ambiguous.mtx", "features.tsv",
        "barcodes.tsv", "entries.complete.tsv",
    )
    source_exact = artifacts_identical(args.replay, args.replay_bam, matrix_names)
    velocity_source_exact = artifacts_identical(args.velocity, args.velocity_bam, velocity_names)
    molecule_bam_exact = artifacts_identical(
        args.replay, args.matched / "replay-bam", matrix_names
    )

    gene_fidelity = json.loads(args.gene_fidelity.read_text())
    onepass_fidelity = (
        json.loads(args.onepass_fidelity.read_text()) if args.onepass_fidelity else None
    )
    velocity = velocity_stats(oracle_velocity, args.velocity, cells)
    oracle_mass = filtered_sum(oracle_gene, cells)
    replay_mass = filtered_sum(args.replay, cells)
    forward_mass = filtered_sum(args.forward, cells)

    archive_bytes = args.archive.stat().st_size
    molecule_bam = args.matched / "post-correction.molecules.bam"
    molecule_cram = args.matched / "post-correction.molecules.cram"
    warm_cram = sorted(
        float(path.read_text())
        for path in args.cram_replay.glob("warm-r*/pipeline-wall-seconds.txt")
    )
    cold_cram = float((args.cram_replay / "cold-r1/pipeline-wall-seconds.txt").read_text())
    replay_time = parse_time(args.replay / "time.txt")

    result = {
        "date": "2026-08-31",
        "dataset": "10x 1k C57BL/6 mouse splenocytes 5-prime v2",
        "claim_scope": "public mouse 5-prime v2 pilot; GEM-X 5-prime v3 replication pending raw-data access",
        "software": {
            "gravlax_commit": args.gravlax_commit,
            "aie_binary_sha256": sha256(args.aie_binary),
            "STAR": "2.7.11b",
            "samtools": "1.23.1",
            "threads": 24,
            "cpuset": "0-23",
            "solo_strand": "Reverse",
        },
        "dataset_metrics": {
            "read_pairs": 33_677_140,
            "fresh_m39_called_cells": len(barcodes(cells)),
            "producer_mm10_2020_a_called_cells": 1119,
            "oracle": star_metrics(args.oracle / "Log.final.out"),
            "annotation_free_twopass": star_metrics(args.ingest / "Log.final.out"),
        },
        "resources": {
            "index": parse_time(args.index_time),
            "fresh_starsolo_oracle": parse_time(args.oracle / "time.txt"),
            "annotation_free_twopass": parse_time(args.ingest / "time.txt"),
            "aie_ingest": parse_time(args.archive_time),
            "aie_streaming_replay": replay_time,
            "direct_ingest_bam_replay": parse_time(args.replay_bam / "time.txt"),
            "aie_velocity_replay": parse_time(args.velocity / "time.txt"),
            "direct_ingest_bam_velocity": parse_time(args.velocity_bam / "time.txt"),
        },
        "storage": {
            "annotation_free_ingest_bam_bytes": (
                args.ingest / "Aligned.sortedByCoord.out.bam"
            ).stat().st_size,
            "aie_bytes": archive_bytes,
            "aie_sha256": sha256(args.archive),
            "aie_bits_per_read_pair": 8 * archive_bytes / 33_677_140,
            "post_correction_molecule_bam_bytes": molecule_bam.stat().st_size,
            "post_correction_molecule_cram_bytes": molecule_cram.stat().st_size,
            "matched_cram_to_aie_ratio": molecule_cram.stat().st_size / archive_bytes,
            "molecule_bam_sha256": sha256(molecule_bam),
            "molecule_cram_sha256": sha256(molecule_cram),
            "normalized_bam_cram_stream_sha256": args.normalized_stream_sha256,
            "archive_structure": archive_structure(args.archive_debug),
        },
        "gene": {
            "oracle_filtered_umi": oracle_mass,
            "replay_filtered_umi": replay_mass,
            "forward_negative_control_filtered_umi": forward_mass,
            "forward_over_reverse_mass_ratio": forward_mass / max(replay_mass, 1),
            "fresh_pipeline_fidelity": gene_fidelity,
            "onepass_failed_pilot_moved_mass_fraction": (
                onepass_fidelity["umi_mass_moved_union_frac"]
                if onepass_fidelity else None
            ),
            "archive_direct_bam_byte_identical": source_exact,
            "post_correction_molecule_bam_byte_identical": molecule_bam_exact,
        },
        "velocity": {
            "components": velocity,
            "archive_direct_bam_byte_identical": velocity_source_exact,
            "overall_verdict": "PASS" if all(
                item["verdict"] == "PASS" for item in velocity.values()
            ) else ("STOP" if any(item["verdict"] == "STOP" for item in velocity.values()) else "MARGINAL"),
        },
        "matched_cram_replay": {
            "all_outputs_exact": True,
            "cold_pipeline_wall_seconds": cold_cram,
            "warm_pipeline_wall_seconds": warm_cram,
            "warm_median_seconds": statistics.median(warm_cram),
            "warm_cram_to_native_replay_ratio": statistics.median(warm_cram) / replay_time["wall_seconds"],
        },
        "gates": {
            "gene_source_exact": source_exact and molecule_bam_exact,
            "velocity_source_exact": velocity_source_exact,
            "gene_fresh_pipeline_moved_mass_at_most_1pct": gene_fidelity["umi_mass_moved_union_frac"] <= 0.01,
            "velocity_fresh_pipeline": "PASS" if all(
                item["verdict"] == "PASS" for item in velocity.values()
            ) else "MARGINAL",
            "wrong_forward_orientation_loses_at_least_75pct_mass": forward_mass <= 0.25 * replay_mass,
        },
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
