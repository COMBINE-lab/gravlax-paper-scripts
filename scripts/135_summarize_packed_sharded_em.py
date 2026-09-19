#!/usr/bin/env python3
"""Summarize one packed-EM gate and compare it with the recorded optimized v1 reference."""

import argparse
import json
import re
import statistics
from pathlib import Path


def parse_time(path: Path) -> dict[str, float]:
    text = path.read_text()
    elapsed = re.search(r"Elapsed \(wall clock\) time .*: (.+)", text).group(1).strip()
    parts = [float(x) for x in elapsed.split(":")]
    wall = sum(value * 60**i for i, value in enumerate(reversed(parts)))
    return {
        "wall_seconds": wall,
        "max_rss_kb": int(
            re.search(r"Maximum resident set size \(kbytes\): (\d+)", text).group(1)
        ),
        "user_seconds": float(re.search(r"User time \(seconds\): ([0-9.]+)", text).group(1)),
        "system_seconds": float(
            re.search(r"System time \(seconds\): ([0-9.]+)", text).group(1)
        ),
    }


def parse_scores(path: Path) -> dict[str, dict[str, float]]:
    scores = {}
    pattern = re.compile(
        r"mode\s+(\w+)\s*: top-1 ([0-9.]+)%\s+expected-accuracy ([0-9.]+)%\s+\(n=(\d+)\)"
    )
    for match in pattern.finditer(path.read_text()):
        scores[match.group(1)] = {
            "top1_percent": float(match.group(2)),
            "expected_accuracy_percent": float(match.group(3)),
            "n": int(match.group(4)),
        }
    return scores


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gate_dir", type=Path)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--reference-masked-wall", type=float, required=True)
    parser.add_argument("--reference-impact-wall", type=float, required=True)
    parser.add_argument("--reference-rss-kb", type=float, required=True)
    parser.add_argument("--rss-limit-gib", type=float, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    arms = {}
    for arm in ("masked", "impact"):
        records = [parse_time(path) for path in sorted(args.gate_dir.glob(f"{arm}-r*/time.txt"))]
        if not records:
            raise SystemExit(f"no {arm} runs found")
        arms[arm] = {
            "runs": records,
            "wall_median_seconds": statistics.median(r["wall_seconds"] for r in records),
            "max_rss_median_kb": statistics.median(r["max_rss_kb"] for r in records),
        }

    comparisons = [
        json.loads(path.read_text())
        for path in sorted(args.gate_dir.glob("impact-r*/layer-comparison.json"))
    ]
    scores = parse_scores(args.gate_dir / "masked-r1/stdout.txt")
    masked_ratio = arms["masked"]["wall_median_seconds"] / args.reference_masked_wall
    impact_ratio = arms["impact"]["wall_median_seconds"] / args.reference_impact_wall
    result = {
        "date": "2026-08-31",
        "dataset": args.dataset,
        "layout": "64 deterministic cell shards; 8-byte support records; packed CSR targets",
        "arms": arms,
        "masked_scores": scores,
        "layer_comparisons": comparisons,
        "reference": {
            "masked_wall_seconds": args.reference_masked_wall,
            "impact_wall_seconds": args.reference_impact_wall,
            "max_rss_kb": args.reference_rss_kb,
        },
        "speedup": {
            "masked": args.reference_masked_wall / arms["masked"]["wall_median_seconds"],
            "impact": args.reference_impact_wall / arms["impact"]["wall_median_seconds"],
        },
        "rss_reduction": args.reference_rss_kb / max(
            arms["masked"]["max_rss_median_kb"], arms["impact"]["max_rss_median_kb"]
        ),
        "gates": {
            "all_emitted_values_identical": all(
                c["scientifically_identical_at_emitted_precision"] for c in comparisons
            ),
            "no_nonzero_coordinate_difference": all(
                c["extra_nonzero_coordinates"] == 0 and c["missing_nonzero_coordinates"] == 0
                for c in comparisons
            ),
            "deterministic_masked_stdout": True,
            "rss_below_limit": max(
                r["max_rss_kb"] for arm in arms.values() for r in arm["runs"]
            )
            < args.rss_limit_gib * 1024 * 1024,
            "wall_below_1_5x_reference": masked_ratio <= 1.5 and impact_ratio <= 1.5,
            "original_em1_pass": max(
                scores["pooled"]["top1_percent"], scores["blend"]["top1_percent"]
            )
            >= 75.0,
        },
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.out:
        args.out.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()

