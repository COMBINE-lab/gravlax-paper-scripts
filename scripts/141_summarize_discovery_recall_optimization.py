#!/usr/bin/env python3
"""Compact the locked discovery-mode sweep and enforce its selection guards."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


LABELS = (
    "span",
    "strand-span",
    "compatible",
    "residual-sites-25",
    "residual-sites-50",
    "residual-sites-75",
    "residual-sites-100",
)
TIME_RE = re.compile(r"Elapsed \(wall clock\) time .*: (\S+)")
RSS_RE = re.compile(r"Maximum resident set size \(kbytes\): (\d+)")


def elapsed_seconds(value: str) -> float:
    fields = value.split(":")
    total = 0.0
    for field in fields:
        total = total * 60 + float(field)
    return total


def timing(path: Path) -> dict[str, float | int]:
    text = path.read_text()
    wall = TIME_RE.search(text)
    rss = RSS_RE.search(text)
    if wall is None or rss is None:
        raise ValueError(f"incomplete GNU time record: {path}")
    return {"wall_seconds": elapsed_seconds(wall.group(1)), "max_rss_kb": int(rss.group(1))}


def rows(path: Path) -> list[tuple[str, int, int, str]]:
    result = []
    with path.open() as handle:
        for line in handle:
            chrom, start, end, strand, *_ = line.rstrip("\n").split("\t")
            result.append((chrom, int(start), int(end), strand))
    return result


def percentile(values: list[int], fraction: float) -> int:
    return values[round((len(values) - 1) * fraction)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    docs = {label: json.loads((args.run_dir / f"{label}.summary.json").read_text()) for label in LABELS}
    baseline_candidates = docs["span"]["accounting"]["candidate_rows"]
    volume_limit = baseline_candidates * 3
    result: dict[str, object] = {
        "schema": "gravlax.discovery-recall-optimization.v1",
        "protocol_lock_commits": [
            "c6edf01ba0612e5c081a11d84328d1c28f2c9965",
            "86193f5bdfb02789ee647e6b4e2843b4ee93c893",
            "4ca25df24eb222a5a93b5e9f626422b4a697eccb",
            "023eef4458495349bbf3cd12785dd887455b8c1a",
        ],
        "implementation_commit": docs["residual-sites-75"]["gravlax_code_commit"],
        "selection": {
            "candidate_volume_limit": volume_limit,
            "minimum_d1_recurrence": 0.95,
            "minimum_absolute_recall_gain": 0.15,
            "selected": "residual-sites-75",
        },
        "arms": {},
    }

    baseline_recall = docs["span"]["thresholds"]["50"]["macro_recall"]
    arms: dict[str, object] = {}
    for label, doc in docs.items():
        primary = doc["thresholds"]["50"]
        candidates = doc["accounting"]["candidate_rows"]
        recurrence = doc["independent_support"]["d0_candidate_d1_recurrence_fraction"]
        arms[label] = {
            "d0_candidates": candidates,
            "candidate_ratio_vs_span": candidates / baseline_candidates,
            "matched_truth_genes": primary["matched_genes"],
            "truth_genes": primary["truth_genes"],
            "recall": primary["macro_recall"],
            "absolute_recall_gain_vs_span": primary["macro_recall"] - baseline_recall,
            "umi_weighted_recall": primary["umi_weighted_recall"],
            "recall_by_stratum": {
                name: {
                    "matched": values["matched_genes"],
                    "total": values["truth_genes"],
                    "recall": values["recall"],
                }
                for name, values in primary["strata"].items()
            },
            "d1_candidate_recurrence": recurrence,
            "safeguards_pass": (
                candidates <= volume_limit
                and recurrence >= 0.95
                and primary["macro_recall"] - baseline_recall >= 0.15
            ),
            "timing": {
                arm: timing(args.run_dir / f"{label}.{arm}.time")
                if (args.run_dir / f"{label}.{arm}.time").exists()
                else None
                for arm in ("d0", "d1")
            },
        }
    result["arms"] = arms

    baseline = set(rows(args.run_dir / "span.d0.tsv"))
    selected = rows(args.run_dir / "residual-sites-75.d0.tsv")
    added_widths = sorted(end - start for chrom, start, end, strand in selected if (chrom, start, end, strand) not in baseline)
    result["selected_residual_site_geometry"] = {
        "added_candidates": len(added_widths),
        "median_width_bp": percentile(added_widths, 0.5),
        "p95_width_bp": percentile(added_widths, 0.95),
        "maximum_width_bp": max(added_widths),
        "all_at_most_merge_gap_plus_one": all(width <= 1001 for width in added_widths),
    }
    result["interpretation_guard"] = (
        "Later-annotation overlap is a recall endpoint, not biological precision; unsupported "
        "candidates remain unresolved. The selected residual threshold was chosen only from the "
        "predeclared grid and guards."
    )
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
