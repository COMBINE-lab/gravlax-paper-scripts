#!/usr/bin/env python3
"""Aggregate the prospectively frozen eight-donor sparse-tail pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path


DONORS = tuple("ABCDEFGH")


def weighted_auc(rows: list[tuple[bool, float, int]]) -> float | None:
    positives = sum(count for label, _, count in rows if label)
    negatives = sum(count for label, _, count in rows if not label)
    if not positives or not negatives:
        return None
    rows.sort(key=lambda row: row[1])
    positive_rank_sum = 0.0
    rank = 1
    start = 0
    while start < len(rows):
        end = start + 1
        while end < len(rows) and rows[end][1] == rows[start][1]:
            end += 1
        tied = sum(count for _, _, count in rows[start:end])
        positive_tied = sum(count for label, _, count in rows[start:end] if label)
        average_rank = (rank + rank + tied - 1) / 2.0
        positive_rank_sum += average_rank * positive_tied
        rank += tied
        start = end
    return (positive_rank_sum - positives * (positives + 1) / 2) / (positives * negatives)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def aucs(rows: list[dict[str, int | str]]) -> tuple[float, float]:
    clip_rows = []
    tail_rows = []
    for row in rows:
        label = row["class"] == "external"
        count = int(row["witnesses"])
        clip_len = int(row["clip_len"])
        tail_bases = int(row["tail_bases"])
        terminal_run = int(row["terminal_run"])
        clip_rows.append((label, float(clip_len), count))
        score = tail_bases / max(clip_len, 1) + terminal_run / 1_000_000_000
        tail_rows.append((label, score, count))
    clip_auc = weighted_auc(clip_rows)
    tail_auc = weighted_auc(tail_rows)
    if clip_auc is None or tail_auc is None:
        raise ValueError("both candidate classes are required for AUC")
    return clip_auc, tail_auc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workspace-tests-passed", action="store_true")
    parser.add_argument("--workspace-test-count", type=int, default=0)
    parser.add_argument("--aie-binary", type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit(f"refusing to overwrite {args.out}")

    audits = []
    hist: dict[tuple[bool, int, int, int], int] = defaultdict(int)
    totals: dict[str, int] = defaultdict(int)
    catalogue = defaultdict(int)
    control = defaultdict(int)
    internal = defaultdict(int)
    non_internal = defaultdict(int)
    storage = defaultdict(int)
    input_digests = {}
    per_donor = []

    for donor in DONORS:
        path = args.run_root / f"donor-{donor}" / "audit.json"
        audit = json.loads(path.read_text())
        if audit["schema"] != "gravlax.sparse-terminal-tail-audit.v1":
            raise SystemExit(f"unexpected schema in {path}")
        audits.append(audit)
        input_digests[str(path)] = sha256(path)
        for key, value in audit["counts"].items():
            if isinstance(value, int):
                totals[key] += value
        for target, key in (
            (catalogue, "external_site_exclusive_events"),
            (control, "shifted_control_exclusive_events"),
            (internal, "internal_priming_site_events"),
            (non_internal, "non_internal_priming_site_events"),
        ):
            for field, value in audit["counts"][key].items():
                target[field] += value
        for key in ("sidecar_bytes", "archive_bytes", "uncompressed_record_payload_bytes"):
            storage[key] += audit["storage"][key]
        for row in audit["candidate_discrimination"]["event_histogram"]:
            key = (
                row["class"] == "external",
                row["clip_len"],
                row["tail_bases"],
                row["terminal_run"],
            )
            hist[key] += row["witnesses"]
        donor_clip_auc, donor_tail_auc = aucs(
            audit["candidate_discrimination"]["event_histogram"]
        )
        per_donor.append(
            {
                "sample": donor,
                "tail_positive_molecules": audit["counts"]["tail_positive_molecules"],
                "external_site_witnesses": audit["counts"][
                    "external_site_tail_positive_molecules"
                ],
                "external_vs_control_odds_ratio": audit["candidate_discrimination"][
                    "tail_positive_odds_ratio_external_vs_control"
                ],
                "tail_score_auc_gain": donor_tail_auc - donor_clip_auc,
                "sidecar_fraction_of_archive": audit["storage"][
                    "sidecar_fraction_of_archive"
                ],
                "elapsed_seconds": audit["elapsed_seconds"],
            }
        )

    pooled_rows = []
    for (label, clip_len, tail_bases, terminal_run), count in hist.items():
        pooled_rows.append(
            {
                "class": "external" if label else "shifted_control",
                "clip_len": clip_len,
                "tail_bases": tail_bases,
                "terminal_run": terminal_run,
                "witnesses": count,
            }
        )
    clip_auc, tail_auc = aucs(pooled_rows)
    odds_ratio = ((catalogue["positive"] + 0.5) * (control["clipped"] - control["positive"] + 0.5)) / (
        (catalogue["clipped"] - catalogue["positive"] + 0.5) * (control["positive"] + 0.5)
    )
    internal_rate = internal["positive"] / max(internal["clipped"], 1)
    non_internal_rate = non_internal["positive"] / max(non_internal["clipped"], 1)
    storage_fraction = storage["sidecar_bytes"] / max(storage["archive_bytes"], 1)
    gates = {
        "implementation": {
            "orientation_packing_and_dedup_tests_passed": args.workspace_tests_passed,
            "workspace_test_count": args.workspace_test_count,
            "pass": args.workspace_tests_passed and args.workspace_test_count > 0,
        },
        "external_enrichment": {
            "threshold_odds_ratio": 5.0,
            "threshold_witnesses": 100,
            "observed_odds_ratio": odds_ratio,
            "observed_witnesses": totals["external_site_tail_positive_molecules"],
            "pass": odds_ratio >= 5 and totals["external_site_tail_positive_molecules"] >= 100,
        },
        "candidate_discrimination": {
            "threshold_auc_gain": 0.05,
            "clip_length_auc": clip_auc,
            "continuous_tail_score_auc": tail_auc,
            "observed_auc_gain": tail_auc - clip_auc,
            "pass": tail_auc - clip_auc >= 0.05,
        },
        "internal_priming": {
            "tail_rate_internal": internal_rate,
            "tail_rate_non_internal": non_internal_rate,
            "pass": internal_rate < non_internal_rate,
        },
        "storage": {
            "threshold_fraction_of_archive": 0.01,
            "observed_fraction_of_archive": storage_fraction,
            "pass": storage_fraction <= 0.01,
        },
    }
    result = {
        "schema": "gravlax.sparse-terminal-tail-cohort-pilot.v1",
        "status": "PASS_PRIMARY_GATES"
        if all(gate["pass"] for gate in gates.values())
        else "FAIL_PRIMARY_GATES",
        "frozen_spec": "docs/pilots/2026-09-01-sparse-terminal-tail.md",
        "samples": list(DONORS),
        "counts": dict(totals),
        "exclusive_event_counts": {
            "external": dict(catalogue),
            "shifted_control": dict(control),
            "internal_priming": dict(internal),
            "non_internal_priming": dict(non_internal),
        },
        "storage": {
            **dict(storage),
            "sidecar_fraction_of_archive": storage_fraction,
        },
        "gates": gates,
        "per_donor": per_donor,
        "input_sha256": input_digests,
        "execution": {
            "aie_binary_sha256": sha256(args.aie_binary) if args.aie_binary else None,
            "workspace_test_count": args.workspace_test_count,
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
