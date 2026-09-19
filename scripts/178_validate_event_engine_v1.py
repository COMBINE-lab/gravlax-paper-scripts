#!/usr/bin/env python3
"""Run and validate the frozen event-engine v1 gate.

Large query JSON stays under the excluded run directory. Only the compact summary is copied into
the tracked results tree after review.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import statistics
import subprocess
import time
from pathlib import Path


D0_SHA256 = "d7c11d4258f7b78b5c3a20ebc5349e3a1f1f52a0279425701d14acbf79948842"
D1_SHA256 = "e8568b72525a84342203e062fc0f544baf84fcb2b346da7d50c8e1ba8e78c834"
D4_SHA256 = "031330c0de6c18d0e93eacf76e2dd617db24dbab6b6d0dbb39c537dfefaf97e9"
AKR1A1_ID = (
    "cassette:chr1:45568177-45568484:45568684-45568926:45568177-45568926"
)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def execute(command: list[str], *, capture: bool = True) -> tuple[bytes, float]:
    started = time.perf_counter()
    result = subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE if capture else subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        env={**os.environ, "RAYON_NUM_THREADS": "24"},
    )
    return result.stdout if capture else b"", time.perf_counter() - started


def profile(command: list[str], timing_path: Path) -> tuple[bytes, float, int]:
    result = subprocess.run(
        ["/usr/bin/time", "-f", "%e\t%M", "-o", str(timing_path), *command],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        env={**os.environ, "RAYON_NUM_THREADS": "24"},
    )
    wall, rss = timing_path.read_text().strip().split("\t")
    return result.stdout, float(wall), int(rss)


def event_command(binary: Path, archive: Path, groups: Path) -> list[str]:
    return [
        str(binary),
        "query",
        str(archive),
        "events",
        "chr1:0-248956422",
        "--min-support",
        "2",
        "--min-informative",
        "10",
        "--groups",
        str(groups),
        "--agg",
        "group",
        "--json",
    ]


def jset_command(binary: Path, archive: Path, groups: Path, event: dict) -> list[str]:
    command = [str(binary), "query", str(archive), "jset"]
    for junction in event["inclusion_junctions"]:
        command.extend(("--include", junction["locus"]))
    for junction in event["exclusion_junctions"]:
        command.extend(("--exclude", junction["locus"]))
    command.extend(("--groups", str(groups), "--agg", "group", "--json"))
    return command


def legacy_bytes(binary: Path, root: Path) -> dict[str, bytes]:
    archive = root / "runs/archive/d0.aie"
    plan = root / "gravlax-paper-scripts/experiments/plans/batched-query-d0.tsv"
    commands = {
        "region": ["query", archive, "region", "chr1:45550000-45571000", "--json"],
        "junction": [
            "query", archive, "junction", "chr1:45568684-45568926", "--top", "0", "--json"
        ],
        "junctions": [
            "query", archive, "junctions", "chr1:45550000-45571000", "--with-cells", "--json"
        ],
        "jset": [
            "query", archive, "jset", "--include", "chr1:45568684-45568926",
            "--exclude", "chr1:45568177-45568926", "--json"
        ],
        "batch": ["query", archive, "batch", "--plan", plan, "--top", "20"],
        "federate": [
            "federate", archive, root / "runs/archive/d1.aie", "chr1:45568684-45568926"
        ],
    }
    outputs = {}
    for name, tail in commands.items():
        raw, _ = execute([str(binary), *(str(value) for value in tail)])
        if name == "federate":
            raw = re.sub(rb" \([0-9.]+s\)\n$", b" (TIME)\n", raw)
        outputs[name] = raw
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--bin", type=Path, required=True)
    parser.add_argument("--reference-bin", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=7)
    args = parser.parse_args()
    if args.repetitions < 1:
        raise SystemExit("--repetitions must be positive")
    args.out.mkdir(parents=True, exist_ok=True)

    paths = {
        "D0": args.root / "runs/archive/d0.aie",
        "D1": args.root / "runs/archive/d1.aie",
        "D4": args.root / "runs/archive/d4.aie",
    }
    expected = {"D0": D0_SHA256, "D1": D1_SHA256, "D4": D4_SHA256}
    observed = {name: digest(path) for name, path in paths.items()}
    assert observed == expected, (observed, expected)
    groups = {
        "D0": args.root / "runs/post-v1/hierarchical-em-groups-r1/groups/D0.real.tsv",
        "D1": args.root / "runs/post-v1/hierarchical-em-groups-r1/groups/D1.real.tsv",
        "D4": args.root / "gravlax-paper-scripts/results/post-v1-em-confirmation-d4/groups-real.tsv",
    }

    raw_d0, d0_wall, d0_rss = profile(
        event_command(args.bin, paths["D0"], groups["D0"]), args.out / "d0.time"
    )
    raw_d0_repeat, _ = execute(event_command(args.bin, paths["D0"], groups["D0"]))
    assert raw_d0 == raw_d0_repeat
    d0 = json.loads(raw_d0)
    raw_d1, d1_wall, d1_rss = profile(
        event_command(args.bin, paths["D1"], groups["D1"]), args.out / "d1.time"
    )
    d1 = json.loads(raw_d1)
    (args.out / "d0-chr1-events.json").write_bytes(raw_d0)
    (args.out / "d1-chr1-events.json").write_bytes(raw_d1)

    checked = min(50, len(d0["events"]))
    for event in d0["events"][:checked]:
        raw, _ = execute(jset_command(args.bin, paths["D0"], groups["D0"], event))
        standalone = json.loads(raw)
        assert event["totals"] == standalone["totals"], event["id"]
        assert event["group_rows"] == standalone["group_rows"], event["id"]

    cohort_command = [
        str(args.bin), "cohort", "events", "chr1:0-248956422",
        "--sample", f"D0={paths['D0']}", "--sample", f"D1={paths['D1']}",
        "--groups", f"D0={groups['D0']}", "--groups", f"D1={groups['D1']}",
        "--min-support", "2", "--min-samples", "2", "--min-informative", "10", "--json",
    ]
    raw_cohort, cohort_wall, cohort_rss = profile(cohort_command, args.out / "cohort.time")
    cohort = json.loads(raw_cohort)
    (args.out / "d01-chr1-cohort.json").write_bytes(raw_cohort)
    individual = {
        "D0": {event["id"]: event for event in d0["events"]},
        "D1": {event["id"]: event for event in d1["events"]},
    }
    cohort_rows_checked = 0
    for event in cohort["events"]:
        for row in event["sample_rows"]:
            source = individual[row["sample"]].get(event["id"])
            if source is None:
                assert row["totals"]["informative_umis"] < 10
                continue
            assert row["totals"] == source["totals"], event["id"]
            assert row.get("group_rows") == source.get("group_rows"), event["id"]
            cohort_rows_checked += 1

    d4_command = [
        str(args.bin), "query", str(paths["D4"]), "events", "chr1:45550000-45571000",
        "--event-type", "cassette", "--min-support", "2", "--min-informative", "1",
        "--groups", str(groups["D4"]), "--agg", "group", "--json",
    ]
    raw_d4, d4_wall, d4_rss = profile(d4_command, args.out / "d4.time")
    d4 = json.loads(raw_d4)
    (args.out / "d4-akr1-events.json").write_bytes(raw_d4)
    locked = next((event for event in d4["events"] if event["id"] == AKR1A1_ID), None)
    held = {
        "event_present": locked is not None,
        "informative_umis": 0 if locked is None else locked["totals"]["informative_umis"],
        "eligible_groups": [],
        "eligible_usage_range": None,
    }
    if locked is not None:
        held["eligible_groups"] = [
            {"group": row["group"], "informative_umis": row["informative_umis"],
             "usage_fraction": row["usage_fraction"]}
            for row in locked["group_rows"] if row["informative_umis"] >= 20
        ]
        usages = [row["usage_fraction"] for row in held["eligible_groups"]]
        if usages:
            held["eligible_usage_range"] = max(usages) - min(usages)
    held["pass"] = bool(
        held["event_present"]
        and held["informative_umis"] >= 100
        and len(held["eligible_groups"]) >= 2
        and held["eligible_usage_range"] is not None
        and held["eligible_usage_range"] >= 0.05
    )

    event_times = []
    standalone_times = []
    panel = d0["events"][:checked]
    for _ in range(args.repetitions):
        _, elapsed = execute(event_command(args.bin, paths["D0"], groups["D0"]), capture=False)
        event_times.append(elapsed)
        started = time.perf_counter()
        for event in panel:
            execute(jset_command(args.bin, paths["D0"], groups["D0"], event), capture=False)
        standalone_times.append(time.perf_counter() - started)
    event_median = statistics.median(event_times)
    standalone_median = statistics.median(standalone_times)

    legacy = {"checked": args.reference_bin is not None, "byte_identical": None, "sha256": {}}
    if args.reference_bin is not None:
        before = legacy_bytes(args.reference_bin, args.root)
        after = legacy_bytes(args.bin, args.root)
        legacy["byte_identical"] = before == after
        assert legacy["byte_identical"]
        legacy["sha256"] = {
            name: hashlib.sha256(payload).hexdigest() for name, payload in after.items()
        }

    result = {
        "schema": "gravlax.gate.event-engine-v1.v1",
        "status": "PASS_ENGINE_FAIL_HELD_OUT_BIOLOGY",
        "inputs_sha256": observed,
        "correctness": {
            "repeated_d0_json_byte_identical": True,
            "standalone_jset_events_checked": checked,
            "standalone_totals_and_group_rows_exact": True,
            "cohort_rows_checked": cohort_rows_checked,
            "cohort_equals_independent_queries": True,
            "legacy": legacy,
        },
        "D0_chr1": {
            **d0["planning"], "wall_seconds_first": d0_wall, "peak_rss_kib": d0_rss,
            "event_types": {
                kind: sum(event["event_type"] == kind for event in d0["events"])
                for kind in ("alt_acceptor", "alt_donor", "cassette")
            },
        },
        "D1_chr1": {
            **d1["planning"], "wall_seconds_first": d1_wall, "peak_rss_kib": d1_rss
        },
        "D0_D1_cohort": {
            **cohort["planning"], "wall_seconds_first": cohort_wall,
            "peak_rss_kib": cohort_rss,
        },
        "performance": {
            "repetitions": args.repetitions,
            "event_seconds": event_times,
            "standalone_50_jset_seconds": standalone_times,
            "event_median_seconds": event_median,
            "standalone_median_seconds": standalone_median,
            "wall_fraction": event_median / standalone_median,
            "gate_max": 0.25,
            "D0_peak_rss_kib_max": 1_048_576,
            "D1_peak_rss_kib_max": 2_097_152,
            "pass": (
                event_median / standalone_median <= 0.25
                and d0_rss <= 1_048_576
                and d1_rss <= 2_097_152
            ),
        },
        "held_out_D4": {
            **held, "wall_seconds_first": d4_wall, "peak_rss_kib": d4_rss
        },
        "interpretation": (
            "Publish the exact event/cohort engine. The locked D4 heterogeneity threshold failed, "
            "so do not claim biological generalization from the AKR1A1 event."
        ),
    }
    (args.out / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
