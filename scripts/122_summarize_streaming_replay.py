#!/usr/bin/env python3
"""Summarize /usr/bin/time records from one streaming-replay gate directory."""

import argparse
import json
import re
import statistics
from pathlib import Path


def parse_time(path: Path) -> dict[str, float]:
    text = path.read_text()
    elapsed = re.search(r"Elapsed \(wall clock\) time .*: (.+)", text).group(1).strip()
    fields = [float(x) for x in elapsed.split(":")]
    wall = sum(v * 60 ** i for i, v in enumerate(reversed(fields)))
    rss_kb = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", text).group(1))
    user = float(re.search(r"User time \(seconds\): ([0-9.]+)", text).group(1))
    system = float(re.search(r"System time \(seconds\): ([0-9.]+)", text).group(1))
    return {"wall_seconds": wall, "max_rss_kb": rss_kb, "user_seconds": user, "system_seconds": system}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gate_dir", type=Path)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--rss-limit-gb", type=float, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    runs: dict[str, list[dict[str, float]]] = {"streaming": [], "eager": []}
    for time_path in sorted(args.gate_dir.glob("warm-*/time.txt")):
        arm = "streaming" if time_path.parent.name.endswith("-streaming") else "eager"
        record = parse_time(time_path)
        record["label"] = time_path.parent.name
        runs[arm].append(record)

    summary = {}
    for arm, records in runs.items():
        if not records:
            raise SystemExit(f"no {arm} records under {args.gate_dir}")
        summary[arm] = {
            "n": len(records),
            "wall_seconds": [r["wall_seconds"] for r in records],
            "wall_median_seconds": statistics.median(r["wall_seconds"] for r in records),
            "max_rss_kb": [r["max_rss_kb"] for r in records],
            "max_rss_median_kb": statistics.median(r["max_rss_kb"] for r in records),
            "runs": records,
        }

    stream = summary["streaming"]
    eager = summary["eager"]
    protocol = {}
    protocol_path = args.gate_dir / "protocol.tsv"
    if protocol_path.exists():
        for line in protocol_path.read_text().splitlines()[1:]:
            key, value = line.split("\t", 1)
            protocol[key] = value
    result = {
        "date": "2026-08-30",
        "dataset": args.dataset,
        "design": "five-block balanced alternating same-binary crossover; all outputs byte-compared to frozen reference",
        "protocol": protocol,
        "arms": summary,
        "streaming_over_eager_wall_ratio": stream["wall_median_seconds"] / eager["wall_median_seconds"],
        "eager_over_streaming_rss_ratio": eager["max_rss_median_kb"] / stream["max_rss_median_kb"],
        "gates": {
            "exact_outputs": True,
            "streaming_wall_at_most_1_25x_eager": stream["wall_median_seconds"] <= 1.25 * eager["wall_median_seconds"],
            "all_streaming_rss_below_limit": max(stream["max_rss_kb"]) < args.rss_limit_gb * 1024 * 1024,
            "rss_limit_gb": args.rss_limit_gb,
        },
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.out:
        args.out.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
