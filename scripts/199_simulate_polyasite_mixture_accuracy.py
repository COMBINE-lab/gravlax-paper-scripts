#!/usr/bin/env python3
"""Deterministic recovery audit for the empirical-kernel PolyASite mixture model."""

import argparse
import csv
import json
import random
from collections import Counter
from pathlib import Path


def load_kernel(path: Path):
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    counts = [int(row["umis"]) for row in rows]
    starts = [int(row["distance_start"]) for row in rows]
    smooth = []
    for index in range(len(counts)):
        smooth.append(1.0 + sum(counts[max(0, index - 2) : min(len(counts), index + 3)]))
    total = sum(smooth)
    return starts, [value / total for value in smooth]


def likelihood(distance, starts, probabilities):
    first = starts[0]
    width = starts[1] - starts[0]
    index = (distance - first) // width
    if index < 0 or index >= len(probabilities):
        return 0.0
    return probabilities[index]


def estimate(endpoints, sites, starts, probabilities):
    theta = [1.0 / len(sites)] * len(sites)
    expected = [0.0] * len(sites)
    for _ in range(50):
        expected = [0.0] * len(sites)
        for endpoint, count in endpoints.items():
            weights = [
                theta[index] * likelihood(site - endpoint, starts, probabilities)
                for index, site in enumerate(sites)
            ]
            denominator = sum(weights)
            if denominator == 0:
                continue
            for index, weight in enumerate(weights):
                expected[index] += count * weight / denominator
        total = sum(expected)
        updated = [value / total for value in expected]
        change = max(abs(a - b) for a, b in zip(theta, updated))
        theta = updated
        if change < 1e-8:
            break
    return theta


def quantile(values, probability):
    values = sorted(values)
    return values[round((len(values) - 1) * probability)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-tsv", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260901)
    args = parser.parse_args()
    starts, probabilities = load_kernel(args.kernel)
    rng = random.Random(args.seed)
    distances = [start + 5 for start in starts]
    separations = [50, 100, 250, 500, 1000]
    depths = [50, 200, 1000]
    truths = [0.1, 0.25, 0.5, 0.75, 0.9]
    rows = []
    for separation in separations:
        for depth in depths:
            errors = []
            biases = []
            for truth in truths:
                for _ in range(args.replicates):
                    endpoints = Counter()
                    for _ in range(depth):
                        distal = rng.random() < truth
                        site = separation if distal else 0
                        distance = rng.choices(distances, weights=probabilities, k=1)[0]
                        endpoints[site - distance] += 1
                    estimate_distal = estimate(endpoints, [0, separation], starts, probabilities)[1]
                    errors.append(abs(estimate_distal - truth))
                    biases.append(estimate_distal - truth)
            rows.append(
                {
                    "separation_bp": separation,
                    "depth": depth,
                    "cases": len(errors),
                    "mae": sum(errors) / len(errors),
                    "median_ae": quantile(errors, 0.5),
                    "p95_ae": quantile(errors, 0.95),
                    "bias": sum(biases) / len(biases),
                }
            )
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    result = {
        "schema": "gravlax.polyasite-mixture-simulation.v1",
        "seed": args.seed,
        "replicates_per_truth": args.replicates,
        "truths": truths,
        "kernel": str(args.kernel),
        "rows": rows,
    }
    args.out_json.write_text(json.dumps(result, indent=2) + "\n")
    with args.out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
