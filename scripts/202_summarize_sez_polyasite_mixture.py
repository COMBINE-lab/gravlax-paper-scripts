#!/usr/bin/env python3
"""Validate the frozen PolyASite-mixture controls and emit publication-sized results."""

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_time(path: Path):
    text = path.read_text()
    wall = re.search(r"Elapsed \(wall clock\) time.*: ([0-9:.]+)", text)
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if not wall or not rss:
        raise ValueError(f"cannot parse GNU time record {path}")
    fields = [float(value) for value in wall.group(1).split(":")]
    seconds = fields[-1]
    if len(fields) >= 2:
        seconds += 60 * fields[-2]
    if len(fields) == 3:
        seconds += 3600 * fields[-3]
    return {"wall_seconds": seconds, "max_rss_kib": int(rss.group(1))}


def load_arm(root: Path, name: str):
    summary = json.loads((root / name / "summary.json").read_text())
    genes = list(csv.DictReader((root / name / "genes.tsv").open(), delimiter="\t"))
    return summary, {row["gene_id"]: row for row in genes}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--simulation", type=Path, required=True)
    parser.add_argument("--ntrk2", type=Path, required=True)
    parser.add_argument("--optimized", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--gravlax-commit", required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-tsv", type=Path, required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.gravlax_commit):
        raise SystemExit("--gravlax-commit must be a full lowercase 40-hex commit")

    arm_names = ["primary-gap24", "sensitivity-gap12", "sensitivity-gap48", "shuffled-gap24"]
    arms = {name: load_arm(args.run_root, name) for name in arm_names}
    times = {name: read_time(args.run_root / f"{name}.time.txt") for name in arm_names}
    legacy_time = read_time(args.run_root / "legacy.time.txt")

    primary_summary, primary_genes = arms["primary-gap24"]
    primary_reported = {
        gene_id: row for gene_id, row in primary_genes.items() if row["reported"] == "1"
    }
    stable = []
    for gene_id, row in primary_reported.items():
        sign = float(row["effect_B_minus_A"]) > 0
        present_and_same = True
        for name in ["sensitivity-gap12", "sensitivity-gap48"]:
            other = arms[name][1].get(gene_id)
            present_and_same &= other is not None and (float(other["effect_B_minus_A"]) > 0) == sign
        if present_and_same:
            stable.append(gene_id)

    heldout = primary_summary["accuracy"]["heldout_predictive_gain_over_uniform_bits_per_umi"]
    ntrk2 = json.loads(args.ntrk2.read_text())
    simulation = json.loads(args.simulation.read_text())
    optimized_time = read_time(Path(str(args.optimized) + ".time.txt"))
    for table in ["sites.tsv", "genes.tsv", "fragment-kernel.tsv"]:
        if sha256(args.optimized / table) != sha256(args.run_root / "primary-gap24" / table):
            raise SystemExit(f"optimized {table} differs from the frozen primary result")
    wall_speedup = legacy_time["wall_seconds"] / optimized_time["wall_seconds"]
    stability_fraction = len(stable) / len(primary_reported)

    gates = {
        "primary_wall_seconds_le_90": optimized_time["wall_seconds"] <= 90,
        "primary_rss_gib_le_8": optimized_time["max_rss_kib"] <= 8 * 1024 * 1024,
        "legacy_wall_speedup_ge_2": wall_speedup >= 2,
        "heldout_kernel_gain_positive_every_donor": all(row["bits_per_umi"] > 0 for row in heldout),
        "shuffle_zero_reported": arms["shuffled-gap24"][0]["inference"]["reported_genes"] == 0,
        "shuffle_global_fail": not arms["shuffled-gap24"][0]["inference"]["registered_global_directional_threshold_pass"],
        "neighboring_gap_sign_stability_ge_80_percent": stability_fraction >= 0.8,
        "legacy_outputs_complete": (args.run_root / "legacy.complete").is_file(),
    }
    if not all(gates.values()):
        raise SystemExit(f"one or more frozen controls failed: {gates}")

    reported_rows = []
    for gene_id, row in sorted(
        primary_reported.items(), key=lambda item: (float(item[1]["q"]), item[1]["gene_name"])
    ):
        reported_rows.append(
            {
                "gene_id": gene_id,
                "gene_name": row["gene_name"],
                "effect_B_minus_A": float(row["effect_B_minus_A"]),
                "q": float(row["q"]),
                "eligible_samples": int(row["eligible_samples"]),
                "concordant": int(row["concordant"]),
                "lodo_stable": int(row["lodo_stable"]),
                "stable_sign_gap12_gap48": int(gene_id in stable),
            }
        )

    result = {
        "schema": "gravlax.sez-polyasite-mixture-result.v1",
        "status": "PASS_METHOD_AMENDED_EXPLORATORY",
        "confirmatory_status": "the preregistered global effect-size gate failed and is not rescued",
        "gates": gates,
        "primary": {
            "recurrent_sites": primary_summary["mixture"]["recurrent_sites"],
            "assigned_expected_umis": primary_summary["mixture"]["assigned_umis"],
            "tested_genes": primary_summary["inference"]["tested_genes"],
            "reported_genes": primary_summary["inference"]["reported_genes"],
            "positive_donor_medians": primary_summary["inference"]["positive_donor_medians"],
            "across_donor_median_effect": primary_summary["inference"]["across_donor_median"],
            "registered_global_directional_threshold_pass": primary_summary["inference"]["registered_global_directional_threshold_pass"],
            "heldout_kernel_gain": heldout,
            "wall_seconds": optimized_time["wall_seconds"],
            "max_rss_kib": optimized_time["max_rss_kib"],
            "pre_hotpath_wall_seconds": times["primary-gap24"]["wall_seconds"],
            "hotpath_speedup": times["primary-gap24"]["wall_seconds"] / optimized_time["wall_seconds"],
        },
        "controls": {
            "reported_genes_by_arm": {
                name: arms[name][0]["inference"]["reported_genes"] for name in arm_names
            },
            "primary_genes_same_sign_both_neighboring_gaps": len(stable),
            "primary_gene_sign_stability_fraction": stability_fraction,
            "legacy_wall_seconds": legacy_time["wall_seconds"],
            "legacy_max_rss_kib": legacy_time["max_rss_kib"],
            "wall_speedup_vs_eight_legacy_scans": wall_speedup,
            "simulation": simulation,
        },
        "biological_story": {
            "headline": "NTRK2 shifts from a proximal shorter-coding terminal group in Astro/NSCs to the distal MANE/full-length-like terminal group in mature neurons",
            "interpretation_guard": "terminal-group usage is not molecule-level isoform phasing",
            "ntrk2": ntrk2,
        },
        "reported_genes": reported_rows,
        "input_digests": {
            "gravlax_commit": args.gravlax_commit,
            "manifest_sha256": sha256(args.manifest),
            "primary_summary_sha256": sha256(args.run_root / "primary-gap24" / "summary.json"),
            "primary_sites_sha256": sha256(args.run_root / "primary-gap24" / "sites.tsv"),
            "primary_genes_sha256": sha256(args.run_root / "primary-gap24" / "genes.tsv"),
            "simulation_sha256": sha256(args.simulation),
            "ntrk2_sha256": sha256(args.ntrk2),
            "optimized_binary_and_outputs_sha256": sha256(Path(str(args.optimized) + ".execution.sha256")),
        },
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    with args.out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(reported_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(reported_rows)


if __name__ == "__main__":
    main()
