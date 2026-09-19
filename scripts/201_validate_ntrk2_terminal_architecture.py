#!/usr/bin/env python3
"""Map NTRK2 PolyASite-mixture counts to GENCODE coding terminal groups.

This audit deliberately stops short of claiming molecule-level isoform phasing.  It asks a
narrower question: do PolyASite candidates adjacent to the MANE/full-length coding terminal
group replace candidates adjacent to the strongest upstream, shorter-coding GENCODE terminal
group in the frozen Astro/NSC-to-mature-neuron contrast?
"""

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from statistics import median


ATTR = re.compile(r'(\w+) "([^"]+)"')


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def attributes(text: str):
    return dict(ATTR.findall(text))


def read_transcripts(gtf: Path, gene_name: str):
    transcripts = {}
    cds_bases = defaultdict(int)
    with gtf.open() as handle:
        for line in handle:
            if line.startswith("#") or f'gene_name "{gene_name}"' not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            attrs = attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue
            if fields[2] == "transcript":
                transcripts[transcript_id] = {
                    "transcript_id": transcript_id,
                    "transcript_name": attrs.get("transcript_name", transcript_id),
                    "transcript_type": attrs.get("transcript_type", ""),
                    "chrom": fields[0],
                    "strand": fields[6],
                    "endpoint": int(fields[4]) if fields[6] == "+" else int(fields[3]),
                    "gencode_primary": 'tag "GENCODE_Primary"' in fields[8],
                    "mane_select": 'tag "MANE_Select"' in fields[8],
                }
            elif fields[2] == "CDS":
                cds_bases[transcript_id] += int(fields[4]) - int(fields[3]) + 1
    coding = []
    for transcript_id, row in transcripts.items():
        if row["transcript_type"] == "protein_coding" and cds_bases[transcript_id] > 0:
            row["cds_bases"] = cds_bases[transcript_id]
            coding.append(row)
    if not coding:
        raise ValueError(f"no coding transcripts found for {gene_name}")
    if len({(row["chrom"], row["strand"]) for row in coding}) != 1:
        raise ValueError(f"{gene_name} coding transcripts do not share one locus and strand")
    return sorted(coding, key=lambda row: row["endpoint"])


def cluster_transcripts(transcripts, gap: int):
    groups = []
    for transcript in transcripts:
        if not groups or transcript["endpoint"] - groups[-1][-1]["endpoint"] > gap:
            groups.append([])
        groups[-1].append(transcript)
    result = []
    for index, members in enumerate(groups):
        result.append(
            {
                "group_index": index,
                "endpoint_min": min(row["endpoint"] for row in members),
                "endpoint_max": max(row["endpoint"] for row in members),
                "median_cds_bases": median(row["cds_bases"] for row in members),
                "contains_gencode_primary": any(row["gencode_primary"] for row in members),
                "contains_mane_select": any(row["mane_select"] for row in members),
                "transcripts": members,
            }
        )
    return result


def select_groups(groups):
    mane = [group for group in groups if group["contains_mane_select"]]
    if len(mane) != 1:
        raise ValueError(f"expected exactly one MANE terminal group, observed {len(mane)}")
    full = mane[0]
    proximal = [
        group
        for group in groups
        if group["endpoint_max"] < full["endpoint_min"]
        and group["contains_gencode_primary"]
        and group["median_cds_bases"] < full["median_cds_bases"]
    ]
    if not proximal:
        raise ValueError("no upstream GENCODE Primary shorter-coding terminal group")
    # The upstream group with the most coding transcript support is the documented truncated
    # receptor family. Ties resolve to the group closest to the MANE/full-length endpoint.
    short = max(proximal, key=lambda group: (len(group["transcripts"]), group["endpoint_max"]))
    return short, full


def assign_sites(path: Path, gene_name: str, groups, radius: int):
    rows = []
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["gene_name"] != gene_name:
                continue
            position = int(row["position"])
            choices = []
            for group in groups:
                distance = min(abs(position - item["endpoint"]) for item in group["transcripts"])
                if distance <= radius:
                    choices.append((distance, group["group_index"]))
            row["position"] = position
            row["terminal_group"] = min(choices)[1] if choices else None
            rows.append(row)
    if not rows:
        raise ValueError(f"no {gene_name} rows in {path}")
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-tsv", type=Path, required=True)
    parser.add_argument("--gene", default="NTRK2")
    parser.add_argument("--transcript-end-gap", type=int, default=5000)
    parser.add_argument("--site-radius", type=int, default=2000)
    args = parser.parse_args()

    transcripts = read_transcripts(args.gtf, args.gene)
    groups = cluster_transcripts(transcripts, args.transcript_end_gap)
    short, full = select_groups(groups)
    sites = assign_sites(args.sites, args.gene, groups, args.site_radius)

    with args.genes.open() as handle:
        gene_rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["gene_name"] == args.gene]
    if len(gene_rows) != 1:
        raise ValueError(f"expected exactly one {args.gene} gene row, observed {len(gene_rows)}")
    gene_test = gene_rows[0]
    summary = json.loads(args.summary.read_text())

    count_columns = [column for column in sites[0] if ":" in column]
    samples = sorted({column.split(":", 1)[0] for column in count_columns})
    conditions = sorted({column.split(":", 1)[1] for column in count_columns})
    if conditions != ["astro_nsc", "mature_neuron"]:
        raise ValueError(f"unexpected conditions: {conditions}")

    rows = []
    for sample in samples:
        values = {}
        for condition in conditions:
            column = f"{sample}:{condition}"
            short_count = sum(
                float(site[column]) for site in sites if site["terminal_group"] == short["group_index"]
            )
            full_count = sum(
                float(site[column]) for site in sites if site["terminal_group"] == full["group_index"]
            )
            total = short_count + full_count
            if total <= 0:
                raise ValueError(f"zero classified terminal evidence for {column}")
            values[condition] = {
                "shorter_coding_expected_umis": short_count,
                "full_length_like_expected_umis": full_count,
                "full_length_like_fraction": full_count / total,
            }
        rows.append(
            {
                "sample": sample,
                "astro_nsc_shorter_coding_expected_umis": values["astro_nsc"]["shorter_coding_expected_umis"],
                "astro_nsc_full_length_like_expected_umis": values["astro_nsc"]["full_length_like_expected_umis"],
                "astro_nsc_full_length_like_fraction": values["astro_nsc"]["full_length_like_fraction"],
                "mature_neuron_shorter_coding_expected_umis": values["mature_neuron"]["shorter_coding_expected_umis"],
                "mature_neuron_full_length_like_expected_umis": values["mature_neuron"]["full_length_like_expected_umis"],
                "mature_neuron_full_length_like_fraction": values["mature_neuron"]["full_length_like_fraction"],
                "fraction_difference": values["mature_neuron"]["full_length_like_fraction"]
                - values["astro_nsc"]["full_length_like_fraction"],
            }
        )

    pooled = {}
    for condition in conditions:
        short_count = sum(row[f"{condition}_shorter_coding_expected_umis"] for row in rows)
        full_count = sum(row[f"{condition}_full_length_like_expected_umis"] for row in rows)
        pooled[condition] = {
            "shorter_coding_expected_umis": short_count,
            "full_length_like_expected_umis": full_count,
            "full_length_like_fraction": full_count / (short_count + full_count),
        }

    selected = {short["group_index"], full["group_index"]}
    result = {
        "schema": "gravlax.ntrk2-terminal-architecture.v1",
        "role": "exploratory structure-aware validation of the protocol-aware PolyASite mixture result",
        "gene": args.gene,
        "semantics": {
            "claim": "PolyASite candidates adjacent to GENCODE coding transcript ends; not molecule-level isoform phasing",
            "terminal_group_rule": f"coding transcript endpoints single-link clustered within {args.transcript_end_gap} bp",
            "site_assignment_rule": f"nearest terminal group within {args.site_radius} bp; ties resolve by distance then genomic order",
            "short_group_rule": "upstream GENCODE Primary group with shorter median CDS and greatest coding-transcript support",
            "full_group_rule": "the unique terminal group containing the GENCODE MANE Select transcript",
        },
        "inputs": {
            "gtf_sha256": digest(args.gtf),
            "sites_sha256": digest(args.sites),
            "genes_sha256": digest(args.genes),
            "summary_sha256": digest(args.summary),
        },
        "analysis_context": {
            "engine_schema": summary["schema"],
            "reported_genes": summary["inference"]["reported_genes"],
            "global_directional_threshold_pass": summary["inference"]["registered_global_directional_threshold_pass"],
            "ntrk2_gene_test": {
                key: gene_test[key]
                for key in ["n_sites", "eligible_samples", "effect_B_minus_A", "t", "p", "p_signflip", "q", "concordant", "lodo_stable", "reported"]
            },
        },
        "terminal_groups": [
            {
                **{key: value for key, value in group.items() if key != "transcripts"},
                "role": (
                    "proximal_shorter_coding" if group["group_index"] == short["group_index"]
                    else "distal_mane_full_length_like" if group["group_index"] == full["group_index"]
                    else "other"
                ),
                "assigned_sites": [site["position"] for site in sites if site["terminal_group"] == group["group_index"]],
                "transcripts": group["transcripts"],
            }
            for group in groups
        ],
        "unassigned_sites": [site["position"] for site in sites if site["terminal_group"] is None],
        "classified_site_count": sum(site["terminal_group"] in selected for site in sites),
        "total_ntrk2_site_count": len(sites),
        "per_sample": rows,
        "pooled": pooled,
        "direction": {
            "positive_samples": sum(row["fraction_difference"] > 0 for row in rows),
            "samples": len(rows),
            "median_fraction_difference": median(row["fraction_difference"] for row in rows),
            "exact_two_sided_signflip_p": 2 / (2 ** len(rows)),
        },
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    with args.out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
