#!/usr/bin/env python3
"""Evaluate the prospectively locked FNBP1 normal-brain confirmation."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from collections import defaultdict
from pathlib import Path


EVENT_ID = (
    "cassette:chr9:129908999-129923843:129924026-129924959:"
    "129908999-129924959"
)
JUNCTIONS = {
    "include_left": ("chr9", 129909000, 129923843),
    "include_right": ("chr9", 129924027, 129924959),
    "skip": ("chr9", 129909000, 129924959),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def event(path: Path) -> dict:
    data = json.loads(path.read_text())
    return next(item for item in data["events"] if item["id"] == EVENT_ID)


def star_sj_counts(directory: Path) -> dict[str, int]:
    selected_rows = {}
    by_coordinate = {coordinate: name for name, coordinate in JUNCTIONS.items()}
    with (directory / "features.tsv").open() as handle:
        for index, line in enumerate(handle, 1):
            fields = line.split()
            coordinate = (fields[0], int(fields[1]), int(fields[2]))
            if coordinate in by_coordinate:
                selected_rows[index] = by_coordinate[coordinate]
    if set(selected_rows.values()) != set(JUNCTIONS):
        raise SystemExit("target STAR-SJ rows are incomplete")
    counts = defaultdict(int)
    dimensions_seen = False
    with (directory / "matrix.mtx").open() as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.split()
            if not dimensions_seen:
                dimensions_seen = True
                continue
            row, _column, value = map(int, fields)
            if row in selected_rows:
                counts[selected_rows[row]] += value
    return {name: counts[name] for name in JUNCTIONS}


def bam_junction_counts(samtools: Path, bam: Path) -> dict[str, int]:
    wanted = {
        (coordinate[1] - 1, coordinate[2]): name
        for name, coordinate in JUNCTIONS.items()
    }
    molecules: dict[str, set[tuple[str, str]]] = {
        name: set() for name in JUNCTIONS
    }
    process = subprocess.Popen(
        [str(samtools), "view", str(bam)], stdout=subprocess.PIPE, text=True
    )
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.rstrip().split("\t")
        flag = int(fields[1])
        if flag & 0x804:
            continue
        cursor = int(fields[3]) - 1
        junctions = []
        for length_text, operation in re.findall(r"(\d+)([MIDNSHP=X])", fields[5]):
            length = int(length_text)
            if operation == "N":
                junctions.append((cursor, cursor + length))
                cursor += length
            elif operation in "MD=X":
                cursor += length
        tags = {
            field[:2]: field[5:]
            for field in fields[11:]
            if len(field) >= 5 and field[2:5] == ":Z:"
        }
        barcode = tags.get("CB", tags.get("CR", "")).removesuffix("-1")
        umi = tags.get("UB", tags.get("UR", ""))
        for coordinate in junctions:
            if coordinate in wanted:
                molecules[wanted[coordinate]].add((barcode, umi))
    if process.wait() != 0:
        raise SystemExit("samtools failed during direct BAM junction audit")
    return {name: len(values) for name, values in molecules.items()}


def fraction(counts: dict[str, int]) -> float:
    return counts["include_right"] / (counts["include_right"] + counts["skip"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--samtools", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    preparation = json.loads((args.run / "preparation.json").read_text())
    cellranger_event = event(args.run / "events.json")
    star_event = event(args.run / "star-events.json")
    star_counts = star_sj_counts(args.run / "star-targeted-t8/Solo.out/SJ/raw")
    bam_counts = bam_junction_counts(args.samtools, args.run / "fnbp1.called.bam")
    archive_totals = star_event["totals"]

    checks = {
        "event_present_at_catalogue_support_2": star_event["present"],
        "called_nuclei_informative_UMIs_at_least_50": (
            archive_totals["informative_umis"] >= 50
        ),
        "archive_cassette_inclusion_at_least_0_75": (
            archive_totals["usage_fraction"] >= 0.75
        ),
        "STAR_SJ_include_right_fraction_at_least_0_50": fraction(star_counts) >= 0.50,
        "archive_and_STAR_direction_agree": (
            archive_totals["usage_fraction"] > 0.5 and fraction(star_counts) > 0.5
        ),
    }
    result = {
        "schema": "gravlax.fnbp1-normal-brain-confirmation.v1",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "locked_event": EVENT_ID,
        "dataset": {
            "name": "10k Human Brain Nuclei, Chromium GEM-X Epi Multiome",
            "provider": "10x Genomics",
            "published": "2026-04-20",
            "url": "https://www.10xgenomics.com/datasets/multiome-gemx-10k-human-brain",
            "called_nuclei": preparation["called_nuclei"],
        },
        "inputs": {
            "filtered_feature_bc_matrix_h5": {
                "bytes": (args.run / "filtered_feature_bc_matrix.h5").stat().st_size,
                "sha256": sha256(args.run / "filtered_feature_bc_matrix.h5"),
            },
            "gex_possorted_bam": {
                "bytes": (args.run / "gex_possorted_bam.bam").stat().st_size,
                "sha256": sha256(args.run / "gex_possorted_bam.bam"),
            },
            "gex_possorted_bam_bai": {
                "bytes": (args.run / "gex_possorted_bam.bam.bai").stat().st_size,
                "sha256": sha256(args.run / "gex_possorted_bam.bam.bai"),
            },
        },
        "preparation": preparation,
        "annotation_free_targeted_STAR_archive": {
            "totals": archive_totals,
            "fully_annotated_label": star_event["annotation"]["fully_annotated"],
        },
        "Cell_Ranger_alignment_archive": {"totals": cellranger_event["totals"]},
        "STAR_SJ_matrix": {
            "UMIs": star_counts,
            "include_right_fraction": fraction(star_counts),
        },
        "Cell_Ranger_BAM_direct": {
            "distinct_CB_UB": bam_counts,
            "include_right_fraction": fraction(bam_counts),
        },
        "gate_checks": checks,
        "interpretation": {
            "claim": "prospective replication of the predeclared FNBP1 normal-brain cassette program",
            "not_claimed": "novelty, cell-type specificity, mechanism, or genome-wide mapping-sensitivity replication",
            "targeted_STAR": (
                "annotation-free realignment of reads selected at the locked locus; validates "
                "alignment direction and reducer semantics, not reads mapped elsewhere by Cell Ranger"
            ),
        },
    }
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
