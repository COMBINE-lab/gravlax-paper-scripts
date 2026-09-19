#!/usr/bin/env python3
"""Prepare the locked FNBP1 normal-brain confirmation from public 10x artifacts.

The public Cell Ranger ARC BAM is filtered to called nuclei at the fixed FNBP1 window.  Original
barcode/UMI tags and forward-orientation cDNA sequence are then reconstructed as a small paired
FASTQ for an annotation-free targeted STAR realignment.  Selecting the window affects which
reads are realigned, so this is an alignment/reducer audit at the predeclared event, not a
whole-transcriptome replication of mapping sensitivity.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import h5py


FNBP1_WINDOW = "chr9:129907500-129926000"


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5", type=Path, required=True)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--samtools", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.h5) as handle:
        called = [value.decode() for value in handle["matrix/barcodes"][:]]
    called_path = args.out / "called-barcodes.txt"
    whitelist_path = args.out / "called-whitelist.txt"
    called_path.write_text("".join(f"{barcode}\n" for barcode in called))
    whitelist_path.write_text(
        "".join(f"{barcode.removesuffix('-1')}\n" for barcode in called)
    )

    target_bam = args.out / "fnbp1.called.bam"
    subprocess.run(
        [
            str(args.samtools), "view", "-@", str(args.threads), "-b",
            "-D", f"CB:{called_path}", "-o", str(target_bam), str(args.bam), FNBP1_WINDOW,
        ],
        check=True,
    )
    subprocess.run(
        [str(args.samtools), "index", "-@", str(args.threads), str(target_bam)], check=True
    )

    r1_path = args.out / "fnbp1.R1.fastq"
    r2_path = args.out / "fnbp1.R2.fastq"
    process = subprocess.Popen(
        [str(args.samtools), "view", str(target_bam)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    emitted: set[str] = set()
    records = 0
    with r1_path.open("w") as r1, r2_path.open("w") as r2:
        for line in process.stdout:
            fields = line.rstrip().split("\t")
            name, flag_text, sequence, quality = fields[0], fields[1], fields[9], fields[10]
            flag = int(flag_text)
            if flag & 0x900 or name in emitted:
                continue
            tags = {
                field[:2]: field[5:]
                for field in fields[11:]
                if len(field) >= 5 and field[2:5] == ":Z:"
            }
            required = {"CR", "CY", "UR", "UY"}
            if required - tags.keys():
                raise SystemExit(f"{name} lacks tags {sorted(required - tags.keys())}")
            barcode_read = tags["CR"] + tags["UR"]
            barcode_quality = tags["CY"] + tags["UY"]
            if len(barcode_read) != 28 or len(barcode_read) != len(barcode_quality):
                raise SystemExit(f"{name} has unexpected barcode/UMI lengths")
            if flag & 0x10:
                sequence = reverse_complement(sequence)
                quality = quality[::-1]
            r1.write(f"@{name}\n{barcode_read}\n+\n{barcode_quality}\n")
            r2.write(f"@{name}\n{sequence}\n+\n{quality}\n")
            emitted.add(name)
            records += 1
    if process.wait() != 0:
        raise SystemExit("samtools view failed while reconstructing FASTQ")

    summary = {
        "schema": "gravlax.fnbp1-confirmation-preparation.v1",
        "called_nuclei": len(called),
        "window": FNBP1_WINDOW,
        "target_bam": target_bam.name,
        "target_primary_reads": records,
        "r1": r1_path.name,
        "r2": r2_path.name,
    }
    (args.out / "preparation.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
