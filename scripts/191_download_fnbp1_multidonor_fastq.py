#!/usr/bin/env python3
"""Acquire the registered GSE234790 paired FASTQs directly from ENA."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter
from pathlib import Path


def digest(path: Path, algorithm: str = "sha256") -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--ena", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--aria2c", default="aria2c")
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--connections", type=int, default=4)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    if manifest["status"] != "PROSPECTIVE_GATE_FROZEN_BEFORE_ANY_FNBP1_SEQUENCE_OR_EVENT_QUERY":
        raise SystemExit("manifest is not the frozen pre-event gate")
    expected_ena = manifest["dataset"]["processed_inputs"]["ena_run_manifest"]
    if digest(args.ena) != expected_ena["sha256"]:
        raise SystemExit("ENA manifest digest mismatch")
    args.out.mkdir(parents=True, exist_ok=True)

    files: list[dict[str, object]] = []
    runs: set[str] = set()
    donors: Counter[str] = Counter()
    with args.ena.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            match = re.search(r", ([A-H]), snRNA-seq$", row["sample_title"])
            if not match:
                raise SystemExit(f"cannot resolve donor from {row['sample_title']!r}")
            donor = match.group(1)
            run = row["run_accession"]
            urls = row["fastq_ftp"].split(";")
            sizes = row["fastq_bytes"].split(";")
            md5s = row["fastq_md5"].split(";")
            if not (len(urls) == len(sizes) == len(md5s) == 2):
                raise SystemExit(f"{run} does not have exactly two registered FASTQs")
            runs.add(run)
            donors[donor] += 1
            for mate, (url, size, md5) in enumerate(zip(urls, sizes, md5s), 1):
                if not re.fullmatch(r"[0-9a-f]{32}", md5):
                    raise SystemExit(f"invalid MD5 for {run} mate {mate}")
                logical = Path(f"donor-{donor}") / f"{run}_{mate}.fastq.gz"
                files.append(
                    {
                        "donor": donor,
                        "run": run,
                        "mate": mate,
                        "url": "https://" + url.removeprefix("ftp://").removeprefix("https://"),
                        "logical_path": str(logical),
                        "bytes": int(size),
                        "md5": md5,
                    }
                )
    if len(runs) != expected_ena["runs"] or sum(int(row["bytes"]) for row in files) != expected_ena["paired_fastq_bytes"]:
        raise SystemExit("ENA run/file totals differ from the frozen manifest")
    if sorted(donors) != list("ABCDEFGH"):
        raise SystemExit(f"unexpected donor set {sorted(donors)}")

    control = args.out / "aria2-input.txt"
    with control.open("w") as handle:
        for row in files:
            destination = args.out / str(row["logical_path"])
            destination.parent.mkdir(exist_ok=True)
            handle.write(f"{row['url']}\n")
            handle.write(f"  dir={destination.parent}\n")
            handle.write(f"  out={destination.name}\n")
            handle.write(f"  checksum=md5={row['md5']}\n")

    subprocess.run(
        [
            args.aria2c,
            "--input-file", str(control),
            "--continue=true",
            "--file-allocation=none",
            "--auto-file-renaming=false",
            "--allow-overwrite=false",
            "--check-integrity=true",
            "--max-tries=8",
            "--retry-wait=5",
            "--timeout=60",
            "--max-concurrent-downloads", str(args.jobs),
            "--split", str(args.connections),
            "--max-connection-per-server", str(args.connections),
            "--summary-interval=30",
        ],
        check=True,
    )

    verified: list[dict[str, object]] = []
    for row in files:
        path = args.out / str(row["logical_path"])
        if path.stat().st_size != row["bytes"]:
            raise SystemExit(f"size mismatch for {path}")
        observed_md5 = digest(path, "md5")
        if observed_md5 != row["md5"]:
            raise SystemExit(f"MD5 mismatch for {path}")
        verified.append(
            {
                "donor": row["donor"],
                "run": row["run"],
                "mate": row["mate"],
                "logical_path": row["logical_path"],
                "bytes": row["bytes"],
                "md5": observed_md5,
            }
        )
    result = {
        "schema": "gravlax.fnbp1-multidonor-fastq-acquisition.v1",
        "status": "PASS",
        "manifest_sha256": digest(args.manifest),
        "ena_manifest_sha256": digest(args.ena),
        "event_evidence_inspected": False,
        "donors": dict(sorted(donors.items())),
        "runs": len(runs),
        "files": verified,
        "total_bytes": sum(int(row["bytes"]) for row in verified),
    }
    (args.out / "acquisition.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in ("status", "donors", "runs", "total_bytes")}, indent=2))


if __name__ == "__main__":
    main()
