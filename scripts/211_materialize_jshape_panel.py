#!/usr/bin/env python3
"""Materialize the frozen 96-junction Gate-B panel after a validated Gate-A PASS.

This script is intentionally unusable without the fail-closed Gate-A result.  It reads only the
rooted archives' chromosome dictionaries and junction catalogues, selects rows by the prospective
domain-separated BLAKE3 rule, and refuses to overwrite any output.  Running it freezes the Gate-B
panel; do not run it during Gate-A development.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import struct
import subprocess
from pathlib import Path
from typing import Any

from archive_root_gate_common import (
    GateError,
    blake3_bytes,
    load_json,
    parse_archive,
    require,
    sha256_file,
    write_json,
)


ROOT_DOMAIN = b"gravlax-aie-directory-root-v2\0"
SELECTION_DOMAIN = b"gravlax-jshape-panel-v1\0"
SAMPLES = tuple("ABCDEFGH")
STRATA = (
    ("one_archive", 1, 1),
    ("two_to_three_archives", 2, 3),
    ("four_to_seven_archives", 4, 7),
    ("eight_archives", 8, 8),
)


def decode_varint(data: bytes, cursor: int) -> tuple[int, int]:
    value = 0
    for shift in range(0, 70, 7):
        require(cursor < len(data), "truncated junction varint")
        byte = data[cursor]
        cursor += 1
        require(shift < 64 or byte == 0, "junction varint exceeds u64")
        value |= (byte & 0x7F) << shift
        if byte < 0x80:
            require(value <= (1 << 64) - 1, "junction varint exceeds u64")
            return value, cursor
    raise GateError("junction varint exceeds ten bytes")


def decompress_section(archive: Path, section: Any, zstd: Path) -> bytes:
    with archive.open("rb") as handle:
        handle.seek(section.payload_offset)
        compressed = handle.read(section.compressed_len)
    require(len(compressed) == section.compressed_len, f"{archive}: short compressed-section read")
    completed = subprocess.run(
        [str(zstd), "-q", "-d", "-c"],
        input=compressed,
        capture_output=True,
        check=False,
    )
    require(completed.returncode == 0, f"zstd failed for {archive}:{section.name}: {completed.stderr.decode(errors='replace')}")
    require(len(completed.stdout) == section.raw_len, f"{archive}:{section.name}: decompressed length differs")
    return completed.stdout


def decode_chroms(raw: bytes, label: str) -> tuple[str, ...]:
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise GateError(f"{label}: chromosome dictionary is not UTF-8") from error
    chroms = tuple(text.splitlines())
    require(bool(chroms) and all(chroms), f"{label}: chromosome dictionary is empty or malformed")
    require(len(chroms) == len(set(chroms)), f"{label}: duplicate chromosome names")
    return chroms


def decode_junctions(raw: bytes, chroms: tuple[str, ...], label: str) -> list[tuple[str, int, int]]:
    rows = []
    cursor = 0
    last_chrom = None
    last_donor = 0
    previous = None
    while cursor < len(raw):
        chrom_id, cursor = decode_varint(raw, cursor)
        donor_delta, cursor = decode_varint(raw, cursor)
        span, cursor = decode_varint(raw, cursor)
        require(chrom_id < len(chroms), f"{label}: junction chromosome id is out of range")
        if chrom_id != last_chrom:
            last_chrom = chrom_id
            last_donor = 0
        donor = last_donor + donor_delta
        acceptor = donor + span
        require(donor <= 0xFFFFFFFF and acceptor <= 0xFFFFFFFF, f"{label}: junction coordinate overflows u32")
        last_donor = donor
        row = (chroms[chrom_id], donor, acceptor)
        coordinate_key = (chrom_id, donor, acceptor)
        require(previous is None or coordinate_key > previous, f"{label}: junction catalogue is not strictly sorted")
        previous = coordinate_key
        rows.append(row)
    return rows


def selection_digest(chrom: str, donor: int, acceptor: int) -> str:
    preimage = SELECTION_DOMAIN + chrom.encode("utf-8") + struct.pack("<II", donor, acceptor)
    return blake3_bytes(preimage)


def atomic_text(path: Path, value: str) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(value)
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gate-a-result", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--digest-out", type=Path, required=True)
    parser.add_argument("--metadata-out", type=Path, required=True)
    parser.add_argument("--zstd", type=Path, default=Path(os.environ.get("GRAVLAX_ZSTD", "zstd")))
    args = parser.parse_args()

    result_path = args.gate_a_result.resolve()
    out = args.out.resolve()
    digest_out = args.digest_out.resolve()
    metadata_out = args.metadata_out.resolve()
    zstd = args.zstd.resolve()
    for path in (out, digest_out, metadata_out):
        require(not path.exists(), f"refusing to overwrite {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
    require(zstd.is_file() and os.access(zstd, os.X_OK), f"zstd executable is unavailable: {zstd}")

    gate_a = load_json(result_path)
    require(gate_a.get("schema") == "gravlax.archive-root-gate-a.v1", "Gate-A result schema differs")
    require(gate_a.get("status") == "PASS", "Gate A has not passed")
    gates = gate_a.get("gates")
    require(type(gates) is dict and gates and all(value is True for value in gates.values()), "Gate-A result contains a failed gate")
    run = Path(gate_a["identity"]["run_dir"]).resolve()
    require(run.is_dir(), f"Gate-A run is unavailable: {run}")
    recorded_archives = gate_a["archive_root"]["archives"]
    require(type(recorded_archives) is dict and set(recorded_archives) >= set(SAMPLES), "Gate-A archive panel differs")

    presence: dict[tuple[str, int, int], set[str]] = {}
    common_chroms = None
    catalogue_counts = {}
    archive_records = []
    for sample in SAMPLES:
        archive = run / "archives" / f"{sample}.v2.aie"
        require(archive.is_file(), f"missing rooted archive {sample}: {archive}")
        recorded = recorded_archives[sample]
        require(archive.stat().st_size == recorded["sealed_bytes"], f"{sample}: rooted archive size differs")
        require(sha256_file(archive) == recorded["sealed_sha256"], f"{sample}: rooted archive SHA-256 differs")
        parsed = parse_archive(archive, ROOT_DOMAIN, verify_payloads=False)
        require(parsed.root == recorded["archive_root"], f"{sample}: rooted archive commitment differs")
        by_name = {entry.name: entry for entry in parsed.sections}
        require("chroms" in by_name and "index.junctions" in by_name, f"{sample}: panel input sections are missing")
        chroms = decode_chroms(decompress_section(archive, by_name["chroms"], zstd), sample)
        if common_chroms is None:
            common_chroms = chroms
        else:
            require(chroms == common_chroms, f"{sample}: chromosome dictionary differs")
        rows = decode_junctions(
            decompress_section(archive, by_name["index.junctions"], zstd), chroms, sample
        )
        catalogue_counts[sample] = len(rows)
        for row in rows:
            chrom, donor, acceptor = row
            if acceptor <= donor:
                continue
            presence.setdefault(row, set()).add(sample)
        archive_records.append(
            {
                "id": sample,
                "path": str(archive),
                "bytes": archive.stat().st_size,
                "sha256": recorded["sealed_sha256"],
                "archive_root": recorded["archive_root"],
                "junctions": len(rows),
            }
        )

    selected_rows = []
    eligible_counts = {}
    for stratum, minimum, maximum in STRATA:
        eligible = []
        for (chrom, donor, acceptor), sample_ids in presence.items():
            count = len(sample_ids)
            if minimum <= count <= maximum:
                digest = selection_digest(chrom, donor, acceptor)
                eligible.append((digest, chrom, donor, acceptor, tuple(sorted(sample_ids))))
        eligible.sort(key=lambda row: (row[0], row[1], row[2], row[3]))
        require(len(eligible) >= 24, f"{stratum}: only {len(eligible)} eligible junctions")
        eligible_counts[stratum] = len(eligible)
        for rank, (digest, chrom, donor, acceptor, sample_ids) in enumerate(eligible[:24], 1):
            selected_rows.append(
                (stratum, rank, digest, chrom, donor, acceptor, len(sample_ids), ",".join(sample_ids))
            )

    require(len(selected_rows) == 96, "panel does not contain exactly 96 rows")
    coordinates = [(row[3], row[4], row[5]) for row in selected_rows]
    require(len(coordinates) == len(set(coordinates)), "panel repeats a junction across strata")
    lines = [
        "stratum\trank\tselection_blake3\tchrom\tdonor\tacceptor\tarchives_present\tarchive_ids"
    ]
    lines.extend("\t".join(map(str, row)) for row in selected_rows)
    atomic_text(out, "\n".join(lines) + "\n")
    panel_sha256 = sha256_file(out)
    atomic_text(digest_out, f"{panel_sha256}  {out.name}\n")
    write_json(
        metadata_out,
        {
            "schema": "gravlax.jshape-prospective-panel.v1",
            "status": "FROZEN_AFTER_GATE_A_PASS_BEFORE_GATE_B_IMPLEMENTATION_OR_MEASUREMENT",
            "gate_a_result": {
                "path": str(result_path),
                "sha256": sha256_file(result_path),
                "gravlax_commit": gate_a["identity"]["gravlax_commit"],
                "paper_scripts_commit": gate_a["identity"]["paper_scripts_commit"],
            },
            "selection": {
                "domain_ascii": "gravlax-jshape-panel-v1\\0",
                "domain_hex": SELECTION_DOMAIN.hex(),
                "preimage": "domain || chrom_utf8 || donor_u32_le || acceptor_u32_le",
                "order": "smallest BLAKE3, then chrom/donor/acceptor",
                "rows_per_stratum": 24,
            },
            "panel": {"path": str(out), "bytes": out.stat().st_size, "sha256": panel_sha256, "rows": 96},
            "digest_file": {"path": str(digest_out), "sha256": sha256_file(digest_out)},
            "eligible_rows": eligible_counts,
            "catalogue_rows": catalogue_counts,
            "archives": archive_records,
            "zstd": {
                "path": str(zstd),
                "version": subprocess.run([str(zstd), "--version"], capture_output=True, text=True, check=True).stdout.strip(),
            },
        },
    )
    print(
        f"frozen 96-row panel at {out} (sha256 {panel_sha256}); Gate-B implementation may begin only after these artifacts are committed"
    )


if __name__ == "__main__":
    try:
        main()
    except GateError as error:
        raise SystemExit(f"junction-shape panel materialization failed: {error}") from error
