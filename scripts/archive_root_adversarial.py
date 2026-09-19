#!/usr/bin/env python3
"""Generate and exercise the complete Gate-A archive-root adversarial fixture panel.

Every large mutation is regenerated from one rooted archive, exercised once, identified by
SHA-256, and removed.  Only the commands, streams, and compact fixture manifest are retained.
The hook deliberately exits successfully even when a case behaves unexpectedly: the independent
Gate-A driver and reducer compare each recorded exit status with its frozen expectation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import struct
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

from archive_root_gate_common import (
    GateError,
    blake3_bytes,
    parse_archive,
    require,
    sha256_file,
    write_json,
)


MAX_V2_RAW_SECTION_BYTES = 512 * 1024 * 1024
COMPRESSION_SLACK = 64 * 1024 * 1024


@dataclass(frozen=True)
class EntryLayout:
    record_start: int
    record_end: int
    name_start: int
    name_len: int
    offset_position: int
    raw_position: int
    compressed_position: int
    digest_position: int
    section_offset: int
    payload_offset: int
    raw_len: int
    compressed_len: int
    name: str


def sha256_bytes(value: bytes | bytearray) -> str:
    return hashlib.sha256(value).hexdigest()


def directory_layout(data: bytes | bytearray) -> tuple[int, int, list[EntryLayout]]:
    require(data[:8] == b"AIE0\x02\x00\x00\x00", "fixture is not a v2 archive")
    require(len(data) >= 44 and data[-4:] == b"AIED", "fixture footer is malformed")
    directory_offset = struct.unpack_from("<Q", data, len(data) - 44)[0]
    footer_offset = len(data) - 44
    require(directory_offset + 4 <= footer_offset, "fixture directory is truncated")
    count = struct.unpack_from("<I", data, directory_offset)[0]
    cursor = directory_offset + 4
    entries = []
    for _ in range(count):
        record_start = cursor
        name_len = data[cursor]
        cursor += 1
        name_start = cursor
        raw_name = bytes(data[cursor : cursor + name_len])
        cursor += name_len
        offset_position = cursor
        section_offset, raw_len, compressed_len = struct.unpack_from("<QQQ", data, cursor)
        cursor += 24
        digest_position = cursor
        cursor += 32
        entries.append(
            EntryLayout(
                record_start=record_start,
                record_end=cursor,
                name_start=name_start,
                name_len=name_len,
                offset_position=offset_position,
                raw_position=offset_position + 8,
                compressed_position=offset_position + 16,
                digest_position=digest_position,
                section_offset=section_offset,
                payload_offset=section_offset + 1 + name_len + 16,
                raw_len=raw_len,
                compressed_len=compressed_len,
                name=raw_name.decode("utf-8"),
            )
        )
    require(cursor == footer_offset, "fixture directory has trailing bytes")
    return directory_offset, footer_offset, entries


def rewrite_root(data: bytearray, domain: bytes) -> None:
    directory_offset, footer_offset, _ = directory_layout(data)
    preimage = (
        domain
        + bytes(data[:8])
        + struct.pack("<Q", directory_offset)
        + bytes(data[directory_offset:footer_offset])
    )
    data[footer_offset + 8 : footer_offset + 40] = bytes.fromhex(blake3_bytes(preimage))


def make_archive(
    domain: bytes,
    sections: Sequence[tuple[bytes, int, bytes]],
) -> bytes:
    header = b"AIE0" + struct.pack("<I", 2)
    body = bytearray(header)
    directory_rows: list[tuple[bytes, int, int, int, bytes]] = []
    for name, raw_len, compressed in sections:
        require(1 <= len(name) <= 255, "synthetic section name is outside the v2 layout")
        offset = len(body)
        body.append(len(name))
        body.extend(name)
        body.extend(struct.pack("<QQ", raw_len, len(compressed)))
        body.extend(compressed)
        directory_rows.append(
            (name, offset, raw_len, len(compressed), bytes.fromhex(blake3_bytes(compressed)))
        )
    body.append(0)
    directory_offset = len(body)
    directory = bytearray(struct.pack("<I", len(directory_rows)))
    for name, offset, raw_len, compressed_len, digest in directory_rows:
        directory.append(len(name))
        directory.extend(name)
        directory.extend(struct.pack("<QQQ", offset, raw_len, compressed_len))
        directory.extend(digest)
    preimage = domain + header + struct.pack("<Q", directory_offset) + directory
    root = bytes.fromhex(blake3_bytes(preimage))
    return bytes(body + directory + struct.pack("<Q", directory_offset) + root + b"AIED")


def set_directory_and_inline_raw(data: bytearray, entry: EntryLayout, raw_len: int) -> None:
    struct.pack_into("<Q", data, entry.raw_position, raw_len)
    inline_raw_position = entry.section_offset + 1 + entry.name_len
    struct.pack_into("<Q", data, inline_raw_position, raw_len)


def mutate_entry_order(data: bytearray, domain: bytes) -> None:
    directory_offset, _, entries = directory_layout(data)
    require(len(entries) >= 2, "entry-order fixture requires two sections")
    first = bytes(data[entries[0].record_start : entries[0].record_end])
    second = bytes(data[entries[1].record_start : entries[1].record_end])
    middle = bytes(data[entries[0].record_end : entries[1].record_start])
    replacement = second + middle + first
    start = entries[0].record_start
    end = entries[1].record_end
    require(len(replacement) == end - start and start >= directory_offset, "entry swap changed directory length")
    data[start:end] = replacement
    rewrite_root(data, domain)


def mutate_swapped_payloads(data: bytearray) -> None:
    _, _, entries = directory_layout(data)
    by_length: dict[int, EntryLayout] = {}
    pair = None
    for entry in entries:
        previous = by_length.get(entry.compressed_len)
        if previous is not None and entry.compressed_len > 0:
            pair = (previous, entry)
            break
        by_length[entry.compressed_len] = entry
    require(pair is not None, "fixture has no equal-length payload pair")
    left, right = pair
    left_payload = bytes(data[left.payload_offset : left.payload_offset + left.compressed_len])
    right_payload = bytes(data[right.payload_offset : right.payload_offset + right.compressed_len])
    data[left.payload_offset : left.payload_offset + left.compressed_len] = right_payload
    data[right.payload_offset : right.payload_offset + right.compressed_len] = left_payload


def run_case(
    *,
    case_id: str,
    expect: str,
    mutation: str,
    data: bytes | bytearray,
    command_builder: Callable[[Path], list[str]],
    out_dir: Path,
) -> dict[str, object]:
    fixture = out_dir / f"{case_id}.aie"
    stdout_path = out_dir / f"{case_id}.stdout.txt"
    stderr_path = out_dir / f"{case_id}.stderr.txt"
    fixture.write_bytes(data)
    fixture_bytes = fixture.stat().st_size
    fixture_sha256 = sha256_file(fixture)
    command = command_builder(fixture)
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        completed = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    fixture.unlink()
    return {
        "id": case_id,
        "expect": expect,
        "command": command,
        "exit_status": completed.returncode,
        "stdout": stdout_path.name,
        "stderr": stderr_path.name,
        "fixture_bytes": fixture_bytes,
        "fixture_sha256": fixture_sha256,
        "mutation": mutation,
        "rejected_before_scientific_output": expect == "reject" and stdout_path.stat().st_size == 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aie", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--root-domain-hex", required=True)
    parser.add_argument("--fixture-archive", type=Path, required=True)
    parser.add_argument("--source-archive", type=Path, required=True)
    args = parser.parse_args()

    aie = args.aie.resolve()
    out_dir = args.out_dir.resolve()
    fixture_archive = args.fixture_archive.resolve()
    source_archive = args.source_archive.resolve()
    require(aie.is_file() and os.access(aie, os.X_OK), "fixture executable is unavailable")
    require(fixture_archive.is_file() and source_archive.is_file(), "fixture archives are unavailable")
    require(not out_dir.exists(), f"refusing to overwrite {out_dir}")
    try:
        domain = bytes.fromhex(args.root_domain_hex)
    except ValueError as error:
        raise GateError("root domain is not hexadecimal") from error
    require(domain == b"gravlax-aie-directory-root-v2\0", "root domain differs from Gate A")
    rooted = parse_archive(fixture_archive, domain, verify_payloads=False)
    require(rooted.version == 2 and len(rooted.sections) >= 2, "fixture archive is not a multi-section v2 archive")
    require(parse_archive(source_archive, domain, verify_payloads=False).version == 1, "source fixture is not v1")

    out_dir.mkdir(parents=True)
    base = fixture_archive.read_bytes()
    _, _, entries = directory_layout(base)
    first = entries[0]
    first_payload = bytes(base[first.payload_offset : first.payload_offset + first.compressed_len])

    inspect = lambda path: [str(aie), "inspect-archive", str(path), "--json"]
    inspect_full = lambda path: [
        str(aie), "inspect-archive", str(path), "--json", "--verify-content"
    ]
    query = lambda path: [
        str(aie),
        "query",
        str(path),
        "junction",
        "chr9:129925157-129927194",
        "--json",
        "--top",
        "0",
    ]
    cases: list[dict[str, object]] = []

    def add(
        case_id: str,
        expect: str,
        mutation: str,
        data: bytes | bytearray,
        command_builder: Callable[[Path], list[str]] = inspect,
    ) -> None:
        cases.append(
            run_case(
                case_id=case_id,
                expect=expect,
                mutation=mutation,
                data=data,
                command_builder=command_builder,
                out_dir=out_dir,
            )
        )

    add("zero_sections", "accept", "canonical v2 archive with zero sections", make_archive(domain, []))
    add(
        "one_section",
        "accept",
        "one copied valid compressed frame",
        make_archive(domain, [(b"meta", first.raw_len, first_payload)]),
        inspect_full,
    )
    add("many_sections", "accept", "unmodified production v2 archive", base)
    add(
        "unknown_optional_section",
        "accept",
        "one valid frame under an unknown optional name",
        make_archive(domain, [(b"unknown.optional", first.raw_len, first_payload)]),
        inspect_full,
    )
    add(
        "maximum_legal_name",
        "accept",
        "one valid frame under a 255-byte UTF-8 section name",
        make_archive(domain, [(b"x" * 255, first.raw_len, first_payload)]),
        inspect_full,
    )

    unselected = bytearray(base)
    target = next((entry for entry in entries if entry.name == "edges"), entries[-1])
    unselected[target.payload_offset] ^= 1
    add(
        "unselected_corruption_lazy_read",
        "accept",
        f"flip first compressed byte of unselected {target.name}; inspect root/directory only",
        unselected,
    )
    add(
        "unselected_corruption_full_verify",
        "reject",
        f"same unselected {target.name} corruption under full verification",
        unselected,
        inspect_full,
    )

    value = bytearray(base)
    value[0] ^= 1
    add("bad_magic", "reject", "flip archive magic byte", value)

    value = bytearray(base)
    struct.pack_into("<I", value, 4, 3)
    add("future_version", "reject", "replace version 2 with future version 3", value)

    value = bytearray(base)
    value[-36] ^= 1
    add("root_mutation", "reject", "flip archive-root byte", value)

    value = bytearray(base)
    directory_offset, _, current_entries = directory_layout(value)
    value[directory_offset + 4] ^= 1
    add("directory_byte_mutation", "reject", "flip committed directory byte without replacing root", value)

    value = bytearray(base)
    mutate_entry_order(value, domain)
    add("entry_order_mutation", "reject", "swap first two directory entries and recompute root", value)

    value = bytearray(base)
    _, _, current_entries = directory_layout(value)
    value[current_entries[0].name_start] ^= 1
    rewrite_root(value, domain)
    add("name_mutation", "reject", "change directory name, retain inline name, recompute root", value)

    value = bytearray(base)
    _, _, current_entries = directory_layout(value)
    struct.pack_into("<Q", value, current_entries[0].offset_position, current_entries[0].section_offset + 1)
    rewrite_root(value, domain)
    add("offset_mutation", "reject", "shift first directory extent by one and recompute root", value)

    value = bytearray(base)
    _, _, current_entries = directory_layout(value)
    struct.pack_into("<Q", value, current_entries[0].raw_position, current_entries[0].raw_len + 1)
    rewrite_root(value, domain)
    add("raw_length_mutation", "reject", "change only directory raw length and recompute root", value)

    value = bytearray(base)
    _, _, current_entries = directory_layout(value)
    struct.pack_into(
        "<Q", value, current_entries[0].compressed_position, current_entries[0].compressed_len + 1
    )
    rewrite_root(value, domain)
    add("compressed_length_mutation", "reject", "change only directory compressed length and recompute root", value)

    value = bytearray(base)
    _, _, current_entries = directory_layout(value)
    value[current_entries[0].digest_position] ^= 1
    rewrite_root(value, domain)
    add("payload_digest_mutation", "reject", "change payload commitment and recompute root", value, inspect_full)

    value = bytearray(base)
    _, _, current_entries = directory_layout(value)
    value[current_entries[0].section_offset + 1] ^= 1
    add("inline_header_mutation", "reject", "change inline name without changing directory", value)

    value = bytearray(base)
    directory_offset, _, _ = directory_layout(value)
    value[directory_offset - 1] = 1
    add("terminator_mutation", "reject", "replace the unique section terminator", value)

    value = bytearray(base)
    value[-1] ^= 1
    add("footer_mutation", "reject", "change footer magic", value)

    value = bytearray(base)
    _, _, current_entries = directory_layout(value)
    value[current_entries[0].payload_offset] ^= 1
    add("payload_mutation", "reject", "flip selected compressed payload byte", value, inspect_full)

    value = bytearray(base)
    mutate_swapped_payloads(value)
    add("swapped_payloads", "reject", "swap two equal-length payloads", value, inspect_full)

    value = bytearray(base)
    _, _, current_entries = directory_layout(value)
    second = current_entries[1]
    struct.pack_into("<Q", value, second.offset_position, current_entries[0].section_offset + 1)
    rewrite_root(value, domain)
    add("overlap", "reject", "move second directory extent inside first extent", value)

    value = bytearray(base)
    _, _, current_entries = directory_layout(value)
    struct.pack_into("<Q", value, current_entries[1].offset_position, current_entries[1].section_offset + 1)
    rewrite_root(value, domain)
    add("gap", "reject", "insert a one-byte logical gap before second extent", value)

    add("truncation", "reject", "remove final footer byte", base[:-1])
    add("trailing_bytes", "reject", "append one byte after footer", base + b"x")

    value = bytearray(base)
    _, _, current_entries = directory_layout(value)
    set_directory_and_inline_raw(value, current_entries[0], current_entries[0].raw_len + 1)
    rewrite_root(value, domain)
    add(
        "decompressed_length_mismatch",
        "reject",
        "change directory and inline raw length consistently, then recompute root",
        value,
        inspect_full,
    )

    value = bytearray(base)
    _, _, current_entries = directory_layout(value)
    set_directory_and_inline_raw(value, current_entries[0], MAX_V2_RAW_SECTION_BYTES + 1)
    rewrite_root(value, domain)
    add("allocation_bomb", "reject", "declare raw section one byte beyond fixed allocation ceiling", value)

    value = bytearray(base)
    _, _, current_entries = directory_layout(value)
    smallest = min(current_entries, key=lambda entry: entry.compressed_len)
    unsafe_raw = smallest.compressed_len * 4096 + COMPRESSION_SLACK + 1
    require(unsafe_raw <= MAX_V2_RAW_SECTION_BYTES, "cannot construct bounded compression-ratio fixture")
    set_directory_and_inline_raw(value, smallest, unsafe_raw)
    rewrite_root(value, domain)
    add("compression_bomb", "reject", "declare raw length one byte beyond fixed ratio+slack bound", value)

    value = bytearray(base)
    _, _, current_entries = directory_layout(value)
    selected = next(entry for entry in current_entries if entry.name == "meta")
    value[selected.payload_offset] ^= 1
    add(
        "selected_corruption_normal_read",
        "reject",
        "corrupt meta payload selected by an ordinary junction query",
        value,
        query,
    )

    write_json(
        out_dir / "result.json",
        {
            "schema": "gravlax.archive-root-adversarial.v1",
            "binary_sha256": sha256_file(aie),
            "source_archive": {
                "path": str(source_archive),
                "bytes": source_archive.stat().st_size,
                "sha256": sha256_file(source_archive),
            },
            "fixture_archive": {
                "path": str(fixture_archive),
                "bytes": fixture_archive.stat().st_size,
                "sha256": sha256_bytes(base),
                "archive_root": rooted.root,
            },
            "root_domain_hex": domain.hex(),
            "cases": cases,
        },
    )


if __name__ == "__main__":
    try:
        main()
    except GateError as error:
        raise SystemExit(f"archive-root adversarial hook failed: {error}") from error
