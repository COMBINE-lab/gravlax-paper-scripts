#!/usr/bin/env python3
"""Run named Rust predicates and independent end-to-end Gate-B corruption fixtures.

Positive geometry/deduplication evidence comes from exact, source-controlled Rust tests. Format,
binding, block-context, reconstruction, and container failures are independently regenerated from
the real routed collection, exercised through the release CLI, identified, and deleted. Only
compact commands, streams, and the result manifest remain in the Gate-B run tree.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import struct
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import zstandard

from archive_root_gate_common import GateError, blake3_bytes, load_json, require, sha256_file, write_json
from jshape_routing_gate_common import (
    EXTERNAL_FIXTURE_CASES,
    PANEL_SHA256,
    REQUIRED_FIXTURE_CASES,
    RUST_FIXTURE_TESTS,
)


EMPTY_SHA256 = hashlib.sha256(b"").hexdigest()
COLLECTION_MAGIC = b"GRVLXCOL"
COLLECTION_VERSION = 4


@dataclass(frozen=True)
class Section:
    name: str
    raw_len: int
    compressed_len: int
    digest: bytes
    payload_offset: int
    compressed: bytes


class Cursor:
    def __init__(self, raw: bytes | bytearray):
        self.raw = raw
        self.position = 0

    def take(self, count: int) -> bytes:
        require(count >= 0 and self.position + count <= len(self.raw), "truncated fixture structure")
        result = bytes(self.raw[self.position : self.position + count])
        self.position += count
        return result

    def byte(self) -> int:
        return self.take(1)[0]

    def varint_location(self) -> tuple[int, int, int]:
        start, value, shift = self.position, 0, 0
        while True:
            byte = self.byte()
            value |= (byte & 0x7F) << shift
            if byte < 0x80:
                return value, start, self.position
            shift += 7
            require(shift < 70, "fixture varint is too long")

    def varint(self) -> int:
        return self.varint_location()[0]

    def string(self) -> str:
        size = self.varint()
        return self.take(size).decode("utf-8")


def put_varint(out: bytearray, value: int) -> None:
    require(type(value) is int and value >= 0, "cannot encode negative fixture varint")
    while value >= 0x80:
        out.append((value & 0x7F) | 0x80)
        value >>= 7
    out.append(value)


def parse_collection(data: bytes) -> list[Section]:
    require(len(data) >= 48 and data[:8] == COLLECTION_MAGIC, "fixture is not a collection")
    version, count = struct.unpack_from("<II", data, 8)
    require(version == COLLECTION_VERSION and 1 <= count <= 1_000_000, "fixture collection version/count differs")
    cursor = 48
    rows = []
    for _ in range(count):
        require(cursor < len(data), "truncated collection directory")
        name_len = data[cursor]
        cursor += 1
        require(name_len > 0 and cursor + name_len + 48 <= len(data), "malformed collection directory")
        name = data[cursor : cursor + name_len].decode("utf-8")
        cursor += name_len
        raw_len, compressed_len = struct.unpack_from("<QQ", data, cursor)
        cursor += 16
        digest = data[cursor : cursor + 32]
        cursor += 32
        rows.append((name, raw_len, compressed_len, digest))
    expected_root = bytes.fromhex(blake3_bytes(data[:16] + data[48:cursor]))
    require(data[16:48] == expected_root, "fixture collection directory root differs")
    payload_offset = cursor
    sections = []
    for name, raw_len, compressed_len, digest in rows:
        end = payload_offset + compressed_len
        require(end <= len(data), "truncated collection payload")
        sections.append(Section(name, raw_len, compressed_len, digest, payload_offset, data[payload_offset:end]))
        payload_offset = end
    require(payload_offset == len(data), "collection has trailing bytes")
    require(sections[0].name == "manifest" and len({row.name for row in sections}) == len(sections), "collection section order/names differ")
    return sections


def raw_section(section: Section) -> bytes:
    try:
        raw = zstandard.ZstdDecompressor().decompress(section.compressed, max_output_size=section.raw_len)
    except zstandard.ZstdError as error:
        raise GateError(f"cannot decompress fixture section {section.name}: {error}") from error
    require(len(raw) == section.raw_len, f"fixture section {section.name} raw length differs")
    require(bytes.fromhex(blake3_bytes(raw)) == section.digest, f"fixture section {section.name} digest differs")
    return raw


def repack(data: bytes, replacements: Mapping[str, bytes]) -> bytes:
    sections = parse_collection(data)
    rows = []
    compressor = zstandard.ZstdCompressor(level=9)
    for section in sections:
        if section.name in replacements:
            raw = replacements[section.name]
            compressed = compressor.compress(raw)
            digest = bytes.fromhex(blake3_bytes(raw))
            rows.append((section.name, len(raw), compressed, digest))
        else:
            rows.append((section.name, section.raw_len, section.compressed, section.digest))
    require(set(replacements) <= {row.name for row in sections}, "replacement names are absent")
    header = COLLECTION_MAGIC + struct.pack("<II", COLLECTION_VERSION, len(rows))
    directory = bytearray()
    for name, raw_len, compressed, digest in rows:
        encoded = name.encode()
        require(1 <= len(encoded) <= 255, "fixture section name is invalid")
        directory.append(len(encoded))
        directory.extend(encoded)
        directory.extend(struct.pack("<QQ", raw_len, len(compressed)))
        directory.extend(digest)
    root = bytes.fromhex(blake3_bytes(header + directory))
    return bytes(header + root + directory + b"".join(row[2] for row in rows))


def binding_locations(manifest: bytes) -> list[dict[str, int]]:
    cursor = Cursor(manifest)
    base = cursor.byte()
    require(base in {0, 1}, "manifest base flag differs")
    if base:
        cursor.string()
        cursor.string()
    for _ in range(2):
        present = cursor.byte()
        require(present in {0, 1}, "manifest optional-string flag differs")
        if present:
            cursor.string()
    for _ in range(cursor.varint()):
        cursor.string()
    cursor.string()
    bindings = []
    for archive_ordinal in range(cursor.varint()):
        cursor.string()
        cursor.string()
        for _ in range(8):
            cursor.varint()
        for _ in range(3):
            cursor.string()
        for _ in range(cursor.varint()):
            for _ in range(7):
                cursor.varint()
        require(cursor.byte() == 1, "routed manifest has an archive without a binding")
        codec, codec_start, codec_end = cursor.varint_location()
        source_root = cursor.position
        cursor.take(32)
        shapes_digest = cursor.position
        cursor.take(32)
        cursor.varint()
        for _ in range(cursor.varint()):
            cursor.varint()
            cursor.varint()
            cursor.string()
        bindings.append(
            {
                "archive_ordinal": archive_ordinal,
                "codec": codec,
                "codec_start": codec_start,
                "codec_end": codec_end,
                "source_root": source_root,
                "shapes_digest": shapes_digest,
            }
        )
    require(bindings, "routed fixture has no bindings")
    return bindings


def decode_block(raw: bytes) -> dict[str, Any]:
    require(len(raw) >= 28 and raw[:4] == b"JSHR", "fixture route block header differs")
    version, archive, n_shapes, first, last, n_spans = struct.unpack_from("<6I", raw, 4)
    cursor = Cursor(raw)
    cursor.position = 28
    spans = []
    previous_span = 0
    for index in range(n_spans):
        code = cursor.varint()
        span = code if index == 0 else previous_span + code
        pairs = []
        previous = None
        for _ in range(cursor.varint()):
            shape_code, donor_code = cursor.varint(), cursor.varint()
            if previous is None:
                pair = [shape_code, donor_code]
            elif shape_code == 0:
                pair = [previous[0], previous[1] + donor_code]
            else:
                pair = [previous[0] + shape_code, donor_code]
            pairs.append(pair)
            previous = pair
        spans.append({"span": span, "pairs": pairs})
        previous_span = span
    require(cursor.position == len(raw), "fixture route block has trailing bytes")
    return {
        "version": version, "archive": archive, "n_shapes": n_shapes,
        "first": first, "last": last, "n_spans": n_spans, "spans": spans,
    }


def encode_block(block: Mapping[str, Any]) -> bytes:
    out = bytearray(b"JSHR")
    out.extend(struct.pack("<6I", block["version"], block["archive"], block["n_shapes"], block["first"], block["last"], block["n_spans"]))
    previous_span = 0
    for index, row in enumerate(block["spans"]):
        put_varint(out, row["span"] if index == 0 else row["span"] - previous_span)
        put_varint(out, len(row["pairs"]))
        previous = None
        for shape, donor in row["pairs"]:
            if previous is None:
                shape_code, donor_code = shape, donor
            elif shape == previous[0]:
                require(donor >= previous[1], "cannot fixture-encode a descending same-shape pair")
                shape_code, donor_code = 0, donor - previous[1]
            else:
                require(shape >= previous[0], "cannot fixture-encode a descending shape id")
                shape_code, donor_code = shape - previous[0], donor
            put_varint(out, shape_code)
            put_varint(out, donor_code)
            previous = [shape, donor]
        previous_span = row["span"]
    return bytes(out)


def mutate_collection(base: bytes, kind: str) -> bytes:
    sections = parse_collection(base)
    manifest_section = sections[0]
    manifest = bytearray(raw_section(manifest_section))
    route_section = next(section for section in sections if section.name.startswith("s."))
    route_raw = raw_section(route_section)
    if kind.startswith("manifest_"):
        binding = binding_locations(manifest)[0]
        if kind == "manifest_codec_version":
            require(binding["codec"] == 1 and binding["codec_end"] == binding["codec_start"] + 1, "fixture codec is not one byte")
            manifest[binding["codec_start"]] = 2
        elif kind == "manifest_source_root":
            manifest[binding["source_root"]] ^= 1
        elif kind == "manifest_shapes_digest":
            manifest[binding["shapes_digest"]] ^= 1
        else:
            raise GateError(f"unknown manifest mutation {kind}")
        return repack(base, {"manifest": bytes(manifest)})
    if kind.startswith("block_"):
        block = decode_block(route_raw)
        if kind == "block_shape_id":
            row = next(row for row in block["spans"] if row["pairs"])
            row["pairs"][-1][0] = block["n_shapes"]
        elif kind == "block_reconstruction_mismatch":
            changed = False
            for row in block["spans"]:
                for index, pair in enumerate(row["pairs"]):
                    following = row["pairs"][index + 1] if index + 1 < len(row["pairs"]) else None
                    if pair[1] + row["span"] < 0xFFFFFFFF and (following is None or following[0] != pair[0] or pair[1] + 1 < following[1]):
                        pair[1] += 1
                        changed = True
                        break
                if changed:
                    break
            require(changed, "could not construct a canonical reconstruction mismatch")
        elif kind == "block_archive_ordinal":
            block["archive"] += 1
        elif kind == "block_span_bucket":
            block["first"] += 1
        else:
            raise GateError(f"unknown block mutation {kind}")
        return repack(base, {route_section.name: encode_block(block)})
    if kind == "container_truncation":
        return base[:-1]
    if kind == "container_trailing_bytes":
        return base + b"x"
    if kind == "route_payload_corruption":
        value = bytearray(base)
        value[route_section.payload_offset + route_section.compressed_len // 2] ^= 1
        return bytes(value)
    raise GateError(f"unknown external fixture mutation {kind}")


def run_command(command: Sequence[str], cwd: Path, stdout: Path, stderr: Path) -> int:
    with stdout.open("wb") as out, stderr.open("wb") as err:
        return subprocess.run(command, cwd=cwd, stdout=out, stderr=err, check=False).returncode


def slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "-", value)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aie", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--fallback-collection", type=Path, required=True)
    parser.add_argument("--route-collection", type=Path, required=True)
    parser.add_argument("--route-chain", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--archive", action="append", default=[])
    args = parser.parse_args()

    aie, out_dir = args.aie.resolve(), args.out_dir.resolve()
    route_collection, panel = args.route_collection.resolve(), args.panel.resolve()
    protocol = load_json(args.protocol.resolve())
    code_repo = Path(protocol["gravlax"]["repo"]).resolve()
    code_commit = protocol["gravlax"]["commit"]
    require(aie.is_file() and os.access(aie, os.X_OK), "fixture executable is unavailable")
    require(route_collection.is_file() and panel.is_file(), "fixture collection/panel is unavailable")
    require(sha256_file(panel) == PANEL_SHA256, "fixture panel differs")
    require(not out_dir.exists(), f"refusing to overwrite {out_dir}")
    head = subprocess.run(["git", "-C", str(code_repo), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(code_repo), "status", "--porcelain"], capture_output=True, text=True, check=True).stdout
    require(head == code_commit and not dirty, "fixture Rust source is not the clean protocol commit")
    out_dir.mkdir(parents=True)

    rust_records = []
    case_rows: dict[str, dict[str, Any]] = {}
    for test_name, case_ids in RUST_FIXTURE_TESTS.items():
        name = slug(test_name)
        stdout, stderr = out_dir / f"rust-{name}.stdout.txt", out_dir / f"rust-{name}.stderr.txt"
        command = ["cargo", "test", "-p", "aie", "--bin", "aie", test_name, "--", "--exact", "--test-threads=1"]
        status = run_command(command, code_repo, stdout, stderr)
        require(status == 0, f"Rust fixture test failed: {test_name}")
        require(re.search(r"test result: ok\. 1 passed;", stdout.read_text()) is not None, f"Rust fixture did not run exactly one test: {test_name}")
        source = code_repo / "crates/aie/src" / ("shaperoute.rs" if test_name.startswith("shaperoute::") else "collectioncmd.rs")
        record = {
            "test_name": test_name, "command": command, "exit_status": status,
            "stdout": stdout.name, "stderr": stderr.name,
            "source": str(source), "source_bytes": source.stat().st_size,
            "source_sha256": sha256_file(source), "cases": list(case_ids),
        }
        rust_records.append(record)
        for case_id in case_ids:
            case_rows[case_id] = {
                "id": case_id,
                "expect": REQUIRED_FIXTURE_CASES[case_id],
                "evidence_kind": "rust_unit_test",
                "test_name": test_name,
                "command": command,
                "exit_status": status,
                "stdout": stdout.name,
                "stderr": stderr.name,
                "fixture_bytes": source.stat().st_size,
                "fixture_sha256": sha256_file(source),
                "mutation": f"source-controlled predicate in {test_name}",
                "rejected_before_scientific_output": REQUIRED_FIXTURE_CASES[case_id] == "reject",
                "reconstruction_exact": REQUIRED_FIXTURE_CASES[case_id] == "accept",
            }

    base = route_collection.read_bytes()
    for case_id, mutation_kind in EXTERNAL_FIXTURE_CASES.items():
        data = mutate_collection(base, mutation_kind)
        fixture = out_dir / f"{case_id}.aicollection"
        stdout, stderr = out_dir / f"{case_id}.stdout.txt", out_dir / f"{case_id}.stderr.txt"
        fixture.write_bytes(data)
        fixture_bytes, fixture_sha = fixture.stat().st_size, sha256_file(fixture)
        command = [str(aie), "collection", "inspect", str(fixture), "--verify-routes"]
        status = run_command(command, code_repo, stdout, stderr)
        fixture.unlink()
        case_rows[case_id] = {
            "id": case_id,
            "expect": "reject",
            "evidence_kind": "cli_mutation",
            "test_name": None,
            "command": command,
            "exit_status": status,
            "stdout": stdout.name,
            "stderr": stderr.name,
            "fixture_bytes": fixture_bytes,
            "fixture_sha256": fixture_sha,
            "mutation": mutation_kind,
            "rejected_before_scientific_output": status != 0 and stdout.stat().st_size == 0,
            "reconstruction_exact": False,
        }

    require(set(case_rows) == set(REQUIRED_FIXTURE_CASES), "fixture evidence coverage differs")
    write_json(
        out_dir / "result.json",
        {
            "schema": "gravlax.jshape-route-adversarial.v2",
            "binary_sha256": sha256_file(aie),
            "route_collection_sha256": sha256_file(route_collection),
            "panel_sha256": sha256_file(panel),
            "code_commit": code_commit,
            "tools": {
                "python": f"{os.sys.version_info.major}.{os.sys.version_info.minor}.{os.sys.version_info.micro}",
                "zstandard_python": zstandard.__version__,
            },
            "rust_tests": rust_records,
            "cases": [case_rows[case] for case in sorted(case_rows)],
        },
    )


if __name__ == "__main__":
    try:
        main()
    except (GateError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"Gate-B adversarial hook failed: {error}") from error
