#!/usr/bin/env python3
"""Shared fail-closed primitives for the prospective archive-root gate.

This module intentionally does not import Gravlax.  In particular, its BLAKE3 and `.aie`
directory parser independently recompute the v2 root and compressed-payload commitments.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import struct
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO, Iterable, Mapping, Sequence


class GateError(RuntimeError):
    """A fail-closed gate or artifact-contract violation."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise GateError(message)


def require_exact_keys(value: Mapping[str, Any], keys: set[str], label: str) -> None:
    require(set(value) == keys, f"{label}: fields changed: {sorted(value)}")


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def load_json(path: Path) -> Any:
    try:
        raw = path.read_bytes()
        require(raw.endswith(b"\n"), f"{path}: JSON lacks a final newline")
        return json.loads(raw)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise GateError(f"cannot load JSON {path}: {error}") from error


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


# ---- Minimal independent BLAKE3 implementation ---------------------------------------------

_IV = (
    0x6A09E667,
    0xBB67AE85,
    0x3C6EF372,
    0xA54FF53A,
    0x510E527F,
    0x9B05688C,
    0x1F83D9AB,
    0x5BE0CD19,
)
_PERMUTATION = (2, 6, 3, 10, 7, 0, 4, 13, 1, 11, 12, 5, 9, 14, 15, 8)
_CHUNK_START = 1
_CHUNK_END = 2
_PARENT = 4
_ROOT = 8
_MASK32 = (1 << 32) - 1


def _rotr32(value: int, count: int) -> int:
    return ((value >> count) | (value << (32 - count))) & _MASK32


def _g(state: list[int], a: int, b: int, c: int, d: int, mx: int, my: int) -> None:
    state[a] = (state[a] + state[b] + mx) & _MASK32
    state[d] = _rotr32(state[d] ^ state[a], 16)
    state[c] = (state[c] + state[d]) & _MASK32
    state[b] = _rotr32(state[b] ^ state[c], 12)
    state[a] = (state[a] + state[b] + my) & _MASK32
    state[d] = _rotr32(state[d] ^ state[a], 8)
    state[c] = (state[c] + state[d]) & _MASK32
    state[b] = _rotr32(state[b] ^ state[c], 7)


def _compress(
    chaining_value: Sequence[int],
    block_words: Sequence[int],
    counter: int,
    block_len: int,
    flags: int,
) -> tuple[int, ...]:
    require(len(chaining_value) == 8 and len(block_words) == 16, "invalid BLAKE3 state")
    state = list(chaining_value) + list(_IV[:4]) + [
        counter & _MASK32,
        (counter >> 32) & _MASK32,
        block_len,
        flags,
    ]
    schedule = list(range(16))
    for _ in range(7):
        _g(state, 0, 4, 8, 12, block_words[schedule[0]], block_words[schedule[1]])
        _g(state, 1, 5, 9, 13, block_words[schedule[2]], block_words[schedule[3]])
        _g(state, 2, 6, 10, 14, block_words[schedule[4]], block_words[schedule[5]])
        _g(state, 3, 7, 11, 15, block_words[schedule[6]], block_words[schedule[7]])
        _g(state, 0, 5, 10, 15, block_words[schedule[8]], block_words[schedule[9]])
        _g(state, 1, 6, 11, 12, block_words[schedule[10]], block_words[schedule[11]])
        _g(state, 2, 7, 8, 13, block_words[schedule[12]], block_words[schedule[13]])
        _g(state, 3, 4, 9, 14, block_words[schedule[14]], block_words[schedule[15]])
        schedule = [schedule[index] for index in _PERMUTATION]
    return tuple(
        [state[index] ^ state[index + 8] for index in range(8)]
        + [state[index + 8] ^ chaining_value[index] for index in range(8)]
    )


def _block_words(block: bytes) -> tuple[int, ...]:
    require(len(block) <= 64, "BLAKE3 block exceeds 64 bytes")
    return struct.unpack("<16I", block.ljust(64, b"\0"))


@dataclass(frozen=True)
class _Output:
    input_cv: tuple[int, ...]
    block_words: tuple[int, ...]
    counter: int
    block_len: int
    flags: int

    def chaining_value(self) -> tuple[int, ...]:
        return _compress(
            self.input_cv, self.block_words, self.counter, self.block_len, self.flags
        )[:8]

    def root_digest(self) -> bytes:
        words = _compress(
            self.input_cv, self.block_words, 0, self.block_len, self.flags | _ROOT
        )
        return struct.pack("<16I", *words)[:32]


def _chunk_output(chunk: bytes, chunk_counter: int) -> _Output:
    require(len(chunk) <= 1024, "BLAKE3 chunk exceeds 1024 bytes")
    blocks = [chunk[index : index + 64] for index in range(0, len(chunk), 64)] or [b""]
    cv = _IV
    for index, block in enumerate(blocks[:-1]):
        flags = _CHUNK_START if index == 0 else 0
        cv = _compress(cv, _block_words(block), chunk_counter, len(block), flags)[:8]
    last_index = len(blocks) - 1
    flags = _CHUNK_END | (_CHUNK_START if last_index == 0 else 0)
    return _Output(tuple(cv), _block_words(blocks[-1]), chunk_counter, len(blocks[-1]), flags)


def _parent_output(left: Sequence[int], right: Sequence[int]) -> _Output:
    return _Output(tuple(_IV), tuple(left) + tuple(right), 0, 64, _PARENT)


class Blake3:
    """Unkeyed streaming BLAKE3 with the standard 32-byte digest."""

    def __init__(self) -> None:
        self._buffer = bytearray()
        self._chunk_count = 0
        self._cv_stack: list[tuple[int, ...]] = []

    def _add_completed_chunk(self, chunk: bytes) -> None:
        cv = _chunk_output(chunk, self._chunk_count).chaining_value()
        self._chunk_count += 1
        total = self._chunk_count
        while total & 1 == 0:
            require(bool(self._cv_stack), "BLAKE3 parent stack underflow")
            cv = _parent_output(self._cv_stack.pop(), cv).chaining_value()
            total >>= 1
        self._cv_stack.append(cv)

    def update(self, data: bytes) -> "Blake3":
        self._buffer.extend(data)
        # Retain one complete chunk: it may be the final Output node, whose full state is needed
        # to set ROOT rather than merely its chaining value.
        while len(self._buffer) > 1024:
            chunk = bytes(self._buffer[:1024])
            del self._buffer[:1024]
            self._add_completed_chunk(chunk)
        return self

    def digest(self) -> bytes:
        output = _chunk_output(bytes(self._buffer), self._chunk_count)
        for left in reversed(self._cv_stack):
            output = _parent_output(left, output.chaining_value())
        return output.root_digest()

    def hexdigest(self) -> str:
        return self.digest().hex()


def blake3_bytes(data: bytes) -> str:
    return Blake3().update(data).hexdigest()


def blake3_file(path: Path, offset: int = 0, length: int | None = None) -> str:
    digest = Blake3()
    with path.open("rb") as handle:
        handle.seek(offset)
        remaining = length
        while remaining is None or remaining > 0:
            request = 1024 * 1024 if remaining is None else min(1024 * 1024, remaining)
            block = handle.read(request)
            if not block:
                break
            digest.update(block)
            if remaining is not None:
                remaining -= len(block)
        require(remaining in (None, 0), f"{path}: short read while hashing")
    return digest.hexdigest()


def self_test_blake3() -> None:
    vectors = {
        b"": "af1349b9f5f9a1a6a0404dea36dcc9499bcb25c9adc112b7cc9a93cae41f3262",
        b"abc": "6437b3ac38465133ffb63b75273a8db548c558465d79db03fd359c6cd5bd9d85",
    }
    for payload, expected in vectors.items():
        require(blake3_bytes(payload) == expected, "independent BLAKE3 self-test failed")
    # Exercise the chunk tree and streaming boundary independently of the published short vectors.
    payload = bytes(index % 251 for index in range(4097))
    one = Blake3().update(payload).hexdigest()
    streaming = Blake3()
    for index in range(0, len(payload), 137):
        streaming.update(payload[index : index + 137])
    require(streaming.hexdigest() == one, "streaming BLAKE3 self-test failed")


# ---- Independent `.aie` v1/v2 audit ---------------------------------------------------------

AIE_MAGIC = b"AIE0"
AIE_DIRECTORY_MAGIC = b"AIED"
V1_FOOTER_BYTES = 12
V2_FOOTER_BYTES = 44
MAX_V2_SECTIONS = 1_000_000
MAX_V2_RAW_SECTION_BYTES = 512 * 1024 * 1024
MAX_V2_COMPRESSED_SECTION_BYTES = MAX_V2_RAW_SECTION_BYTES + 16 * 1024 * 1024
MAX_V2_DIRECTORY_BYTES = 4 + MAX_V2_SECTIONS * (1 + 255 + 24 + 32)
MAX_V2_COMPRESSION_RATIO = 4_096
V2_COMPRESSION_SLACK = 64 * 1024 * 1024


@dataclass(frozen=True)
class SectionEntry:
    name: str
    offset: int
    raw_len: int
    compressed_len: int
    committed_blake3: str | None
    payload_offset: int


@dataclass(frozen=True)
class ArchiveDirectory:
    path: Path
    file_bytes: int
    version: int
    directory_offset: int
    directory_bytes: bytes
    root: str | None
    sections: tuple[SectionEntry, ...]


def _read_exact(handle: BinaryIO, size: int, label: str) -> bytes:
    value = handle.read(size)
    require(len(value) == size, f"truncated {label}")
    return value


def parse_archive(path: Path, root_domain: bytes, verify_payloads: bool) -> ArchiveDirectory:
    path = path.resolve()
    file_bytes = path.stat().st_size
    require(file_bytes >= 20, f"{path}: truncated archive")
    with path.open("rb") as handle:
        header = _read_exact(handle, 8, "archive header")
        require(header[:4] == AIE_MAGIC, f"{path}: wrong archive magic")
        version = struct.unpack("<I", header[4:])[0]
        require(version in (1, 2), f"{path}: unsupported audit version {version}")
        footer_bytes = V2_FOOTER_BYTES if version == 2 else V1_FOOTER_BYTES
        require(file_bytes >= 8 + footer_bytes, f"{path}: truncated footer")
        footer_offset = file_bytes - footer_bytes
        handle.seek(footer_offset)
        footer = _read_exact(handle, footer_bytes, "archive footer")
        directory_offset = struct.unpack("<Q", footer[:8])[0]
        root = footer[8:40].hex() if version == 2 else None
        require(footer[-4:] == AIE_DIRECTORY_MAGIC, f"{path}: wrong directory magic")
        require(9 <= directory_offset <= footer_offset - 4, f"{path}: invalid directory offset")
        handle.seek(directory_offset)
        directory_len = footer_offset - directory_offset
        if version == 2:
            require(directory_len <= MAX_V2_DIRECTORY_BYTES, f"{path}: v2 directory exceeds safety limit")
        directory_bytes = _read_exact(handle, directory_len, "directory")

        require(len(directory_bytes) >= 4, f"{path}: truncated directory count")
        count = struct.unpack("<I", directory_bytes[:4])[0]
        if version == 2:
            require(count <= MAX_V2_SECTIONS, f"{path}: v2 section count exceeds safety limit")
        entry_floor = 25 + (32 if version == 2 else 0)
        require(count <= (len(directory_bytes) - 4) // entry_floor, f"{path}: impossible section count")
        cursor = 4
        names: set[str] = set()
        sections: list[SectionEntry] = []
        for index in range(count):
            require(cursor < len(directory_bytes), f"{path}: directory entry {index} is truncated")
            name_len = directory_bytes[cursor]
            cursor += 1
            require(name_len > 0, f"{path}: empty section name")
            end = cursor + name_len
            require(end + 24 + (32 if version == 2 else 0) <= len(directory_bytes), f"{path}: truncated directory entry")
            try:
                name = directory_bytes[cursor:end].decode("utf-8")
            except UnicodeDecodeError as error:
                raise GateError(f"{path}: invalid UTF-8 section name") from error
            cursor = end
            require(name not in names, f"{path}: duplicate section {name}")
            names.add(name)
            offset, raw_len, compressed_len = struct.unpack("<QQQ", directory_bytes[cursor : cursor + 24])
            cursor += 24
            if version == 2:
                require(raw_len <= MAX_V2_RAW_SECTION_BYTES, f"{path}: {name} raw length exceeds safety limit")
                require(
                    compressed_len <= MAX_V2_COMPRESSED_SECTION_BYTES,
                    f"{path}: {name} compressed length exceeds safety limit",
                )
                permitted_raw = compressed_len * MAX_V2_COMPRESSION_RATIO + V2_COMPRESSION_SLACK
                require(raw_len <= permitted_raw, f"{path}: {name} declares an unsafe compression ratio")
            committed = None
            if version == 2:
                committed = directory_bytes[cursor : cursor + 32].hex()
                cursor += 32
            payload_offset = offset + 1 + name_len + 16
            sections.append(
                SectionEntry(name, offset, raw_len, compressed_len, committed, payload_offset)
            )
        require(cursor == len(directory_bytes), f"{path}: trailing directory bytes")

        expected_offset = 8
        for section in sections:
            require(section.offset == expected_offset, f"{path}: noncanonical extent before {section.name}")
            handle.seek(section.offset)
            inline_name_len = _read_exact(handle, 1, f"{section.name} inline name length")[0]
            inline_name = _read_exact(handle, inline_name_len, f"{section.name} inline name")
            lengths = _read_exact(handle, 16, f"{section.name} inline lengths")
            inline_raw, inline_compressed = struct.unpack("<QQ", lengths)
            require(inline_name == section.name.encode(), f"{path}: inline name mismatch for {section.name}")
            require(inline_raw == section.raw_len, f"{path}: inline raw length mismatch for {section.name}")
            require(inline_compressed == section.compressed_len, f"{path}: inline compressed length mismatch for {section.name}")
            expected_offset = section.payload_offset + section.compressed_len
            require(expected_offset < directory_offset, f"{path}: section exceeds section area")
        require(expected_offset + 1 == directory_offset, f"{path}: gap or overlap before directory")
        handle.seek(expected_offset)
        require(_read_exact(handle, 1, "section terminator") == b"\0", f"{path}: missing section terminator")

    if version == 2:
        preimage = root_domain + header + struct.pack("<Q", directory_offset) + directory_bytes
        observed_root = blake3_bytes(preimage)
        require(observed_root == root, f"{path}: archive root mismatch")
    if verify_payloads:
        for section in sections:
            observed = blake3_file(path, section.payload_offset, section.compressed_len)
            if version == 2:
                require(observed == section.committed_blake3, f"{path}: payload digest mismatch for {section.name}")
    return ArchiveDirectory(path, file_bytes, version, directory_offset, directory_bytes, root, tuple(sections))


def compare_v1_v2(v1: Path, v2: Path, root_domain: bytes) -> dict[str, Any]:
    before = parse_archive(v1, root_domain, verify_payloads=False)
    # Payload equality is checked byte-for-byte below. The Gravlax full-inspection command is
    # separately mandatory in the driver and verifies every committed BLAKE3 at native speed;
    # the independent Python BLAKE3 here is intentionally reserved for the small root preimage.
    after = parse_archive(v2, root_domain, verify_payloads=False)
    require(before.version == 1 and after.version == 2, "comparison is not v1 -> v2")
    require(before.directory_offset == after.directory_offset, "migration changed the section area")
    require(len(before.sections) == len(after.sections), "migration changed section count")
    rows: list[dict[str, Any]] = []
    with before.path.open("rb") as left, after.path.open("rb") as right:
        for old, new in zip(before.sections, after.sections, strict=True):
            require(
                (old.name, old.offset, old.raw_len, old.compressed_len)
                == (new.name, new.offset, new.raw_len, new.compressed_len),
                f"migration changed section metadata at {old.name}",
            )
            left.seek(old.payload_offset)
            right.seek(new.payload_offset)
            remaining = old.compressed_len
            sha = hashlib.sha256()
            independent_blake3 = Blake3()
            while remaining:
                size = min(1024 * 1024, remaining)
                old_bytes = _read_exact(left, size, f"{old.name} v1 payload")
                new_bytes = _read_exact(right, size, f"{new.name} v2 payload")
                require(old_bytes == new_bytes, f"migration changed compressed payload {old.name}")
                sha.update(old_bytes)
                independent_blake3.update(old_bytes)
                remaining -= size
            payload_blake3 = independent_blake3.hexdigest()
            require(
                payload_blake3 == new.committed_blake3,
                f"migration committed the wrong compressed-payload BLAKE3 for {old.name}",
            )
            rows.append(
                {
                    "name": old.name,
                    "raw_bytes": old.raw_len,
                    "compressed_bytes": old.compressed_len,
                    "compressed_sha256": sha.hexdigest(),
                    "compressed_blake3": payload_blake3,
                    "committed_compressed_blake3": new.committed_blake3,
                }
            )
    expected_growth = 32 * len(rows) + 32
    require(after.file_bytes - before.file_bytes == expected_growth, "v2 migration growth is not the fixed layout cost")
    return {
        "schema": "gravlax.archive-root-section-comparison.v1",
        "source": str(before.path),
        "sealed": str(after.path),
        "source_bytes": before.file_bytes,
        "sealed_bytes": after.file_bytes,
        "added_bytes": expected_growth,
        "section_count": len(rows),
        "directory_offset": after.directory_offset,
        "archive_root": after.root,
        "sections": rows,
    }


# ---- Commands, GNU time, and complete-artifact manifests ------------------------------------

def parse_wall_seconds(value: str) -> float:
    fields = [float(field) for field in value.split(":")]
    require(1 <= len(fields) <= 3, f"invalid GNU-time wall value {value!r}")
    seconds = fields[-1]
    if len(fields) >= 2:
        seconds += 60 * fields[-2]
    if len(fields) == 3:
        seconds += 3600 * fields[-3]
    return seconds


def parse_gnu_time(path: Path, require_success: bool = True) -> dict[str, Any]:
    text = path.read_text()
    patterns = {
        "command": r'^\s*Command being timed: "(.*)"$',
        "user_seconds": r"^\s*User time \(seconds\): ([0-9.]+)$",
        "system_seconds": r"^\s*System time \(seconds\): ([0-9.]+)$",
        "wall": r"^\s*Elapsed \(wall clock\) time.*: ([0-9:.]+)$",
        "max_rss_kib": r"^\s*Maximum resident set size \(kbytes\): (\d+)$",
        "fs_inputs": r"^\s*File system inputs: (\d+)$",
        "fs_outputs": r"^\s*File system outputs: (\d+)$",
        "exit_status": r"^\s*Exit status: (\d+)$",
    }
    found = {key: re.findall(pattern, text, re.MULTILINE) for key, pattern in patterns.items()}
    for key, values in found.items():
        require(len(values) == 1, f"{path}: expected one {key} record, found {len(values)}")
    record = {
        "command": found["command"][0],
        "user_seconds": float(found["user_seconds"][0]),
        "system_seconds": float(found["system_seconds"][0]),
        "wall_seconds": parse_wall_seconds(found["wall"][0]),
        "max_rss_kib": int(found["max_rss_kib"][0]),
        "filesystem_inputs": int(found["fs_inputs"][0]),
        "filesystem_outputs": int(found["fs_outputs"][0]),
        "exit_status": int(found["exit_status"][0]),
    }
    if require_success:
        require(record["exit_status"] == 0, f"{path}: timed command failed")
    require(record["wall_seconds"] >= 0 and record["max_rss_kib"] > 0, f"{path}: invalid resource record")
    return record


def run_timed(
    label_dir: Path,
    command: Sequence[str],
    *,
    time_binary: Path,
    cwd: Path,
    environment: Mapping[str, str],
) -> None:
    require(not label_dir.exists(), f"refusing to overwrite command directory {label_dir}")
    label_dir.mkdir(parents=True)
    command = [str(token) for token in command]
    require(bool(command) and all("\0" not in token for token in command), "invalid command argv")
    write_json(
        label_dir / "command.json",
        {
            "schema": "gravlax.archive-root-command.v1",
            "argv": command,
            "cwd": str(cwd.resolve()),
            "environment": dict(sorted(environment.items())),
        },
    )
    with (label_dir / "stdout.txt").open("wb") as stdout, (label_dir / "stderr.txt").open("wb") as stderr:
        completed = subprocess.run(
            [str(time_binary), "-v", "-o", str(label_dir / "time.txt"), *command],
            cwd=cwd,
            env={**os.environ, **environment},
            stdout=stdout,
            stderr=stderr,
            check=False,
        )
    require(completed.returncode == 0, f"command failed ({completed.returncode}): {command}")
    parse_gnu_time(label_dir / "time.txt")


def tree_digest(root: Path) -> dict[str, Any]:
    require(root.is_dir(), f"missing output directory {root}")
    rows = []
    for path in sorted(root.rglob("*")):
        require(not path.is_symlink(), f"symlink not allowed in artifact tree: {path}")
        if path.is_file():
            relative = path.relative_to(root).as_posix()
            require("\t" not in relative and "\n" not in relative, "invalid artifact path")
            rows.append({"path": relative, "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    canonical = b"".join(
        row["path"].encode() + b"\0" + str(row["bytes"]).encode() + b"\0" + row["sha256"].encode() + b"\n"
        for row in rows
    )
    return {"files": rows, "tree_sha256": hashlib.sha256(canonical).hexdigest()}


def write_artifact_manifest(run_dir: Path) -> dict[str, Any]:
    manifest_path = run_dir / "artifact-manifest.tsv"
    require(not manifest_path.exists(), "artifact manifest already exists")
    rows = []
    paths = sorted(
        run_dir.rglob("*"),
        key=lambda path: path.relative_to(run_dir).as_posix(),
    )
    for path in paths:
        require(not path.is_symlink(), f"symlink not allowed in run tree: {path}")
        if path.is_file():
            relative = path.relative_to(run_dir).as_posix()
            require(relative != manifest_path.name, "artifact manifest recursion")
            require("\t" not in relative and "\n" not in relative, "invalid artifact path")
            rows.append((relative, path.stat().st_size, sha256_file(path)))
    with manifest_path.open("w") as handle:
        handle.write("relative_path\tbytes\tsha256\n")
        for relative, size, digest in rows:
            handle.write(f"{relative}\t{size}\t{digest}\n")
    return {
        "files": len(rows),
        "bytes": sum(size for _, size, _ in rows),
        "manifest_sha256": sha256_file(manifest_path),
    }


def validate_artifact_manifest(run_dir: Path) -> dict[str, Any]:
    manifest_path = run_dir / "artifact-manifest.tsv"
    lines = manifest_path.read_text().splitlines()
    require(lines and lines[0] == "relative_path\tbytes\tsha256", "artifact manifest header changed")
    recorded: dict[str, tuple[int, str]] = {}
    for line_number, line in enumerate(lines[1:], 2):
        fields = line.split("\t")
        require(len(fields) == 3, f"artifact manifest line {line_number} is malformed")
        relative, raw_size, digest = fields
        require(relative not in recorded and relative != manifest_path.name, "duplicate/recursive artifact")
        require(re.fullmatch(r"[0-9a-f]{64}", digest) is not None, "invalid artifact SHA-256")
        require(re.fullmatch(r"0|[1-9][0-9]*", raw_size) is not None, "invalid artifact size")
        recorded[relative] = (int(raw_size), digest)
    require(list(recorded) == sorted(recorded), "artifact manifest is not path-sorted")
    actual: set[str] = set()
    for path in run_dir.rglob("*"):
        require(not path.is_symlink(), f"symlink not allowed in run tree: {path}")
        if path.is_file() and path != manifest_path:
            actual.add(path.relative_to(run_dir).as_posix())
    require(actual == set(recorded), f"run artifact set changed; missing={sorted(set(recorded)-actual)}, extra={sorted(actual-set(recorded))}")
    for relative, (size, digest) in recorded.items():
        path = run_dir / relative
        require(path.stat().st_size == size, f"artifact size changed: {relative}")
        require(sha256_file(path) == digest, f"artifact digest changed: {relative}")
    return {
        "files": len(recorded),
        "bytes": sum(size for size, _ in recorded.values()),
        "manifest_sha256": sha256_file(manifest_path),
    }


def median(values: Iterable[float]) -> float:
    ordered = sorted(values)
    require(bool(ordered), "median of empty values")
    middle = len(ordered) // 2
    return ordered[middle] if len(ordered) % 2 else (ordered[middle - 1] + ordered[middle]) / 2
