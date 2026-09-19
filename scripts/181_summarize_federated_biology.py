#!/usr/bin/env python3
"""Summarize the exploratory federation screen and audit selected loci against STAR SJ matrices."""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import re
import statistics
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable


PBMC = ("D0", "D1", "D3", "D4")
BRAIN = ("D2", "D2p")
ALL = ("D0", "D1", "D2", "D2p", "D3", "D4")
SJ_DIRS = {
    "D0": "runs/oracle/full/v49/Solo.out/SJ/raw",
    "D1": "runs/oracle/d1/v49/Solo.out/SJ/raw",
    "D2": "runs/oracle/d2/v49/Solo.out/SJ/raw",
    "D2p": "runs/oracle/d2p/v49/Solo.out/SJ/raw",
    "D3": "runs/ingest/d3/Solo.out/SJ/raw",
    "D4": "runs/ingest/d4/Solo.out/SJ/raw",
}
SELECTED_EVENTS = {
    "FYB1": "cassette:chr5:39124278-39125997:39126135-39127740:39124278-39127740",
    "CD47": "alt_acceptor:chr3:108047292-108049618:108047292-108057476",
    "FNBP1": "cassette:chr9:129908999-129923843:129924026-129924959:129908999-129924959",
    "MYL6": "cassette:chr12:56160320-56160625:56160670-56161386:56160320-56161386",
    "PPP1R12A": "cassette:chr12:79798584-79805591:79805768-79806165:79798584-79806165",
    "RPS24": "cassette:chr10:78037304-78040203:78040225-78040614:78037304-78040614",
}
SJ_TARGETS = {
    "FYB1_skip": ("chr5", 39124279, 39127740),
    "FYB1_include_right": ("chr5", 39126136, 39127740),
    "FNBP1_skip": ("chr9", 129909000, 129924959),
    "FNBP1_include_right": ("chr9", 129924027, 129924959),
    "MYL6_skip": ("chr12", 56160321, 56161386),
    "MYL6_include_right": ("chr12", 56160671, 56161386),
}
CHROMS = tuple(f"chr{number}" for number in range(1, 23)) + ("chrX", "chrY")


class StreamingJsonReader:
    """Incremental JSON decoder that bounds memory to one value plus a read buffer."""

    def __init__(self, path: Path, chunk_size: int = 1024 * 1024) -> None:
        self.path = path
        self.handle = path.open()
        self.chunk_size = chunk_size
        self.buffer = ""
        self.position = 0
        self.absolute_offset = 0
        self.eof = False
        self.decoder = json.JSONDecoder()

    def close(self) -> None:
        self.handle.close()

    def _fill(self) -> bool:
        if self.position:
            self.absolute_offset += self.position
            self.buffer = self.buffer[self.position :]
            self.position = 0
        chunk = self.handle.read(self.chunk_size)
        if not chunk:
            self.eof = True
            return False
        self.buffer += chunk
        return True

    def _skip_whitespace(self) -> None:
        while True:
            while self.position < len(self.buffer) and self.buffer[self.position].isspace():
                self.position += 1
            if self.position < len(self.buffer) or not self._fill():
                return

    def peek(self) -> str:
        self._skip_whitespace()
        if self.position >= len(self.buffer):
            raise ValueError(f"unexpected end of JSON in {self.path}")
        return self.buffer[self.position]

    def expect(self, token: str) -> None:
        actual = self.peek()
        if actual != token:
            offset = self.absolute_offset + self.position
            raise ValueError(
                f"expected {token!r} at character {offset} in {self.path}, found {actual!r}"
            )
        self.position += 1

    def decode(self) -> Any:
        self._skip_whitespace()
        start_offset = self.absolute_offset + self.position
        while True:
            try:
                value, end = self.decoder.raw_decode(self.buffer, self.position)
                self.position = end
                return value
            except json.JSONDecodeError as error:
                if self.eof or not self._fill():
                    raise ValueError(
                        f"invalid JSON value starting near character {start_offset} in "
                        f"{self.path}: {error.msg}"
                    ) from error

    def ensure_finished(self) -> None:
        self._skip_whitespace()
        if self.position != len(self.buffer) or (not self.eof and self._fill()):
            self._skip_whitespace()
            if self.position != len(self.buffer):
                offset = self.absolute_offset + self.position
                raise ValueError(f"trailing content at character {offset} in {self.path}")


def visit_event_document(path: Path, visitor: Callable[[dict[str, Any]], None]) -> dict:
    """Visit ``events`` one at a time and return the small top-level metadata."""

    reader = StreamingJsonReader(path)
    metadata: dict[str, Any] = {}
    saw_events = False
    try:
        reader.expect("{")
        first_field = True
        while reader.peek() != "}":
            if first_field:
                first_field = False
            else:
                reader.expect(",")
            key = reader.decode()
            if not isinstance(key, str):
                raise ValueError(f"non-string top-level key in {path}")
            reader.expect(":")
            if key != "events":
                metadata[key] = reader.decode()
                continue
            if saw_events:
                raise ValueError(f"duplicate top-level events field in {path}")
            saw_events = True
            reader.expect("[")
            first_event = True
            while reader.peek() != "]":
                if first_event:
                    first_event = False
                else:
                    reader.expect(",")
                event = reader.decode()
                if not isinstance(event, dict):
                    raise ValueError(f"non-object event in {path}")
                visitor(event)
            reader.expect("]")
        reader.expect("}")
        reader.ensure_finished()
    except (OSError, ValueError) as error:
        raise SystemExit(str(error)) from error
    finally:
        reader.close()
    if not saw_events:
        raise SystemExit(f"missing top-level events field in {path}")
    return metadata


def chromosome_paths(screen: Path, prefix: str) -> list[Path]:
    pattern = re.compile(rf"{re.escape(prefix)}\.(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y))\.json")
    paths: dict[str, Path] = {}
    for path in screen.iterdir():
        match = pattern.fullmatch(path.name)
        if match is None:
            continue
        chrom = match.group(1)
        if chrom in paths:
            raise SystemExit(f"duplicate {prefix} result for {chrom}")
        paths[chrom] = path
    missing = [chrom for chrom in CHROMS if chrom not in paths]
    if missing:
        raise SystemExit(f"missing {prefix} chromosome results: {', '.join(missing)}")
    return [paths[chrom] for chrom in CHROMS]


@dataclass
class ArmSummary:
    emitted_events: int
    catalogue_recurrent_events: int
    min_row_informative: int | None
    rows: list[dict]
    selected: dict[str, dict]


def group_counts(event: dict, group: str) -> dict[str, dict]:
    result = {}
    for sample in event["sample_rows"]:
        rows = {row["group"]: row for row in sample.get("group_rows", [])}
        if group in rows:
            result[sample["sample"]] = rows[group]
    return result


def tm_values(event: dict) -> dict[str, dict]:
    values = {}
    for sample in event["sample_rows"]:
        rows = {row["group"]: row for row in sample.get("group_rows", [])}
        if {"T", "M"} - rows.keys():
            continue
        t, m = rows["T"], rows["M"]
        if t["usage_fraction"] is None or m["usage_fraction"] is None:
            continue
        values[sample["sample"]] = {
            "T_usage": t["usage_fraction"],
            "M_usage": m["usage_fraction"],
            "delta_T_minus_M": t["usage_fraction"] - m["usage_fraction"],
            "T_informative": t["informative_umis"],
            "M_informative": m["informative_umis"],
            "T_usage_denominator": t["include_only"] + t["exclude_only"],
            "M_usage_denominator": m["include_only"] + m["exclude_only"],
        }
    return values


def genes(event: dict) -> str:
    return ",".join(gene["gene_name"] for gene in event["annotation"]["genes"])


def selected_summary(event: dict, grouped_tm: bool) -> dict:
    return {
        "event_id": event["id"],
        "annotation": event["annotation"],
        "samples": tm_values(event)
        if grouped_tm
        else {
            sample: {
                "usage": row["usage_fraction"],
                "informative": row["informative_umis"],
            }
            for sample, row in group_counts(event, "called").items()
        },
    }


def metadata_row_threshold(metadata: dict, path: Path) -> int | None:
    planning = metadata.get("planning", {})
    nested = planning.get("row_filter", {}) if isinstance(planning, dict) else {}
    candidates = (
        metadata.get("min_row_informative"),
        planning.get("min_row_informative") if isinstance(planning, dict) else None,
        nested.get("min_row_informative") if isinstance(nested, dict) else None,
    )
    values = {int(value) for value in candidates if value is not None and int(value) > 0}
    if len(values) > 1:
        raise SystemExit(f"conflicting min-row-informative metadata in {path}: {sorted(values)}")
    return next(iter(values), None)


def planning_candidate_count(metadata: dict, path: Path) -> int:
    planning = metadata.get("planning")
    if not isinstance(planning, dict) or not isinstance(planning.get("candidate_events"), int):
        raise SystemExit(f"missing integer planning.candidate_events in {path}")
    return planning["candidate_events"]


def summarize_pbmc(screen: Path) -> ArmSummary:
    rows: list[dict] = []
    selected: dict[str, dict] = {}
    selected_ids = {SELECTED_EVENTS[name]: name for name in ("FYB1", "CD47")}
    seen: set[str] = set()
    emitted_events = 0
    catalogue_recurrent_events = 0
    thresholds: set[int | None] = set()

    def visit(event: dict) -> None:
        nonlocal emitted_events
        event_id = event.get("id")
        if not isinstance(event_id, str):
            raise SystemExit("event is missing a string id")
        if event_id in seen:
            raise SystemExit(f"duplicate event {event_id}")
        seen.add(event_id)
        emitted_events += 1
        if event_id in selected_ids:
            selected[selected_ids[event_id]] = selected_summary(event, grouped_tm=True)

        values = tm_values(event)
        if set(values) != set(PBMC):
            return
        if any(
            min(value["T_usage_denominator"], value["M_usage_denominator"]) < 10
            for value in values.values()
        ):
            return
        deltas = [value["delta_T_minus_M"] for value in values.values()]
        if not (all(delta > 0 for delta in deltas) or all(delta < 0 for delta in deltas)):
            return
        effects = [abs(delta) for delta in deltas]
        rows.append(
            {
                "event_id": event_id,
                "event_type": event["event_type"],
                "genes": genes(event),
                "fully_annotated": event["annotation"]["fully_annotated"],
                "min_abs_delta": min(effects),
                "median_abs_delta": statistics.median(effects),
                "total_usage_denominator": sum(
                    value["T_usage_denominator"] + value["M_usage_denominator"]
                    for value in values.values()
                ),
                "samples": values,
            }
        )

    for path in chromosome_paths(screen, "pbmc4-tm"):
        metadata = visit_event_document(path, visit)
        catalogue_recurrent_events += planning_candidate_count(metadata, path)
        thresholds.add(metadata_row_threshold(metadata, path))
    if len(thresholds) != 1:
        raise SystemExit(f"inconsistent PBMC min-row-informative thresholds: {thresholds}")
    missing_selected = set(selected_ids.values()) - set(selected)
    if missing_selected:
        raise SystemExit(f"selected PBMC events absent: {', '.join(sorted(missing_selected))}")
    rows.sort(
        key=lambda row: (
            row["min_abs_delta"],
            row["median_abs_delta"],
            row["total_usage_denominator"],
        ),
        reverse=True,
    )
    return ArmSummary(
        emitted_events=emitted_events,
        catalogue_recurrent_events=catalogue_recurrent_events,
        min_row_informative=next(iter(thresholds)),
        rows=rows,
        selected=selected,
    )


def summarize_tissue(screen: Path) -> ArmSummary:
    rows: list[dict] = []
    selected: dict[str, dict] = {}
    selected_ids = {
        SELECTED_EVENTS[name]: name for name in ("FNBP1", "MYL6", "PPP1R12A", "RPS24")
    }
    seen: set[str] = set()
    emitted_events = 0
    catalogue_recurrent_events = 0
    thresholds: set[int | None] = set()

    def visit(event: dict) -> None:
        nonlocal emitted_events
        event_id = event.get("id")
        if not isinstance(event_id, str):
            raise SystemExit("event is missing a string id")
        if event_id in seen:
            raise SystemExit(f"duplicate event {event_id}")
        seen.add(event_id)
        emitted_events += 1
        if event_id in selected_ids:
            selected[selected_ids[event_id]] = selected_summary(event, grouped_tm=False)

        values = group_counts(event, "called")
        if set(values) != set(ALL):
            return
        if any(value["include_only"] + value["exclude_only"] < 20 for value in values.values()):
            return
        pbmc_usage = [values[sample]["usage_fraction"] for sample in PBMC]
        brain_usage = [values[sample]["usage_fraction"] for sample in BRAIN]
        if any(value is None for value in pbmc_usage + brain_usage):
            return
        pbmc_median = statistics.median(pbmc_usage)
        shifts = [value - pbmc_median for value in brain_usage]
        if shifts[0] * shifts[1] <= 0:
            return
        pbmc_range = max(pbmc_usage) - min(pbmc_usage)
        brain_range = max(brain_usage) - min(brain_usage)
        min_shift = min(abs(value) for value in shifts)
        if pbmc_range > 0.15 or brain_range > 0.20 or min_shift < 0.15:
            return
        rows.append(
            {
                "event_id": event_id,
                "event_type": event["event_type"],
                "genes": genes(event),
                "fully_annotated": event["annotation"]["fully_annotated"],
                "pbmc_median_usage": pbmc_median,
                "pbmc_range": pbmc_range,
                "brain_range": brain_range,
                "min_brain_abs_delta": min_shift,
                "samples": {
                    sample: {
                        "usage": values[sample]["usage_fraction"],
                        "informative": values[sample]["informative_umis"],
                        "usage_denominator": (
                            values[sample]["include_only"] + values[sample]["exclude_only"]
                        ),
                    }
                    for sample in ALL
                },
            }
        )

    for path in chromosome_paths(screen, "six-called"):
        metadata = visit_event_document(path, visit)
        catalogue_recurrent_events += planning_candidate_count(metadata, path)
        thresholds.add(metadata_row_threshold(metadata, path))
    if len(thresholds) != 1:
        raise SystemExit(f"inconsistent tissue min-row-informative thresholds: {thresholds}")
    missing_selected = set(selected_ids.values()) - set(selected)
    if missing_selected:
        raise SystemExit(f"selected tissue events absent: {', '.join(sorted(missing_selected))}")
    rows.sort(key=lambda row: (row["min_brain_abs_delta"], -row["pbmc_range"]), reverse=True)
    return ArmSummary(
        emitted_events=emitted_events,
        catalogue_recurrent_events=catalogue_recurrent_events,
        min_row_informative=next(iter(thresholds)),
        rows=rows,
        selected=selected,
    )


def event_count_fields(summary: ArmSummary) -> dict[str, int | str | None]:
    if summary.min_row_informative is None:
        # Preserve the v1 result schema and meaning exactly for unfiltered archives.
        return {"events_screened": summary.emitted_events}
    return {
        "catalogue_recurrent_events": summary.catalogue_recurrent_events,
        "denominator_qualified_events": summary.emitted_events,
        "legacy_summed_minimum_events": None,
        "legacy_summed_minimum_status": "not_evaluated_due_to_row_pushdown",
        "min_row_informative": summary.min_row_informative,
    }


def parse_matrix_market_selected(
    matrix: Path,
    row_names: dict[int, str],
    column_labels: dict[int, tuple[str, ...]],
) -> dict[tuple[str, str], int]:
    totals: dict[tuple[str, str], int] = defaultdict(int)
    with matrix.open() as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.split()
            if len(fields) != 3:
                continue
            row, column, value = map(int, fields)
            if row not in row_names or column not in column_labels:
                continue
            for label in column_labels[column]:
                totals[(row_names[row], label)] += value
    return totals


def sj_audit(root: Path, screen: Path) -> dict[str, dict]:
    audit = {}
    wanted = {coordinate: name for name, coordinate in SJ_TARGETS.items()}
    for dataset in ALL:
        sj_dir = root / SJ_DIRS[dataset]
        called = {line.split()[0] for line in (screen / f"{dataset}.called.tsv").open()}
        tm = {}
        tm_path = screen / f"{dataset}.TM.tsv"
        if tm_path.exists():
            tm = {fields[0]: fields[1] for fields in map(str.split, tm_path.open())}
        column_labels = {}
        with (sj_dir / "barcodes.tsv").open() as handle:
            for index, line in enumerate(handle, 1):
                barcode = line.strip()
                labels = []
                if barcode in called:
                    labels.append("called")
                if barcode in tm:
                    labels.append(tm[barcode])
                if labels:
                    column_labels[index] = tuple(labels)
        row_names = {}
        with (sj_dir / "features.tsv").open() as handle:
            for index, line in enumerate(handle, 1):
                fields = line.split()
                coordinate = (fields[0], int(fields[1]), int(fields[2]))
                if coordinate in wanted:
                    row_names[index] = wanted[coordinate]
        required = {
            name for name in SJ_TARGETS
            if not name.startswith("FYB1_") or dataset in PBMC
        }
        missing = required - set(row_names.values())
        if missing:
            raise SystemExit(f"{dataset} SJ matrix lacks targets {sorted(missing)}")
        totals = parse_matrix_market_selected(sj_dir / "matrix.mtx", row_names, column_labels)
        record: dict[str, object] = {"selected_called_cells": len(called)}
        for locus, skip_name, include_name in (
            ("FYB1", "FYB1_skip", "FYB1_include_right"),
            ("FNBP1", "FNBP1_skip", "FNBP1_include_right"),
            ("MYL6", "MYL6_skip", "MYL6_include_right"),
        ):
            if locus == "FYB1" and dataset not in PBMC:
                continue
            labels = ("T", "M") if locus == "FYB1" else ("called",)
            record[locus] = {}
            for label in labels:
                skip = totals[(skip_name, label)]
                include = totals[(include_name, label)]
                denominator = skip + include
                record[locus][label] = {
                    "skip_umis": skip,
                    "include_right_umis": include,
                    "skip_fraction": None if denominator == 0 else skip / denominator,
                    "include_right_fraction": None if denominator == 0 else include / denominator,
                }
        audit[dataset] = record
    return audit


def timing_summary(
    screen: Path, prefix: str, concurrent_driver: bool = False
) -> dict[str, float | int]:
    walls = []
    rss_values = []
    for chrom in CHROMS:
        path = screen / f"{prefix}.{chrom}.time.txt"
        if not path.is_file():
            raise SystemExit(f"missing timing record {path}")
        text = path.read_text()
        wall_match = re.search(r"wall_seconds=([0-9.]+)", text) or re.search(
            r"wall=([0-9.]+)", text
        )
        if wall_match is None:
            verbose = re.search(
                r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([0-9.:]+)", text
            )
            if verbose is not None:
                parts = [float(value) for value in verbose.group(1).split(":")]
                walls.append(sum(value * 60**power for power, value in enumerate(reversed(parts))))
        else:
            walls.append(float(wall_match.group(1)))
        rss_match = re.search(r"peak_rss_kib=([0-9]+)", text) or re.search(
            r"rss_kb=([0-9]+)", text
        ) or re.search(r"Maximum resident set size \(kbytes\):\s*([0-9]+)", text)
        if rss_match is not None:
            rss_values.append(int(rss_match.group(1)))
    if len(walls) != 24 or len(rss_values) != 24:
        raise SystemExit(
            f"expected 24 timing records for {prefix}, found {len(walls)} wall/{len(rss_values)} RSS"
        )
    wall_key = (
        "summed_command_wall_seconds"
        if concurrent_driver
        else "sequential_chromosome_wall_seconds"
    )
    return {wall_key: sum(walls), "maximum_peak_rss_kib": max(rss_values)}


def atomic_write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = ""
    try:
        with tempfile.NamedTemporaryFile(
            "w", newline="", dir=path.parent, prefix=f".{path.name}.", delete=False
        ) as handle:
            temporary_name = handle.name
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow({column: row[column] for column in columns})
        os.replace(temporary_name, path)
    finally:
        if temporary_name:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = ""
    try:
        with tempfile.NamedTemporaryFile(
            "w", dir=path.parent, prefix=f".{path.name}.", delete=False
        ) as handle:
            temporary_name = handle.name
            handle.write(text)
        os.replace(temporary_name, path)
    finally:
        if temporary_name:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-pbmc", type=Path, required=True)
    parser.add_argument("--out-tissue", type=Path, required=True)
    parser.add_argument("--driver-started-utc")
    parser.add_argument("--driver-start-ns", type=int)
    parser.add_argument("--driver-workers", type=int)
    parser.add_argument("--driver-total-threads", type=int)
    parser.add_argument("--driver-per-process-threads", type=int)
    args = parser.parse_args()

    driver_values = (
        args.driver_started_utc,
        args.driver_start_ns,
        args.driver_workers,
        args.driver_total_threads,
        args.driver_per_process_threads,
    )
    if any(value is not None for value in driver_values) and any(
        value is None for value in driver_values
    ):
        parser.error("all driver metadata arguments must be supplied together")
    if args.driver_workers is not None and (
        args.driver_start_ns < 1
        or args.driver_workers < 1
        or args.driver_total_threads < 1
        or args.driver_per_process_threads < 1
        or args.driver_workers * args.driver_per_process_threads > args.driver_total_threads
    ):
        parser.error("driver concurrency metadata violates the fixed total thread budget")

    pbmc = summarize_pbmc(args.screen)
    tissue = summarize_tissue(args.screen)

    atomic_write_tsv(
        args.out_pbmc,
        pbmc.rows,
        [
            "event_id",
            "event_type",
            "genes",
            "fully_annotated",
            "min_abs_delta",
            "median_abs_delta",
            "total_usage_denominator",
        ],
    )
    atomic_write_tsv(
        args.out_tissue,
        tissue.rows,
        [
            "event_id",
            "event_type",
            "genes",
            "fully_annotated",
            "pbmc_median_usage",
            "pbmc_range",
            "brain_range",
            "min_brain_abs_delta",
        ],
    )

    selected = {
        name: (pbmc.selected if name in {"FYB1", "CD47"} else tissue.selected)[name]
        for name in SELECTED_EVENTS
    }

    result = {
        "schema": "gravlax.federated-biology-screen.v1",
        "status": "EXPLORATORY_POSITIVE_CONTROLS_PASS_FNBP1_CONFIRMATION_LOCKED",
        "interpretation": {
            "screen": "post-hoc exploratory; no confirmatory p-values are assigned",
            "positive_controls": "literature and RT-PCR/RNA-FISH validation were checked after ranking",
            "tissue_screen": "descriptive and potentially confounded by cell composition, cancer, and nuclei/cell protocol",
            "SJ_audit": "same alignments but an independent STAR SJ matrix aggregation; reducer audit, not independent biological validation",
        },
        "groups": json.loads((args.screen / "groups.json").read_text()),
        "pbmc4": {
            **event_count_fields(pbmc),
            "all_four_testable_strictly_direction_consistent": len(pbmc.rows),
            "min_abs_delta_at_least_0_10": sum(
                row["min_abs_delta"] >= 0.10 for row in pbmc.rows
            ),
            "min_abs_delta_at_least_0_15": sum(
                row["min_abs_delta"] >= 0.15 for row in pbmc.rows
            ),
            "candidate_table": args.out_pbmc.name,
            "performance": timing_summary(
                args.screen, "pbmc4-tm", concurrent_driver=pbmc.min_row_informative is not None
            ),
        },
        "six_shard_tissue": {
            **event_count_fields(tissue),
            "strict_candidates": len(tissue.rows),
            "criteria": {
                "each_sample_usage_denominator_min": 20,
                "PBMC_usage_range_max": 0.15,
                "brain_usage_range_max": 0.20,
                "each_brain_abs_delta_from_PBMC_median_min": 0.15,
            },
            "candidate_table": args.out_tissue.name,
            "performance": timing_summary(
                args.screen,
                "six-called",
                concurrent_driver=tissue.min_row_informative is not None,
            ),
        },
        "selected_events": selected,
        "STAR_SJ_matrix_audit": sj_audit(args.root, args.screen),
        "prospective_confirmation": {
            "primary": "FNBP1",
            "event_id": SELECTED_EVENTS["FNBP1"],
            "locked_gate_manifest": "experiments/manifests/federated-biological-screen.yaml",
        },
    }
    if args.driver_workers is not None:
        finished_ns = time.time_ns()
        result["driver"] = {
            "started_utc": args.driver_started_utc,
            "finished_utc": datetime.datetime.now(datetime.timezone.utc).strftime(
                "%Y-%m-%dT%H:%M:%SZ"
            ),
            "wall_seconds": (finished_ns - args.driver_start_ns) / 1_000_000_000,
            "workers": args.driver_workers,
            "total_thread_budget": args.driver_total_threads,
            "threads_per_process": args.driver_per_process_threads,
            "maximum_concurrent_rayon_threads": (
                args.driver_workers * args.driver_per_process_threads
            ),
        }
    serialized = json.dumps(result, indent=2, sort_keys=True)
    atomic_write_text(args.out_json, serialized + "\n")
    print(serialized)


if __name__ == "__main__":
    main()
