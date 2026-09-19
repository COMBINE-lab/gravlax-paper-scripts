#!/usr/bin/env python3
"""Fail-closed equivalence check for sparse cohort-event output versus matched JSON."""

import argparse
import csv
import hashlib
import io
import json
import re
import shlex
import subprocess
from pathlib import Path


EVENT_HEADER_JSON = [
    "event_id",
    "event_type",
    "chrom",
    "inclusion_junctions",
    "exclusion_junctions",
    "annotation_genes_json",
    "strand",
    "fully_annotated",
]
EVENT_HEADER_LEGACY = [
    "event_id",
    "event_type",
    "chrom",
    "inclusion_junctions",
    "exclusion_junctions",
    "gene_ids",
    "gene_names",
    "strand",
    "fully_annotated",
]
PRESENCE_HEADER = ["event_id", "sample"]
COUNT_HEADER = [
    "event_id",
    "sample",
    "aggregation",
    "group",
    "include_only",
    "exclude_only",
    "both",
    "cells",
    "selected_cells",
]
SPARSE_FILES = {
    "events": "events.tsv.zst",
    "presence": "presence.tsv.zst",
    "counts": "counts.tsv.zst",
}


class ValidationError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValidationError(f"cannot read JSON {path}: {error}") from error


def load_zstd_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        process = subprocess.run(
            ["zstd", "-q", "-dc", str(path)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        detail = getattr(error, "stderr", b"").decode("utf-8", "replace").strip()
        raise ValidationError(f"cannot decompress {path}: {detail or error}") from error
    try:
        text = process.stdout.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValidationError(f"{path} is not UTF-8: {error}") from error
    require(text.endswith("\n"), f"{path} lacks a final newline")
    reader = csv.reader(io.StringIO(text), delimiter="\t", strict=True)
    try:
        header = next(reader)
    except StopIteration as error:
        raise ValidationError(f"{path} is empty") from error
    require(len(header) == len(set(header)), f"{path} has duplicate header fields")
    rows = []
    for line_number, values in enumerate(reader, 2):
        require(
            len(values) == len(header),
            f"{path}:{line_number}: expected {len(header)} fields, got {len(values)}",
        )
        rows.append(dict(zip(header, values)))
    return header, rows


def parse_uint(value: str, label: str) -> int:
    require(re.fullmatch(r"0|[1-9][0-9]*", value) is not None, f"{label}: invalid uint {value!r}")
    return int(value)


def parse_bool(value: str, label: str):
    require(value in {"true", "false", "NA"}, f"{label}: invalid boolean {value!r}")
    return {"true": True, "false": False, "NA": None}[value]


def comma_list(value: str) -> list[str]:
    return [] if value == "" else value.split(",")


def first_difference(actual: list, expected: list, label: str) -> None:
    if actual == expected:
        return
    limit = min(len(actual), len(expected))
    index = next((i for i in range(limit) if actual[i] != expected[i]), limit)
    left = actual[index] if index < len(actual) else "<missing>"
    right = expected[index] if index < len(expected) else "<missing>"
    raise ValidationError(
        f"{label} mismatch at row {index + 1}: actual={left!r}, expected={right!r}; "
        f"row counts {len(actual)} versus {len(expected)}"
    )


def base_counts(row: dict, label: str) -> dict[str, int]:
    counts = {key: row[key] for key in ("include_only", "exclude_only", "both")}
    require(all(type(value) is int and value >= 0 for value in counts.values()), f"{label}: bad count")
    informative = counts["include_only"] + counts["exclude_only"]
    require(row["informative_umis"] == informative, f"{label}: informative_umis is inconsistent")
    usage = None if informative == 0 else counts["include_only"] / informative
    require(row["usage_fraction"] == usage, f"{label}: usage_fraction is inconsistent")
    return counts


def validate_dense(dense: dict) -> tuple[list[dict], list[dict]]:
    expected_keys = {
        "coordinates",
        "events",
        "locus",
        "min_informative",
        "min_row_informative",
        "min_samples",
        "min_support",
        "planning",
        "schema",
        "semantics",
    }
    require(set(dense) == expected_keys, f"dense JSON top-level fields changed: {sorted(dense)}")
    require(dense["schema"] == "gravlax.cohort.events.v1", "unsupported dense JSON schema")
    events = dense["events"]
    require(type(events) is list and events, "dense JSON has no events")
    event_ids = [event.get("id") for event in events]
    require(len(event_ids) == len(set(event_ids)), "dense JSON has duplicate event IDs")
    require(dense["planning"]["retained_events"] == len(events), "dense retained_events mismatch")

    dimensions = None
    for event_index, event in enumerate(events):
        label = f"dense event {event_index + 1}"
        require(
            set(event)
            == {
                "annotation",
                "chrom",
                "event_type",
                "exclusion_junctions",
                "id",
                "inclusion_junctions",
                "sample_rows",
            },
            f"{label}: fields changed",
        )
        annotation = event["annotation"]
        require(
            set(annotation) == {"fully_annotated", "genes", "strand"},
            f"{label}: annotation fields changed",
        )
        for gene in annotation["genes"]:
            require(set(gene) == {"gene_id", "gene_name"}, f"{label}: gene fields changed")
        sample_dimensions = []
        seen_samples = set()
        for sample in event["sample_rows"]:
            sample_id = sample["sample"]
            require(sample_id not in seen_samples, f"{label}: duplicate sample {sample_id}")
            seen_samples.add(sample_id)
            require(
                set(sample)
                == {"archive", "cells", "group_rows", "present", "sample", "scope", "totals"},
                f"{label}/{sample_id}: fields changed",
            )
            require(type(sample["cells"]) is int and sample["cells"] >= 0, f"{label}/{sample_id}: bad cells")
            base_counts(sample["totals"], f"{label}/{sample_id}/total")
            scope = sample["scope"]
            require(
                set(scope) == {"aggregation", "groups", "selected_cells", "source"},
                f"{label}/{sample_id}: scope fields changed",
            )
            require(type(scope["selected_cells"]) is int and scope["selected_cells"] >= 0,
                    f"{label}/{sample_id}: bad selected_cells")
            if scope["aggregation"] == "group":
                require(len(sample["group_rows"]) == len(scope["groups"]),
                        f"{label}/{sample_id}: group dimension mismatch")
                for group_row, group_dimension in zip(sample["group_rows"], scope["groups"]):
                    require(
                        set(group_row)
                        == {
                            "both",
                            "cells",
                            "exclude_only",
                            "group",
                            "include_only",
                            "informative_umis",
                            "selected_cells",
                            "usage_fraction",
                        },
                        f"{label}/{sample_id}: group-row fields changed",
                    )
                    require(group_row["group"] == group_dimension["name"],
                            f"{label}/{sample_id}: group name mismatch")
                    require(group_row["selected_cells"] == group_dimension["selected_cells"],
                            f"{label}/{sample_id}: group selected_cells mismatch")
                    require(type(group_row["cells"]) is int and group_row["cells"] >= 0,
                            f"{label}/{sample_id}: bad group cells")
                    base_counts(group_row, f"{label}/{sample_id}/{group_row['group']}")
            else:
                require(scope["aggregation"] == "bulk", f"{label}/{sample_id}: bad aggregation")
                require(sample["group_rows"] == [], f"{label}/{sample_id}: bulk row has groups")
            sample_dimensions.append(
                {"sample": sample_id, "archive": sample["archive"], "scope": scope}
            )
        if dimensions is None:
            dimensions = sample_dimensions
        else:
            require(sample_dimensions == dimensions, f"{label}: sample dimensions changed")
    return events, dimensions


def event_rows(events: list[dict], header: list[str], rows: list[dict]) -> tuple[list[dict], str]:
    require(header in (EVENT_HEADER_JSON, EVENT_HEADER_LEGACY), f"unsupported events header: {header}")
    layout = "annotation_genes_json" if header == EVENT_HEADER_JSON else "legacy_gene_id_name_lists"
    actual = []
    for index, row in enumerate(rows, 2):
        label = f"events.tsv.zst:{index}"
        if layout == "annotation_genes_json":
            try:
                genes = json.loads(row["annotation_genes_json"])
            except json.JSONDecodeError as error:
                raise ValidationError(f"{label}: invalid annotation_genes_json: {error}") from error
            require(type(genes) is list, f"{label}: annotation_genes_json is not a list")
        else:
            gene_ids = comma_list(row["gene_ids"])
            gene_names = comma_list(row["gene_names"])
            require(len(gene_ids) == len(gene_names), f"{label}: gene ID/name counts differ")
            genes = [
                {"gene_id": gene_id, "gene_name": gene_name}
                for gene_id, gene_name in zip(gene_ids, gene_names)
            ]
        actual.append(
            {
                "id": row["event_id"],
                "event_type": row["event_type"],
                "chrom": row["chrom"],
                "inclusion_junctions": comma_list(row["inclusion_junctions"]),
                "exclusion_junctions": comma_list(row["exclusion_junctions"]),
                "annotation": {
                    "genes": genes,
                    "strand": None if row["strand"] == "NA" else row["strand"],
                    "fully_annotated": parse_bool(row["fully_annotated"], label),
                },
            }
        )
    expected = [
        {
            key: event[key]
            for key in (
                "id",
                "event_type",
                "chrom",
                "inclusion_junctions",
                "exclusion_junctions",
                "annotation",
            )
        }
        for event in events
    ]
    first_difference(actual, expected, "event catalogue")
    return actual, layout


def expected_presence(events: list[dict]) -> list[dict[str, str]]:
    return [
        {"event_id": event["id"], "sample": sample["sample"]}
        for event in events
        for sample in event["sample_rows"]
        if sample["present"]
    ]


def count_fact(event_id: str, sample: str, aggregation: str, group: str,
               counts: dict, cells: int, selected_cells: int) -> dict:
    return {
        "event_id": event_id,
        "sample": sample,
        "aggregation": aggregation,
        "group": group,
        "include_only": counts["include_only"],
        "exclude_only": counts["exclude_only"],
        "both": counts["both"],
        "cells": cells,
        "selected_cells": selected_cells,
    }


def expected_counts(events: list[dict]) -> list[dict]:
    facts = []
    for event in events:
        for sample in event["sample_rows"]:
            scope = sample["scope"]
            if scope["aggregation"] == "group":
                candidates = [
                    count_fact(
                        event["id"], sample["sample"], "total", "total", sample["totals"],
                        sample["cells"], scope["selected_cells"],
                    )
                ]
                candidates.extend(
                    count_fact(
                        event["id"], sample["sample"], "group", row["group"], row,
                        row["cells"], row["selected_cells"],
                    )
                    for row in sample["group_rows"]
                )
            else:
                candidates = [
                    count_fact(
                        event["id"], sample["sample"], "bulk", "bulk", sample["totals"],
                        sample["cells"], scope["selected_cells"],
                    )
                ]
            facts.extend(
                fact
                for fact in candidates
                if any(fact[key] != 0 for key in ("include_only", "exclude_only", "both", "cells"))
            )
    return facts


def parsed_count_rows(rows: list[dict[str, str]]) -> list[dict]:
    keys = ("include_only", "exclude_only", "both", "cells", "selected_cells")
    output = []
    for index, row in enumerate(rows, 2):
        parsed = {key: row[key] for key in ("event_id", "sample", "aggregation", "group")}
        parsed.update({key: parse_uint(row[key], f"counts.tsv.zst:{index}/{key}") for key in keys})
        require(
            any(parsed[key] != 0 for key in ("include_only", "exclude_only", "both", "cells")),
            f"counts.tsv.zst:{index}: explicit all-zero fact violates sparse semantics",
        )
        output.append(parsed)
    return output


def validate_metadata(metadata: dict, dense: dict, dimensions: list[dict],
                      sparse_dir: Path, row_counts: dict[str, int]) -> int:
    expected_keys = {
        "coordinates",
        "dimensions",
        "locus",
        "min_informative",
        "min_row_informative",
        "min_samples",
        "min_support",
        "output",
        "planning",
        "schema",
        "semantics",
        "source_schema",
    }
    require(set(metadata) == expected_keys, f"sparse metadata fields changed: {sorted(metadata)}")
    require(metadata["schema"] == "gravlax.cohort.events.sparse.v1", "unsupported sparse schema")
    require(metadata["source_schema"] == dense["schema"], "source_schema mismatch")
    for key in (
        "coordinates",
        "locus",
        "min_informative",
        "min_row_informative",
        "min_samples",
        "min_support",
        "planning",
    ):
        require(metadata[key] == dense[key], f"metadata {key} differs from dense JSON")
    require(
        metadata["dimensions"] == {"events": len(dense["events"]), "samples": dimensions},
        "sparse dimensions do not reproduce dense sample/event dimensions",
    )
    semantics = metadata["semantics"]
    require(semantics["candidate_union"] == dense["semantics"]["candidate_union"],
            "candidate_union semantics differ")
    require(semantics["statistics"] == dense["semantics"]["statistics"],
            "statistics semantics differ")
    require(semantics["both_in_usage_denominator"] is False,
            "both-in-denominator semantics changed")
    require(
        semantics["missing_presence_row"]
        == "catalogue present=false; evidence is not imputed",
        "missing-presence semantics changed",
    )
    require(
        semantics["missing_count_row"]
        == "exact logical zero for include_only, exclude_only, both, and cells over the explicit event/sample/scope dimensions",
        "missing-count semantics changed",
    )
    require(
        semantics["derived_fields"]
        == {
            "informative_umis": "include_only + exclude_only",
            "usage_fraction": "include_only / informative_umis, or null when informative_umis is zero",
        },
        "derived-field semantics changed",
    )
    output = metadata["output"]
    require(output["files"] == SPARSE_FILES, "sparse file manifest changed")
    require(output["rows"] == row_counts, "metadata row counts differ from decoded tables")
    compressed_bytes = sum((sparse_dir / name).stat().st_size for name in SPARSE_FILES.values())
    require(output["compressed_table_bytes"] == compressed_bytes,
            "compressed_table_bytes differs from file sizes")
    return compressed_bytes


def read_time(path: Path) -> dict:
    text = path.read_text()
    wall_match = re.search(r"Elapsed \(wall clock\) time.*: ([0-9:.]+)", text)
    rss_match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    status_match = re.search(r"Exit status: (\d+)", text)
    command_match = re.search(r'^\s*Command being timed: "(.*)"$', text, re.MULTILINE)
    require(all((wall_match, rss_match, status_match, command_match)), f"incomplete GNU time record: {path}")
    fields = [float(field) for field in wall_match.group(1).split(":")]
    require(len(fields) <= 3, f"invalid wall time in {path}")
    wall = fields[-1] + (60 * fields[-2] if len(fields) >= 2 else 0) + (
        3600 * fields[-3] if len(fields) == 3 else 0
    )
    return {
        "command": command_match.group(1),
        "wall_seconds": wall,
        "max_rss_kib": int(rss_match.group(1)),
        "exit_status": int(status_match.group(1)),
        "sha256": sha256(path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dense-json", type=Path, required=True)
    parser.add_argument("--sparse-dir", type=Path, required=True)
    parser.add_argument("--dense-time", type=Path, required=True)
    parser.add_argument("--sparse-time", type=Path, required=True)
    parser.add_argument("--metadata-stdout", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--gravlax-commit", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    require(
        re.fullmatch(r"[0-9a-f]{40}", args.gravlax_commit) is not None,
        "--gravlax-commit must be a full lowercase 40-hex commit",
    )

    dense = load_json(args.dense_json)
    events, dimensions = validate_dense(dense)
    metadata_path = args.sparse_dir / "metadata.json"
    metadata = load_json(metadata_path)
    require(metadata == load_json(args.metadata_stdout), "stdout metadata differs from metadata.json")

    event_header, event_table = load_zstd_tsv(args.sparse_dir / SPARSE_FILES["events"])
    presence_header, presence_table = load_zstd_tsv(args.sparse_dir / SPARSE_FILES["presence"])
    count_header, count_table = load_zstd_tsv(args.sparse_dir / SPARSE_FILES["counts"])
    require(presence_header == PRESENCE_HEADER, f"unsupported presence header: {presence_header}")
    require(count_header == COUNT_HEADER, f"unsupported counts header: {count_header}")
    _, event_layout = event_rows(events, event_header, event_table)
    first_difference(presence_table, expected_presence(events), "catalogue presence")
    parsed_counts = parsed_count_rows(count_table)
    first_difference(parsed_counts, expected_counts(events), "nonzero count facts")

    row_counts = {
        "events": len(event_table),
        "presence": len(presence_table),
        "nonzero_counts": len(count_table),
    }
    compressed_bytes = validate_metadata(metadata, dense, dimensions, args.sparse_dir, row_counts)
    dense_time = read_time(args.dense_time)
    sparse_time = read_time(args.sparse_time)
    require(dense_time["exit_status"] == sparse_time["exit_status"] == 0,
            "one or more matched commands failed")
    dense_tokens = shlex.split(dense_time["command"])
    sparse_tokens = shlex.split(sparse_time["command"])
    require(dense_tokens and sparse_tokens, "one or more timing commands are empty")
    require(
        Path(dense_tokens[0]).resolve() == args.binary.resolve()
        and Path(sparse_tokens[0]).resolve() == args.binary.resolve(),
        "timing command executable differs from --binary",
    )
    require(dense_tokens[-1:] == ["--json"], "dense timing command is not the JSON arm")
    require(
        len(sparse_tokens) >= 2 and sparse_tokens[-2] == "--sparse-dir",
        "sparse timing command is not the sparse arm",
    )
    require(
        Path(sparse_tokens[-1]).resolve() == args.sparse_dir.resolve(),
        "timed sparse output path differs from --sparse-dir",
    )
    require(
        dense_tokens[:-1] == sparse_tokens[:-2],
        "dense and sparse commands differ before output selection",
    )

    metadata_bytes = metadata_path.stat().st_size
    sparse_bundle_bytes = compressed_bytes + metadata_bytes
    dense_bytes = args.dense_json.stat().st_size
    result = {
        "schema": "gravlax.cohort.events.sparse-validation.v1",
        "status": "PASS",
        "gates": {
            "event_catalogue_exact": True,
            "catalogue_presence_exact": True,
            "nonzero_count_facts_exact": True,
            "logical_zeros_reconstruct_dense_rows": True,
            "sample_and_scope_dimensions_exact": True,
            "metadata_exact": True,
            "matched_commands_exit_zero": True,
        },
        "event_table_layout": event_layout,
        "dimensions": {
            "events": len(events),
            "samples": len(dimensions),
            "event_sample_rows_reconstructed": len(events) * len(dimensions),
            "presence_rows": len(presence_table),
            "nonzero_count_rows": len(count_table),
        },
        "storage": {
            "dense_json_bytes": dense_bytes,
            "sparse_compressed_table_bytes": compressed_bytes,
            "sparse_metadata_bytes": metadata_bytes,
            "sparse_bundle_bytes": sparse_bundle_bytes,
            "dense_to_sparse_bundle_ratio": dense_bytes / sparse_bundle_bytes,
            "sparse_bundle_fraction_of_dense": sparse_bundle_bytes / dense_bytes,
        },
        "runtime": {
            "dense_json": dense_time,
            "sparse": sparse_time,
            "sparse_to_dense_wall_ratio": sparse_time["wall_seconds"] / dense_time["wall_seconds"],
            "sparse_to_dense_rss_ratio": sparse_time["max_rss_kib"] / dense_time["max_rss_kib"],
        },
        "digests": {
            "gravlax_commit": args.gravlax_commit,
            "dense_json_sha256": sha256(args.dense_json),
            "metadata_sha256": sha256(metadata_path),
            "metadata_stdout_sha256": sha256(args.metadata_stdout),
            "events_sha256": sha256(args.sparse_dir / SPARSE_FILES["events"]),
            "presence_sha256": sha256(args.sparse_dir / SPARSE_FILES["presence"]),
            "counts_sha256": sha256(args.sparse_dir / SPARSE_FILES["counts"]),
            "binary_sha256": sha256(args.binary),
        },
        "inputs": {
            "dense_json": str(args.dense_json),
            "sparse_dir": str(args.sparse_dir),
            "metadata_stdout": str(args.metadata_stdout),
            "binary": str(args.binary),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except ValidationError as error:
        raise SystemExit(f"sparse cohort validation failed: {error}") from error
