#!/usr/bin/env python3
"""Compare two coordinate MatrixMarket EM layers, including explicit-zero differences."""

import argparse
import json
from pathlib import Path


def read_matrix(path: Path) -> tuple[tuple[int, int, int], dict[tuple[int, int], float]]:
    dims = None
    values: dict[tuple[int, int], float] = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("%") or not line.strip():
                continue
            fields = line.split()
            if dims is None:
                dims = tuple(map(int, fields))
                continue
            row, col, value = int(fields[0]), int(fields[1]), float(fields[2])
            values[(row, col)] = value
    if dims is None:
        raise ValueError(f"missing MatrixMarket dimensions: {path}")
    return dims, values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("packed", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    packed_dims, packed = read_matrix(args.packed)
    ref_dims, reference = read_matrix(args.reference)
    if packed_dims[:2] != ref_dims[:2]:
        raise SystemExit(f"matrix shapes differ: {packed_dims[:2]} vs {ref_dims[:2]}")
    keys = set(packed) | set(reference)
    absolute = [abs(packed.get(key, 0.0) - reference.get(key, 0.0)) for key in keys]
    ref_mass = sum(abs(value) for value in reference.values())
    extra = set(packed) - set(reference)
    missing = set(reference) - set(packed)
    result = {
        "packed_declared_nnz": packed_dims[2],
        "reference_declared_nnz": ref_dims[2],
        "coordinate_symmetric_difference": len(extra) + len(missing),
        "extra_coordinates": len(extra),
        "missing_coordinates": len(missing),
        "extra_nonzero_coordinates": sum(packed[key] != 0.0 for key in extra),
        "missing_nonzero_coordinates": sum(reference[key] != 0.0 for key in missing),
        "changed_values_on_union": sum(
            packed.get(key, 0.0) != reference.get(key, 0.0) for key in keys
        ),
        "max_absolute_error": max(absolute, default=0.0),
        "relative_l1": sum(absolute) / ref_mass if ref_mass else 0.0,
        "moved_mass": sum(absolute) / 2.0,
        "packed_mass": sum(packed.values()),
        "reference_mass": sum(reference.values()),
        "scientifically_identical_at_emitted_precision": all(value == 0.0 for value in absolute),
        "byte_identical": args.packed.read_bytes() == args.reference.read_bytes(),
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.out:
        args.out.write_text(rendered)
    else:
        print(rendered, end="")
    if not result["scientifically_identical_at_emitted_precision"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

