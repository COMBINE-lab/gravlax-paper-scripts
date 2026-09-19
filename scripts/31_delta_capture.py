#!/usr/bin/env python
"""Gate -1B's decisive metric: does replay capture the CHANGE between two annotations?

Overall agreement can hide failure on exactly the molecules the capability exists for (source plan
§11, hard gate C). So compare per-gene deltas: delta_oracle = oracle(B) - oracle(A) against
delta_replay = replay(B) - replay(A), restricted to the annotation-sensitive genes (oracle count
changes by more than --sens-frac between A and B). All sums over the filtered cell set; genes
joined unversioned; genes absent from an annotation count as zero there.

Usage: 31_delta_capture.py <oracleA_raw> <oracleB_raw> <replayA_dir> <replayB_dir> \
           --cells filtered_barcodes.tsv [--sens-frac 0.10] [--min-count 20] [--out json]
"""
import argparse
import json
from pathlib import Path

import numpy as np
import scipy.io


def load_sums(d: Path, cells: set):
    m = scipy.io.mmread(str(d / "matrix.mtx")).tocsr()
    genes = [l.split("\t")[0].split(".", 1)[0] for l in (d / "features.tsv").read_text().splitlines() if l]
    bcs = [l.strip() for l in (d / "barcodes.tsv").read_text().splitlines() if l.strip()]
    cols = [i for i, b in enumerate(bcs) if b in cells]
    s = np.asarray(m[:, cols].sum(axis=1)).ravel()
    return dict(zip(genes, s))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("oracle_a", type=Path)
    ap.add_argument("oracle_b", type=Path)
    ap.add_argument("replay_a", type=Path)
    ap.add_argument("replay_b", type=Path)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--sens-frac", type=float, default=0.10)
    ap.add_argument("--min-count", type=float, default=20)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    cells = {l.strip() for l in args.cells.read_text().splitlines() if l.strip()}
    oa, ob = load_sums(args.oracle_a, cells), load_sums(args.oracle_b, cells)
    ra, rb = load_sums(args.replay_a, cells), load_sums(args.replay_b, cells)

    genes = sorted(set(oa) | set(ob))
    sens = []
    for g in genes:
        a, b = oa.get(g, 0.0), ob.get(g, 0.0)
        if max(a, b) < args.min_count:
            continue
        if abs(b - a) > args.sens_frac * max(a, 1.0):
            sens.append(g)

    d_o = np.array([ob.get(g, 0.0) - oa.get(g, 0.0) for g in sens])
    d_r = np.array([rb.get(g, 0.0) - ra.get(g, 0.0) for g in sens])

    sign_ok = int(np.sum(np.sign(d_o) == np.sign(d_r)))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.abs(d_r - d_o) / np.maximum(np.abs(d_o), 1.0)
    within10 = int(np.sum(rel <= 0.10))
    within25 = int(np.sum(rel <= 0.25))
    l1_ratio = float(np.abs(d_r - d_o).sum() / max(np.abs(d_o).sum(), 1.0))
    # Pearson on deltas reported as a supplement only; the doc forbids leaning on correlation.
    corr = float(np.corrcoef(d_o, d_r)[0, 1]) if len(sens) > 2 else float("nan")

    res = {
        "sensitive_genes": len(sens),
        "delta_mass_oracle": float(np.abs(d_o).sum()),
        "delta_mass_replay": float(np.abs(d_r).sum()),
        "sign_agreement": sign_ok, "sign_agreement_frac": sign_ok / max(len(sens), 1),
        "delta_within_10pct": within10, "delta_within_10pct_frac": within10 / max(len(sens), 1),
        "delta_within_25pct": within25, "delta_within_25pct_frac": within25 / max(len(sens), 1),
        "delta_error_l1_ratio": l1_ratio,
        "delta_pearson_supplementary": corr,
    }
    print(f"annotation-sensitive genes (>{args.sens_frac:.0%} change, >={args.min_count:g} UMIs): {len(sens)}")
    print(f"oracle delta mass: {np.abs(d_o).sum():.0f}   replay delta mass: {np.abs(d_r).sum():.0f}")
    print(f"sign agreement:        {sign_ok}/{len(sens)}  ({100*sign_ok/max(len(sens),1):.2f}%)")
    print(f"delta within 10%:      {within10}/{len(sens)}  ({100*within10/max(len(sens),1):.2f}%)")
    print(f"delta within 25%:      {within25}/{len(sens)}  ({100*within25/max(len(sens),1):.2f}%)")
    print(f"delta L1 error ratio:  {100*l1_ratio:.3f}%")
    print(f"(supplementary Pearson on deltas: {corr:.5f})")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
