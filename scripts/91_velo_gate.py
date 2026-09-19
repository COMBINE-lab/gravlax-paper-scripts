#!/usr/bin/env python
"""Gate VELO-1: velocity replay vs fresh STARsolo Velocyto, per component, filtered cells.

Frozen (docs/decisions/2026-08-30-capability-gates.md): PASS per-component mass-weighted rel. L1
<=1.5% and total mass within 5%; MARGINAL <=5%; STOP >10%.

Usage: 91_velo_gate.py <oracle_velocyto_raw> <replay_dir> <filtered_barcodes.tsv>
"""
import sys
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse

oracle, replay, cells_f = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
cells = {l.strip() for l in open(cells_f)}

def load(d, name):
    m = scipy.io.mmread(str(d / name)).tocsr()
    bcs = [l.strip() for l in open(d / "barcodes.tsv")]
    genes = [l.split("\t")[0] for l in open(d / "features.tsv")]
    cols = [i for i, b in enumerate(bcs) if b in cells]
    return m[:, cols], genes

verdicts = []
for comp in ["spliced.mtx", "unspliced.mtx", "ambiguous.mtx"]:
    A, ga = load(oracle, comp)
    B, gb = load(replay, comp)
    assert ga == gb, "gene universes differ"
    ta, tb = A.sum(), B.sum()
    l1 = abs(A - B).sum()
    rel = l1 / max(ta, 1)
    mass = abs(tb - ta) / max(ta, 1)
    v = "PASS" if rel <= 0.015 and mass <= 0.05 else ("MARGINAL" if rel <= 0.05 else "STOP" if rel > 0.10 else "MARGINAL")
    verdicts.append(v)
    print(f"{comp:15s} oracle {int(ta):>9d}  replay {int(tb):>9d}  L1 {int(l1):>8d}  rel {100*rel:6.3f}%  mass-delta {100*mass:6.3f}%  [{v}]")
print("VELO-1:", "PASS" if all(v == "PASS" for v in verdicts) else ("STOP" if "STOP" in verdicts else "MARGINAL"))
