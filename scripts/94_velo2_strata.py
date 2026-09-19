#!/usr/bin/env python
"""VELO-2: attribute the residual velocity disagreement to the two-representative substitution.

Hypothesis: (cell,gene) entries whose counted classes all carry COMPLETE read evidence (no chain
with 2 reps standing in for >2 reads) should agree with the oracle far better than incomplete
entries. Entry flags come from replay's entries.complete.tsv.
"""
import sys
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse

oracle, replay, cells_f = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
cells = {l.strip() for l in open(cells_f)}

# Completeness mask keyed by (gene_row, cell_col) in replay coordinates.
flags = {}
for line in open(replay / "entries.complete.tsv"):
    g, c, f = line.split()
    flags[(int(g) - 1, int(c) - 1)] = f == "1"

def load(d, name):
    m = scipy.io.mmread(str(d / name)).tocsr()
    bcs = [l.strip() for l in open(d / "barcodes.tsv")]
    return m, bcs

for comp in ["spliced.mtx", "unspliced.mtx", "ambiguous.mtx"]:
    A, bca = load(oracle, comp)
    B, bcb = load(replay, comp)
    assert bca == bcb
    keep = [i for i, b in enumerate(bca) if b in cells]
    keepset = set(keep)
    Ac, Bc = A.tocoo(), B.tocoo()
    # Union of entries over filtered cells, split by replay-side completeness (entries absent
    # from replay count as incomplete=False? — attribute them separately as 'replay-missing').
    av = {}
    for g, c, v in zip(Ac.row, Ac.col, Ac.data):
        if c in keepset and v:
            av[(g, c)] = av.get((g, c), 0) + v
    bv = {}
    for g, c, v in zip(Bc.row, Bc.col, Bc.data):
        if c in keepset and v:
            bv[(g, c)] = bv.get((g, c), 0) + v
    keys = set(av) | set(bv)
    stats = {True: [0, 0], False: [0, 0], None: [0, 0]}  # flag -> [l1, oracle_mass]
    for k in keys:
        a, b = av.get(k, 0), bv.get(k, 0)
        f = flags.get(k)
        stats[f][0] += abs(a - b)
        stats[f][1] += a
    for label, f in [("complete", True), ("incomplete", False), ("replay-absent", None)]:
        l1, mass = stats[f]
        if mass or l1:
            print(f"{comp:15s} {label:14s}: oracle mass {int(mass):>9d}  L1 {int(l1):>8d}  rel {100*l1/max(mass,1):6.3f}%")
