#!/usr/bin/env python
"""Gate EM-0: aie em --star vs STARsolo --soloMultiMappers EM (UniqueAndMult-EM), filtered cells."""
import sys
from pathlib import Path
import scipy.io

oracle, replay, cells_f = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
cells = {l.strip() for l in open(cells_f)}

def load(d, name):
    m = scipy.io.mmread(str(d / name)).tocsr()
    bcs = [l.strip() for l in open(d / "barcodes.tsv")]
    cols = [i for i, b in enumerate(bcs) if b in cells]
    return m[:, cols]

A = load(oracle, "UniqueAndMult-EM.mtx")
B = load(replay, "UniqueAndMult-EM.mtx")
ta, tb = A.sum(), B.sum()
l1 = abs(A - B).sum()
rel = l1 / max(ta, 1)
v = "PASS" if rel <= 0.01 else ("MARGINAL" if rel <= 0.03 else ("STOP" if rel > 0.10 else "MARGINAL"))
print(f"EM-0: oracle {ta:.1f}  replay {tb:.1f}  L1 {l1:.1f}  rel {100*rel:.3f}%  [{v}]")
