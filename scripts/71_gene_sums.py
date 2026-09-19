#!/usr/bin/env python
"""Per-gene UMI sums from a STARsolo raw matrix, for the withholding sampler.

Usage: 71_gene_sums.py Solo.out/Gene/raw > gene_sums.tsv
Gene ids are unversioned to match scripts/70_withhold_annotation.py --expr.
"""
import sys
from pathlib import Path

import numpy as np
import scipy.io

raw = Path(sys.argv[1])
m = scipy.io.mmread(str(raw / "matrix.mtx")).tocsr()
genes = [l.split("\t")[0].split(".", 1)[0] for l in (raw / "features.tsv").read_text().splitlines() if l]
sums = np.asarray(m.sum(axis=1)).ravel()
for g, s in zip(genes, sums):
    if s > 0:
        print(f"{g}\t{s:.0f}")
