#!/usr/bin/env python
"""APA-D demonstration: differential 3' usage between cell populations, from the archive alone.

Groups: T cells (CD3E+ LYZ-) vs monocytes (LYZ+ CD3E-) from the D0 oracle matrix (annotation used
only to LABEL cells; the 3' evidence is annotation-free). Per gene: weighted mean 3'-site position
per group over sites with >=MIN UMIs in both groups; shift = |mean_A - mean_B|. Frozen null: the
same statistic between two random halves of the T-cell group; report genes whose real shift
exceeds the null's 95th percentile.
"""
import os
import random
import subprocess
from pathlib import Path

import numpy as np
import scipy.io

P = Path(os.environ["GRAVLAX_PROJECT_ROOT"])
AIE = P / "src/target/release/aie"
ARCH = P / "runs/archive/d0.aie"
MIN_UMIS = 30
random.seed(11)

d = P / "runs/oracle/full/v49/Solo.out/Gene/filtered"
m = scipy.io.mmread(str(d / "matrix.mtx")).tocsr()
genes = [l.split("\t")[1] for l in open(d / "features.tsv")]
bcs = [l.strip() for l in open(d / "barcodes.tsv")]
gi = {g: i for i, g in enumerate(genes)}
cd3e = np.asarray(m[gi["CD3E"], :].todense()).ravel()
lyz = np.asarray(m[gi["LYZ"], :].todense()).ravel()
tcells = [b for b, c, l in zip(bcs, cd3e, lyz) if c > 0 and l == 0]
monos = [b for b, c, l in zip(bcs, cd3e, lyz) if l > 3 and c == 0]
print(f"groups: {len(tcells)} T cells, {len(monos)} monocytes")

groups_real = P / "runs/replay/apa-groups-real.tsv"
groups_null = P / "runs/replay/apa-groups-null.tsv"
with open(groups_real, "w") as f:
    for b in tcells:
        f.write(f"{b}\tT\n")
    for b in monos:
        f.write(f"{b}\tM\n")
half = random.sample(tcells, len(tcells) // 2)
halfset = set(half)
with open(groups_null, "w") as f:
    for b in tcells:
        f.write(f"{b}\t{'A' if b in halfset else 'B'}\n")

# Top expressed genes with spans.
spans = {}
for line in open(P / "annotations/gencode.v49.annotation.gtf"):
    if line.startswith("#"):
        continue
    f2 = line.split("\t")
    if f2[2] != "gene":
        continue
    i = f2[8].find('gene_name "')
    if i < 0:
        continue
    name = f2[8][i + 11 : f2[8].index('"', i + 11)]
    spans.setdefault(name, (f2[0], int(f2[3]) - 1, int(f2[4]), f2[6]))
expr = np.asarray(m.sum(axis=1)).ravel()
top = [genes[i] for i in np.argsort(-expr)[:250] if genes[i] in spans]

def shifts(groups_file, ga, gb):
    out = {}
    for g in top:
        ch, s, e, st = spans[g]
        r = subprocess.run(
            [str(AIE), "query", str(ARCH), "apa", f"{ch}:{s}-{e}", "--strand", st, "--tsv",
             "--groups", str(groups_file)],
            capture_output=True, text=True)
        pos = {ga: [], gb: []}
        hdr = None
        for line in r.stdout.splitlines():
            f3 = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                hdr = f3
                continue
            mid = (int(f3[1]) + int(f3[2])) / 2
            ca, cb = int(f3[hdr.index(ga)]), int(f3[hdr.index(gb)])
            pos[ga].append((mid, ca))
            pos[gb].append((mid, cb))
        def wmean(pairs):
            tot = sum(c for _, c in pairs)
            return (sum(p * c for p, c in pairs) / tot, tot) if tot >= MIN_UMIS else (None, tot)
        (ma, ta), (mb, tb) = wmean(pos[ga]), wmean(pos[gb])
        if ma is not None and mb is not None:
            out[g] = abs(ma - mb)
    return out

real = shifts(groups_real, "T", "M")
null = shifts(groups_null, "A", "B")
common = sorted(set(real) & set(null))
nv = np.array([null[g] for g in common])
rv = np.array([real[g] for g in common])
p95 = np.percentile(nv, 95)
sig = [(g, real[g]) for g in common if real[g] > p95]
sig.sort(key=lambda x: -x[1])
print(f"genes measured in both: {len(common)}; null median shift {np.median(nv):.0f} bp (P95 {p95:.0f}); real median {np.median(rv):.0f} bp")
print(f"genes with real shift > null P95: {len(sig)} ({100*len(sig)/max(len(common),1):.0f}%)  [null expectation 5%]")
for g, s in sig[:10]:
    print(f"  {g}: {s:.0f} bp shift (null {null[g]:.0f})")
