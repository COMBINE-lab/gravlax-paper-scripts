#!/usr/bin/env python
"""Gates DISCR-1/2: the discover→GTF→replay loop, judged on the 45 cleanly-withheld genes
(no residual-w1 overlap, Gene-oracle >= 50 UMIs), vs the GeneFull oracle (correct comparator
per D24). Frozen: DISCR-1 PASS = 45/45 quantified >=1 UMI; DISCR-2 PASS = median rel err <=35%.
"""
import os
import bisect
import numpy as np
import scipy.io

P = os.environ["GRAVLAX_PROJECT_ROOT"]

def gene_spans(gtf, want=None):
    out = {}
    for line in open(gtf):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        i = f[8].find('gene_id "')
        g = f[8][i + 9 : f[8].index('"', i + 9)]
        if want is None or g in want:
            out[g] = (f[0], int(f[3]) - 1, int(f[4]))
    return out

cells = {l.strip() for l in open(P + "/runs/oracle/full/v49/Solo.out/Gene/filtered/barcodes.tsv")}

def load_sums(d):
    m = scipy.io.mmread(d + "/matrix.mtx").tocsr()
    genes = [l.split("\t")[0] for l in open(d + "/features.tsv")]
    bcs = [l.strip() for l in open(d + "/barcodes.tsv")]
    cols = [i for i, b in enumerate(bcs) if b in cells]
    return dict(zip(genes, np.asarray(m[:, cols].sum(axis=1)).ravel()))

wh = [l.split("\t")[1] for l in open(P + "/annotations/withheld/w1.manifest.tsv") if l.startswith("gene\t")]
spans = gene_spans(P + "/annotations/gencode.v49.annotation.gtf", set(wh))
o_gene = load_sums(P + "/runs/oracle/full/v49/Solo.out/Gene/raw")
o_full = load_sums(P + "/runs/oracle/full/v49/Solo.out/GeneFull/raw")

w1_spans = {}
for g, (ch, s, e) in gene_spans(P + "/annotations/withheld/w1.gtf").items():
    w1_spans.setdefault(ch, []).append((s, e))
for ch in w1_spans:
    w1_spans[ch].sort()

def overlaps_w1(ch, s, e):
    lst = w1_spans.get(ch, [])
    i = bisect.bisect_left(lst, (e, 0))
    return any(lst[j][0] < e and lst[j][1] > s for j in range(max(0, i - 200), min(len(lst), i + 1)))

# Candidate loci spans from the emitted GTF, with their replayed counts.
novel = load_sums(P + "/runs/replay/d0-novel")
nspans = gene_spans(P + "/runs/replay/d0-novel.gtf")
by_chrom = {}
for g, (ch, s, e) in nspans.items():
    by_chrom.setdefault(ch, []).append((s, e, g))

hit = tot = 0
errs = []
for g in wh:
    if o_gene.get(g, 0) < 50:
        continue
    ch, s, e = spans[g]
    if overlaps_w1(ch, s, e):
        continue
    tot += 1
    got = sum(novel.get(cg, 0) for (cs, ce, cg) in by_chrom.get(ch, []) if cs < e and ce > s)
    if got > 0:
        hit += 1
        oracle = o_full.get(g, 0)
        errs.append(abs(got - oracle) / max(oracle, 1))

errs = np.array(errs)
print(f"DISCR-1: {hit}/{tot} cleanly-withheld genes quantified through the loop  [PASS = {tot}/{tot}]")
print(f"DISCR-2: median rel err vs GeneFull {100*np.median(errs):.1f}% (n={len(errs)}), corr n/a  [PASS <=35%]")
