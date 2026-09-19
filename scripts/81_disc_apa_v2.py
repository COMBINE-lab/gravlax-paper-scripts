#!/usr/bin/env python
"""Corrected DISC/APA gate analyses.

DISC-2 correction: discover counts every molecule in a locus (introns included — under w1 the gene
does not exist, so its pre-mRNA is all unclaimed), while the Gene oracle counts exonic-concordant
molecules only. The like-for-like oracle is **GeneFull**. The first analysis compared against Gene
and read 76% median error; that was a semantics mismatch, not a capability failure.

DISC-1 diagnosis: recall misses are expected where a withheld gene overlaps a gene still present
in w1 — discover claims generously (either strand, any transcript span), so such molecules are
never 'unclaimed'. Report recall stratified by overlap status.

APA-2 refined (rerun on a stable archive after the binary/format race): per-gene correlation of
in-window 3'-site mass against the per-gene w1-oracle deficit.
"""
import os
import bisect
import subprocess
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
            out[g] = (f[0], int(f[3]) - 1, int(f[4]), f[6])
    return out

def load_sums(d, cells):
    m = scipy.io.mmread(d + "/matrix.mtx").tocsr()
    genes = [l.split("\t")[0] for l in open(d + "/features.tsv")]
    bcs = [l.strip() for l in open(d + "/barcodes.tsv")]
    cols = [i for i, b in enumerate(bcs) if b in cells]
    return dict(zip(genes, np.asarray(m[:, cols].sum(axis=1)).ravel()))

cells = {l.strip() for l in open(P + "/runs/oracle/full/v49/Solo.out/Gene/filtered/barcodes.tsv")}
wh = [l.split("\t")[1] for l in open(P + "/annotations/withheld/w1.manifest.tsv") if l.startswith("gene\t")]
utr = [l.split("\t")[1] for l in open(P + "/annotations/withheld/w1.manifest.tsv") if l.startswith("utr3\t")]
spans = gene_spans(P + "/annotations/gencode.v49.annotation.gtf", set(wh) | set(utr))
o_gene = load_sums(P + "/runs/oracle/full/v49/Solo.out/Gene/raw", cells)
o_full = load_sums(P + "/runs/oracle/full/v49/Solo.out/GeneFull/raw", cells)
w1_gene = load_sums(P + "/runs/oracle/full/w1/Solo.out/Gene/raw", cells)

# w1 spans for overlap diagnosis.
w1_spans_by_chrom = {}
for g, (ch, s, e, st) in gene_spans(P + "/annotations/withheld/w1.gtf").items():
    w1_spans_by_chrom.setdefault(ch, []).append((s, e))
for ch in w1_spans_by_chrom:
    w1_spans_by_chrom[ch].sort()

def overlaps_w1(ch, s, e):
    lst = w1_spans_by_chrom.get(ch, [])
    i = bisect.bisect_left(lst, (e, 0))
    return any(lst[j][0] < e and lst[j][1] > s for j in range(max(0, i - 200), min(len(lst), i + 1)))

cands = []
for l in open(P + "/results/discover-w1.tsv"):
    ch, s, e, st, u, c = l.rstrip("\n").split("\t")
    cands.append((ch, int(s), int(e), st, int(u), int(c)))
by_chrom = {}
for c in cands:
    by_chrom.setdefault(c[0], []).append(c)

print("=== DISC-1 stratified recall (withheld genes, Gene-oracle >= 50 UMIs) ===")
free_hit = free_tot = ov_hit = ov_tot = 0
quant = []
for g in wh:
    if o_gene.get(g, 0) < 50:
        continue
    ch, s, e, st = spans[g]
    got = sum(u for (cch, cs, ce, cst, u, cc) in by_chrom.get(ch, []) if cs < e and ce > s)
    if overlaps_w1(ch, s, e):
        ov_tot += 1
        ov_hit += got > 0
    else:
        free_tot += 1
        free_hit += got > 0
        if got > 0:
            quant.append((o_full.get(g, 0), got))
print(f"  not overlapping any w1 gene: {free_hit}/{free_tot} ({100*free_hit/max(free_tot,1):.1f}%)  [the fair recall]")
print(f"  overlapping a w1 gene:       {ov_hit}/{ov_tot} ({100*ov_hit/max(ov_tot,1):.1f}%)  [claimed-by-design misses live here]")

q = np.array(quant, float)
rel = np.abs(q[:, 1] - q[:, 0]) / np.maximum(q[:, 0], 1)
r = np.corrcoef(q[:, 0], q[:, 1])[0, 1]
print(f"=== DISC-2 vs GeneFull (like-for-like): median rel err {100*np.median(rel):.1f}%, corr {r:.4f} (n={len(q)}) ===")

print("=== APA-2 refined (stable archive) ===")
defs, masses = [], []
for g in utr:
    ch, s, e, st = spans[g]
    win = (e - 500, e) if st == "+" else (s, s + 500)
    out = subprocess.run(
        [P + "/src/target/release/aie", "query", P + "/runs/archive/d0.aie", "apa",
         f"{ch}:{s}-{e}", "--strand", st, "--tsv"],
        capture_output=True, text=True).stdout
    m = 0
    for row in out.splitlines():
        _, ss, ee, _, u, _ = row.split("\t")
        if int(ee) >= win[0] and int(ss) < win[1]:
            m += int(u)
    defs.append(max(0.0, o_gene.get(g, 0) - w1_gene.get(g, 0)))
    masses.append(m)
defs, masses = np.array(defs), np.array(masses)
keep = defs >= 5
r2 = np.corrcoef(defs[keep], masses[keep])[0, 1]
print(f"  per-gene corr(in-window site mass, w1 deficit) = {r2:.3f} over {int(keep.sum())} genes (deficit>=5)")
print(f"  genes with deficit>=5 and in-window mass>0: {int(((defs>=5)&(masses>0)).sum())}/{int(keep.sum())}")
