#!/usr/bin/env python
"""Pool-style reference optimization as a REPLAY (Pool et al. 2023, Nat Methods 20:1506):
extend every gene's 3' terminal region downstream (up to EXT bp, clipped before the next
annotated gene on the same strand), then replay the index under the optimized annotation and
measure recovered expression — no reads touched."""
import os
import sys
EXT = 3000
P = os.environ["GRAVLAX_PROJECT_ROOT"]
src = P + "/annotations/gencode.v49.annotation.gtf"
out = P + "/annotations/gencode.v49.ext3p.gtf"

# gene spans per chrom for clipping
genes = {}
rows = []
for line in open(src):
    if line.startswith("#"):
        rows.append(line); continue
    f = line.rstrip("\n").split("\t")
    rows.append(f)
    if f[2] == "gene":
        genes.setdefault(f[0], []).append((int(f[3]) - 1, int(f[4]), f[6]))
starts = {ch: sorted(s for s, e, st in v) for ch, v in genes.items()}
ends = {ch: sorted(e for s, e, st in v) for ch, v in genes.items()}
import bisect
def clip_plus(ch, gene_end):
    lst = starts.get(ch, [])
    i = bisect.bisect_right(lst, gene_end)
    nxt = lst[i] if i < len(lst) else gene_end + EXT + 1
    return min(gene_end + EXT, max(gene_end, nxt - 1))
def clip_minus(ch, gene_start):
    lst = ends.get(ch, [])
    i = bisect.bisect_left(lst, gene_start)
    prv = lst[i - 1] if i > 0 else gene_start - EXT - 1
    return max(gene_start - EXT, min(gene_start, prv + 1))

# First pass: per gene_id, its span & strand; per transcript, terminal exon line index.
gid_of = lambda a: a[a.find('gene_id "') + 9 : a.index('"', a.find('gene_id "') + 9)]
gene_new_end = {}
for f in rows:
    if isinstance(f, str) or f[2] != "gene":
        continue
    g = gid_of(f[8])
    if f[6] == "+":
        gene_new_end[g] = ("+", clip_plus(f[0], int(f[4])))
    else:
        gene_new_end[g] = ("-", clip_minus(f[0], int(f[3]) - 1))

# Second pass: for each transcript, find its 3'-most exon; extend exon/transcript/gene ends.
# Group feature lines by transcript to find terminal exons.
tx_exons = {}
for idx, f in enumerate(rows):
    if isinstance(f, str) or f[2] != "exon":
        continue
    a = f[8]
    t = a[a.find('transcript_id "') + 15 : a.index('"', a.find('transcript_id "') + 15)]
    tx_exons.setdefault(t, []).append(idx)
n_ext = 0
for t, idxs in tx_exons.items():
    f0 = rows[idxs[0]]
    g = gid_of(f0[8])
    if g not in gene_new_end:
        continue
    strand, newend = gene_new_end[g]
    if strand == "+":
        term = max(idxs, key=lambda i: int(rows[i][4]))
        if newend > int(rows[term][4]):
            rows[term][4] = str(newend); n_ext += 1
    else:
        term = min(idxs, key=lambda i: int(rows[i][3]))
        if newend + 1 < int(rows[term][3]):
            rows[term][3] = str(newend + 1); n_ext += 1
for f in rows:
    if isinstance(f, str):
        continue
    if f[2] in ("gene", "transcript"):
        g = gid_of(f[8])
        if g in gene_new_end:
            strand, newend = gene_new_end[g]
            if strand == "+" and newend > int(f[4]):
                f[4] = str(newend)
            elif strand == "-" and newend + 1 < int(f[3]):
                f[3] = str(newend + 1)
with open(out, "w") as w:
    for f in rows:
        w.write(f if isinstance(f, str) else "\t".join(f) + "\n")
print(f"extended {n_ext} terminal exons -> {out}")
