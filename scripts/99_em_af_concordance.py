#!/usr/bin/env python
"""EM external concordance: does adding the pooled-EM recovered layer move gene totals TOWARD an
independent sequence-aware quantifier (alevin-fry, selective alignment + its own EM)?"""
import os
import numpy as np, scipy.io
P=os.environ["GRAVLAX_PROJECT_ROOT"]
AF="/fs/cbcb-lab/rob/pbmc_1k_v3_simpleaf/af_quant/alevin"
m=scipy.io.mmread(AF+"/quants_mat.mtx").tocsr()
cols=[l.strip() for l in open(AF+"/quants_mat_cols.txt")]
af=np.asarray(m.sum(axis=0)).ravel()
af_gene={}
for c,v in zip(cols,af):
    g=c.split("-")[0].split(".")[0]
    af_gene[g]=af_gene.get(g,0.0)+v
def sums(d,mtx="matrix.mtx"):
    mm=scipy.io.mmread(f"{d}/{mtx}").tocsr()
    genes=[l.split("\t")[0].split(".")[0] for l in open(d+"/features.tsv")]
    return dict(zip(genes,np.asarray(mm.sum(axis=1)).ravel()))
uniq=sums(P+"/runs/replay/full-aie")            # exact Gene replay (v49)
em=sums(P+"/runs/replay/d0-em-layer","em.mtx")  # recovered layer
top=sorted(em,key=lambda g:-em.get(g,0))[:300]  # multimapper-heavy genes
rows=[(g,uniq.get(g,0),uniq.get(g,0)+em.get(g,0),af_gene.get(g,0)) for g in top if af_gene.get(g,0)>0]
u=np.array([r[1] for r in rows]); ue=np.array([r[2] for r in rows]); a=np.array([r[3] for r in rows])
def l1rel(x): return np.median(np.abs(np.log2((x+1)/(a+1))))
cu=np.corrcoef(np.log1p(u),np.log1p(a))[0,1]; cue=np.corrcoef(np.log1p(ue),np.log1p(a))[0,1]
print(f"top-{len(rows)} EM-recipient genes vs alevin-fry:")
print(f"  unique-only:   log-corr {cu:.4f}, median |log2 ratio| {l1rel(u):.3f}")
print(f"  unique + EM:   log-corr {cue:.4f}, median |log2 ratio| {l1rel(ue):.3f}")
closer=int(np.sum(np.abs(np.log2((ue+1)/(a+1)))<np.abs(np.log2((u+1)/(a+1)))))
print(f"  genes moved CLOSER to alevin-fry by the EM layer: {closer}/{len(rows)} ({100*closer/len(rows):.0f}%)")
