#!/usr/bin/env python
"""Prospective discovery: given the real 2020 annotation (v32), does discover+requant find and
correctly quantify the genes GENCODE only added by v49? External truth = annotation history."""
import os
import numpy as np, scipy.io
P=os.environ["GRAVLAX_PROJECT_ROOT"]
def gene_ids(gtf):
    out={}
    for line in open(gtf):
        if line.startswith("#"): continue
        f=line.split("\t")
        if len(f)<9 or f[2]!="gene": continue
        i=f[8].find('gene_id "'); g=f[8][i+9:f[8].index('"',i+9)]
        out[g.split(".")[0]]=(f[0],int(f[3])-1,int(f[4]))
    return out
v32_spans=gene_ids(P+"/annotations/gencode.v32.annotation.gtf")
v32=set(v32_spans)
import bisect
by32={}
for g,(ch,s2,e2) in v32_spans.items():
    by32.setdefault(ch,[]).append((s2,e2))
for ch in by32: by32[ch].sort()
def ov32(ch,s2,e2):
    lst=by32.get(ch,[])
    i=bisect.bisect_left(lst,(e2,0))
    return any(lst[j][0]<e2 and lst[j][1]>s2 for j in range(max(0,i-2000),min(len(lst),i+1)))
v49=gene_ids(P+"/annotations/gencode.v49.annotation.gtf")
new=set(v49)-v32
cells={l.strip() for l in open(P+"/runs/oracle/full/v49/Solo.out/Gene/filtered/barcodes.tsv")}
def sums(d):
    m=scipy.io.mmread(d+"/matrix.mtx").tocsr()
    genes=[l.split("\t")[0].split(".")[0] for l in open(d+"/features.tsv")]
    bcs=[l.strip() for l in open(d+"/barcodes.tsv")]
    cols=[i for i,b in enumerate(bcs) if b in cells]
    return dict(zip(genes,np.asarray(m[:,cols].sum(axis=1)).ravel()))
o_gene=sums(P+"/runs/oracle/full/v49/Solo.out/Gene/raw")
o_full=sums(P+"/runs/oracle/full/v49/Solo.out/GeneFull/raw")
novel=sums(P+"/runs/replay/d0-novel-v32")
nspans={}
for line in open(P+"/runs/replay/d0-novel-v32.gtf"):
    f=line.split("\t")
    if f[2]!="gene": continue
    i=f[8].find('gene_id "'); g=f[8][i+9:f[8].index('"',i+9)]
    nspans.setdefault(f[0],[]).append((int(f[3])-1,int(f[4]),g))
stats={True:[0,0,[]],False:[0,0,[]]}  # overlaps_v32 -> [tot,hit,errs]
for g in new:
    if o_gene.get(g,0)<50: continue
    ch,s,e=v49[g]
    key=ov32(ch,s,e)
    stats[key][0]+=1
    got=sum(novel.get(cg,0) for (cs,ce,cg) in nspans.get(ch,[]) if cs<e and ce>s)
    if got>0:
        stats[key][1]+=1
        oracle=o_full.get(g,0)
        stats[key][2].append(abs(got-oracle)/max(oracle,1))
for lab,key in [("clean (no v32 overlap)",False),("overlapping v32 genes (claimed-by-design)",True)]:
    tot,hit,errs=stats[key]; errs=np.array(errs)
    med=100*np.median(errs) if len(errs) else float('nan')
    print(f"{lab}: {hit}/{tot} recovered ({100*hit/max(tot,1):.1f}%), median err {med:.1f}% (n={len(errs)})")
