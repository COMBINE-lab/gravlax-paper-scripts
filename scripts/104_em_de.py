#!/usr/bin/env python
"""Does the EM layer change differential-expression conclusions, and is it right?
T-vs-monocyte pseudobulk log2FC per gene: unique-only vs unique+EM, with alevin-fry
(sequence-aware) as the independent arbiter."""
import numpy as np, scipy.io
P=os.environ["GRAVLAX_PROJECT_ROOT"]
AF="/fs/cbcb-lab/rob/pbmc_1k_v3_simpleaf/af_quant/alevin"

groups={}
for l in open(P+"/runs/replay/apa-groups-d0.tsv"):
    b,g=l.split()
    groups[b]=g
def grouped(mtx_dir, mtx="matrix.mtx", bc="barcodes.tsv", feat="features.tsv", gene_col=0):
    m=scipy.io.mmread(f"{mtx_dir}/{mtx}").tocsc()
    bcs=[l.strip() for l in open(f"{mtx_dir}/{bc}")]
    genes=[l.split("\t")[gene_col].split(".")[0].split("-")[0] for l in open(f"{mtx_dir}/{feat}")]
    tcols=[i for i,b in enumerate(bcs) if groups.get(b)=="T"]
    mcols=[i for i,b in enumerate(bcs) if groups.get(b)=="M"]
    t=np.asarray(m[:,tcols].sum(axis=1)).ravel(); mo=np.asarray(m[:,mcols].sum(axis=1)).ravel()
    dt={};dm={}
    for g,a,b2 in zip(genes,t,mo):
        dt[g]=dt.get(g,0)+a; dm[g]=dm.get(g,0)+b2
    return dt,dm
uT,uM=grouped(P+"/runs/replay/full-aie")
eT,eM=grouped(P+"/runs/replay/d0-em-layer",mtx="em.mtx")
class AFdir: pass
import os
afT,afM=grouped(AF, mtx="quants_mat.mtx", bc="quants_mat_rows.txt", feat="quants_mat_cols.txt")
# library-size normalize each side
def norm(d):
    s=sum(d.values()); return {g:v/s*1e6 for g,v in d.items()}
def fc(dT,dM,g):
    return np.log2((dT.get(g,0)+1)/(dM.get(g,0)+1))
nuT,nuM=norm(uT),norm(uM)
ceT={g:uT.get(g,0)+eT.get(g,0) for g in set(uT)|set(eT)}
ceM={g:uM.get(g,0)+eM.get(g,0) for g in set(uM)|set(eM)}
nceT,nceM=norm(ceT),norm(ceM)
nafT,nafM=norm(afT),norm(afM)
# genes: top EM recipients with decent expression both sides
top=sorted(eT, key=lambda g: -(eT.get(g,0)+eM.get(g,0)))[:400]
rows=[]
for g in top:
    if min(uT.get(g,0),uM.get(g,0))<20 or g not in nafT: continue
    f_u=fc(nuT,nuM,g); f_e=fc(nceT,nceM,g); f_a=fc(nafT,nafM,g)
    rows.append((g,f_u,f_e,f_a))
f_u=np.array([r[1] for r in rows]); f_e=np.array([r[2] for r in rows]); f_a=np.array([r[3] for r in rows])
moved=np.abs(f_e-f_u)>0.5
closer=np.abs(f_e-f_a)<np.abs(f_u-f_a)
print(f"genes analyzed (EM-heavy, expressed both groups, in alevin-fry): {len(rows)}")
print(f"genes whose T-vs-mono log2FC moves >0.5 with the EM layer: {moved.sum()}")
print(f"median |FC - FC_alevinfry|: unique-only {np.median(np.abs(f_u-f_a)):.3f} -> with EM {np.median(np.abs(f_e-f_a)):.3f}")
print(f"genes moved CLOSER to alevin-fry's FC: {closer.sum()}/{len(rows)} ({100*closer.mean():.0f}%)")
amoved=[(abs(e-u),g,u,e,a) for (g,u,e,a) in rows]
amoved.sort(reverse=True)
print("largest FC corrections (gene: unique-FC -> EM-FC | alevin-fry-FC):")
for d,g,u,e,a in amoved[:6]:
    print(f"  {g}: {u:+.2f} -> {e:+.2f} | af {a:+.2f}")
