set -x
P=${GRAVLAX_PROJECT_ROOT:?set GRAVLAX_PROJECT_ROOT to the analysis project root (see README)}
PY=${GRAVLAX_PYTHON:-python3}
AIE=$P/src/target/release/aie
cd $P
# 1. APA replication
$PY scripts/100_apa_replicate.py run d0 runs/oracle/full/v49/Solo.out/Gene/filtered runs/archive/d0.aie results/apa-shifts-d0.tsv
$PY scripts/100_apa_replicate.py run d1 runs/oracle/d1/v49/Solo.out/Gene/filtered runs/archive/d1.aie results/apa-shifts-d1.tsv
$PY scripts/100_apa_replicate.py analyze
echo APA_REPLICATION_DONE
# 2. Intergenic-locus replication (discover vs v49 on both PBMC datasets)
$AIE query runs/archive/d0.aie discover --gtf annotations/gencode.v49.annotation.gtf --tsv > results/discover-v49-d0.tsv
$AIE query runs/archive/d1.aie discover --gtf annotations/gencode.v49.annotation.gtf --tsv > results/discover-v49-d1.tsv
$PY - <<'PYEOF'
import bisect
P=${GRAVLAX_PROJECT_ROOT:?set GRAVLAX_PROJECT_ROOT to the analysis project root (see README)}
v49={}
for line in open(P+"/annotations/gencode.v49.annotation.gtf"):
    if line.startswith("#"): continue
    f=line.split("\t")
    if f[2]!="gene": continue
    v49.setdefault(f[0],[]).append((int(f[3])-1,int(f[4]),f[6]))
for ch in v49: v49[ch].sort()
def near_or_ov(ch,s,e,d=5000):
    lst=v49.get(ch,[])
    i=bisect.bisect_left(lst,(e+d,0,""))
    for j in range(max(0,i-2000),min(len(lst),i+1)):
        gs,ge,st=lst[j]
        if gs<e and ge>s: return True
        if st=="+" and 0<=s-ge<d: return True
        if st=="-" and 0<=gs-e<d: return True
    return False
def intergenic(path,minu=10):
    out=[]
    for l in open(path):
        ch,s,e,st,u,c=l.rstrip("\n").split("\t")
        s,e,u=int(s),int(e),int(u)
        if u>=minu and not near_or_ov(ch,s,e): out.append((ch,s,e,u))
    return out
d0=intergenic(P+"/results/discover-v49-d0.tsv")
d1=intergenic(P+"/results/discover-v49-d1.tsv")
by1={}
for ch,s,e,u in d1: by1.setdefault(ch,[]).append((s,e))
for ch in by1: by1[ch].sort()
rep=0; repmass=0; totmass=0
for ch,s,e,u in d0:
    totmass+=u
    lst=by1.get(ch,[])
    i=bisect.bisect_left(lst,(e,0))
    if any(lst[j][0]<e and lst[j][1]>s for j in range(max(0,i-500),min(len(lst),i+1))):
        rep+=1; repmass+=u
print(f"intergenic candidates: D0 {len(d0)}, D1 {len(d1)}")
print(f"D0 intergenic loci replicated in D1: {rep}/{len(d0)} ({100*rep/max(len(d0),1):.1f}% of loci, {100*repmass/max(totmass,1):.1f}% of mass)")
PYEOF
echo INTERGENIC_REPLICATION_DONE
# 3. EM recovered-fraction replication
$AIE em runs/archive/d1.aie --gtf annotations/gencode.v49.annotation.gtf --mask 0 \
  --emit runs/replay/d1-em-layer --barcodes runs/oracle/d1/v49/Solo.out/Gene/raw/barcodes.tsv 2>&1 | tail -1
$PY - <<'PYEOF'
import numpy as np, scipy.io
P=${GRAVLAX_PROJECT_ROOT:?set GRAVLAX_PROJECT_ROOT to the analysis project root (see README)}
def sums(d,mtx="matrix.mtx"):
    m=scipy.io.mmread(f"{d}/{mtx}").tocsr()
    genes=[l.split("\t")[0].split(".")[0] for l in open(d+"/features.tsv")]
    return dict(zip(genes,np.asarray(m.sum(axis=1)).ravel()))
u0=sums(P+"/runs/replay/full-aie"); e0=sums(P+"/runs/replay/d0-em-layer","em.mtx")
u1=sums(P+"/runs/replay/d1-v15"); e1=sums(P+"/runs/replay/d1-em-layer","em.mtx")
rows=[]
for g in set(e0)&set(e1):
    t0=u0.get(g,0)+e0.get(g,0); t1=u1.get(g,0)+e1.get(g,0)
    if t0>=100 and t1>=100 and (e0.get(g,0)>0 or e1.get(g,0)>0):
        rows.append((e0.get(g,0)/t0, e1.get(g,0)/t1))
f0=np.array([r[0] for r in rows]); f1=np.array([r[1] for r in rows])
print(f"EM recovered-fraction replication (D0 vs D1, {len(rows)} genes with total>=100 in both):")
print(f"  Pearson r = {np.corrcoef(f0,f1)[0,1]:.4f}; median |delta fraction| = {np.median(np.abs(f0-f1)):.4f}")
PYEOF
echo EM_REPLICATION_DONE
echo REPLICATION_PACK_DONE
