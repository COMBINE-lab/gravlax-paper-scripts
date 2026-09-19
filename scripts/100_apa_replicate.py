#!/usr/bin/env python
"""APA replication + annotated-3'-end concordance.
Part 1: run the T-vs-mono contrast on a dataset, emit per-gene (real, null) shifts.
Part 2 (d0+d1 done): overlap of significant genes; for D0 significants, do the two populations'
dominant sites coincide with distinct annotated transcript 3' ends (<=200 bp)?"""
import os
import random, subprocess, sys
from pathlib import Path
import numpy as np, scipy.io

P = Path(os.environ["GRAVLAX_PROJECT_ROOT"])
AIE = P/"src/target/release/aie"
MIN_UMIS = 30
random.seed(11)

def gene_spans():
    spans={}
    for line in open(P/"annotations/gencode.v49.annotation.gtf"):
        if line.startswith("#"): continue
        f=line.split("\t")
        if f[2]!="gene": continue
        i=f[8].find('gene_name "')
        if i<0: continue
        name=f[8][i+11:f[8].index('"',i+11)]
        spans.setdefault(name,(f[0],int(f[3])-1,int(f[4]),f[6]))
    return spans

def tx_ends():
    ends={}
    for line in open(P/"annotations/gencode.v49.annotation.gtf"):
        if line.startswith("#"): continue
        f=line.split("\t")
        if f[2]!="transcript": continue
        i=f[8].find('gene_name "')
        if i<0: continue
        name=f[8][i+11:f[8].index('"',i+11)]
        e = int(f[4]) if f[6]=="+" else int(f[3])-1
        ends.setdefault(name,set()).add(e)
    return ends

def run_dataset(tag, filt_dir, archive, out_tsv):
    d=Path(filt_dir)
    m=scipy.io.mmread(str(d/"matrix.mtx")).tocsr()
    genes=[l.split("\t")[1] for l in open(d/"features.tsv")]
    bcs=[l.strip() for l in open(d/"barcodes.tsv")]
    gi={g:i for i,g in enumerate(genes)}
    cd3e=np.asarray(m[gi["CD3E"],:].todense()).ravel()
    lyz=np.asarray(m[gi["LYZ"],:].todense()).ravel()
    t=[b for b,c,l in zip(bcs,cd3e,lyz) if c>0 and l==0]
    mo=[b for b,c,l in zip(bcs,cd3e,lyz) if l>3 and c==0]
    gr=P/f"runs/replay/apa-groups-{tag}.tsv"; gn=P/f"runs/replay/apa-null-{tag}.tsv"
    with open(gr,"w") as f:
        for b in t: f.write(f"{b}\tT\n")
        for b in mo: f.write(f"{b}\tM\n")
    half=set(random.sample(t,len(t)//2))
    with open(gn,"w") as f:
        for b in t: f.write(f"{b}\t{'A' if b in half else 'B'}\n")
    spans=gene_spans()
    expr=np.asarray(m.sum(axis=1)).ravel()
    top=[genes[i] for i in np.argsort(-expr)[:250] if genes[i] in spans]
    def shifts(gfile,ga,gb):
        out={}
        for g in top:
            ch,s,e,st=spans[g]
            r=subprocess.run([str(AIE),"query",archive,"apa",f"{ch}:{s}-{e}","--strand",st,"--tsv","--groups",str(gfile)],capture_output=True,text=True)
            hdr=None; pa=[]; pb=[]
            for line in r.stdout.splitlines():
                f3=line.rstrip("\n").split("\t")
                if line.startswith("#"): hdr=f3; continue
                mid=(int(f3[1])+int(f3[2]))/2
                pa.append((mid,int(f3[hdr.index(ga)]))); pb.append((mid,int(f3[hdr.index(gb)])))
            def wm(pp):
                tot=sum(c for _,c in pp)
                return (sum(p*c for p,c in pp)/tot,tot) if tot>=MIN_UMIS else (None,tot)
            (ma,_),(mb,_)=wm(pa),wm(pb)
            if ma is not None and mb is not None: out[g]=abs(ma-mb)
        return out
    real=shifts(gr,"T","M"); null=shifts(gn,"A","B")
    with open(out_tsv,"w") as f:
        for g in sorted(set(real)&set(null)):
            f.write(f"{g}\t{real[g]:.1f}\t{null[g]:.1f}\n")
    print(f"{tag}: {len(t)} T / {len(mo)} mono; {len(set(real)&set(null))} genes measured")

def analyze():
    def load(tag):
        d={}
        for l in open(P/f"results/apa-shifts-{tag}.tsv"):
            g,r,n=l.split("\t"); d[g]=(float(r),float(n))
        return d
    d0=load("d0"); d1=load("d1")
    def sig(d):
        p95=np.percentile([v[1] for v in d.values()],95)
        return {g for g,(r,_) in d.items() if r>p95}, p95
    s0,p0=sig(d0); s1,p1=sig(d1)
    common=set(d0)&set(d1)
    s0c=s0&common; s1c=s1&common
    inter=s0c&s1c
    print(f"D0: {len(s0)}/{len(d0)} significant (P95={p0:.0f}bp); D1: {len(s1)}/{len(d1)} (P95={p1:.0f}bp)")
    exp=len(s0c)*len(s1c)/max(len(common),1)
    print(f"replication: {len(inter)} genes significant in BOTH (of {len(s0c)} D0-sig testable in D1; expected by chance {exp:.1f})")
    print("replicated genes:", ", ".join(sorted(inter)[:20]))
    # 3'-end concordance for D0 significants
    spans=gene_spans(); ends=tx_ends()
    conc=both_near=tot=0
    for g in sorted(s0):
        if g not in spans or g not in ends: continue
        ch,s,e,st=spans[g]
        r=subprocess.run([str(AIE),"query",str(P/"runs/archive/d0.aie"),"apa",f"{ch}:{s}-{e}","--strand",st,"--tsv","--groups",str(P/"runs/replay/apa-groups-d0.tsv")],capture_output=True,text=True)
        hdr=None; best={"T":(0,None),"M":(0,None)}
        for line in r.stdout.splitlines():
            f3=line.rstrip("\n").split("\t")
            if line.startswith("#"): hdr=f3; continue
            mid=(int(f3[1])+int(f3[2]))//2
            for grp in ("T","M"):
                c=int(f3[hdr.index(grp)])
                if c>best[grp][0]: best[grp]=(c,mid)
        if best["T"][1] is None or best["M"][1] is None: continue
        tot+=1
        def near(p): return min((abs(p-x),x) for x in ends[g])
        dT,eT=near(best["T"][1]); dM,eM=near(best["M"][1])
        if dT<=200 and dM<=200:
            both_near+=1
            if eT!=eM: conc+=1
    print(f"3'-end concordance (D0 significants): {both_near}/{tot} genes have BOTH dominant sites within 200bp of annotated ends; {conc}/{both_near} use DISTINCT annotated ends")

if sys.argv[1]=="run":
    run_dataset(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
else:
    analyze()
