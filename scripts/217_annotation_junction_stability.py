#!/usr/bin/env python
# Archived from the analysis host as scripts/152_annotation_junction_stability.py (2026-09-11); renumbered here to follow the repository sequence.
"""How much do transcript sets and splice-junction sets change between GENCODE releases,
and how much observed read support falls on junctions that a given release lacks?

Outputs: annotation-stability/summary.json and summary.md
"""
import os
import json, sys, gzip, collections
from pathlib import Path
P=Path(os.environ["GRAVLAX_PROJECT_ROOT"])
OUT=P/"runs/paper-followup-20260911/annotation-stability"; OUT.mkdir(parents=True, exist_ok=True)
VERSIONS={"v32":P/"annotations/gencode.v32.annotation.gtf","v48":P/"annotations/gencode.v48.annotation.gtf","v49":P/"annotations/gencode.v49.annotation.gtf"}
SJ={"D0":P/"runs/ingest/full/SJ.out.tab","D1":P/"runs/ingest/d1/SJ.out.tab","D2p":P/"runs/ingest/d2p/SJ.out.tab"}

def strip(x): return x.split(".")[0] if not x.endswith("_PAR_Y") else x.split(".")[0]+"_PAR_Y"

def parse(gtf):
    tx=collections.defaultdict(list); tgene={}; tstrand={}; tchr={}; ttype={}
    with open(gtf) as f:
        for line in f:
            if line[0]=="#": continue
            c=line.split("\t",9)
            if c[2]!="exon": continue
            attrs=c[8]
            i=attrs.find('transcript_id "'); tid=attrs[i+15:attrs.find('"',i+15)]
            j=attrs.find('gene_id "'); gid=attrs[j+9:attrs.find('"',j+9)]
            k=attrs.find('transcript_type "'); tt=attrs[k+17:attrs.find('"',k+17)] if k>=0 else "NA"
            tid=strip(tid); gid=strip(gid)
            tx[tid].append((int(c[3]),int(c[4]))); tgene[tid]=gid; tstrand[tid]=c[6]; tchr[tid]=c[0]; ttype[tid]=tt
    junctions=set(); jtx=collections.Counter(); structs={}; single=0
    for tid,ex in tx.items():
        ex=sorted(set(ex)); chr_,st=tchr[tid],tstrand[tid]
        introns=tuple((ex[i][1]+1, ex[i+1][0]-1) for i in range(len(ex)-1) if ex[i+1][0]-1>=ex[i][1]+1)
        if not introns:
            single+=1; structs[tid]=("single",chr_,st,ex[0][0],ex[-1][1])
        else:
            structs[tid]=("chain",chr_,st,introns)
        for a,b in introns:
            junctions.add((chr_,a,b,st)); jtx[(chr_,a,b,st)]+=1
    return dict(transcripts=set(tx), genes=set(tgene.values()), structs=structs, junctions=junctions,
                jtx=jtx, single=single, ttype=ttype, tgene=tgene)

A={v:parse(p) for v,p in VERSIONS.items()}
res={"versions":{},"pairs":{},"observed":{}}
for v,d in A.items():
    res["versions"][v]=dict(genes=len(d["genes"]),transcripts=len(d["transcripts"]),single_exon_transcripts=d["single"],
                            distinct_intron_chains=len({s for s in d["structs"].values() if s[0]=="chain"}),
                            junctions=len(d["junctions"]),
                            junctions_in_one_transcript=sum(1 for j,n in d["jtx"].items() if n==1),
                            protein_coding_transcripts=sum(1 for t,tt in d["ttype"].items() if tt=="protein_coding"))
def jac(a,b): return len(a&b)/len(a|b)
for x,y in [("v32","v48"),("v48","v49"),("v32","v49")]:
    a,b=A[x],A[y]
    chains_a={s for s in a["structs"].values() if s[0]=="chain"}; chains_b={s for s in b["structs"].values() if s[0]=="chain"}
    structs_a=set(a["structs"].values()); structs_b=set(b["structs"].values())
    r={}
    for name,sa,sb in [("genes_by_id",a["genes"],b["genes"]),("transcripts_by_id",a["transcripts"],b["transcripts"]),
                       ("transcripts_by_structure",structs_a,structs_b),("intron_chains",chains_a,chains_b),("junctions",a["junctions"],b["junctions"])]:
        r[name]=dict(old=len(sa),new=len(sb),shared=len(sa&sb),removed=len(sa-sb),added=len(sb-sa),
                     jaccard=round(jac(sa,sb),4), frac_of_new_present_in_old=round(len(sa&sb)/len(sb),4),
                     frac_of_old_retained=round(len(sa&sb)/len(sa),4))
    # transcripts whose ID persists but whose structure changed
    both=a["transcripts"]&b["transcripts"]
    changed=sum(1 for t in both if a["structs"][t]!=b["structs"][t])
    r["transcripts_by_id"]["shared_with_changed_structure"]=changed
    res["pairs"][f"{x}->{y}"]=r
# observed junction support (annotation-free STAR two-pass, unique reads)
for ds,path in SJ.items():
    obs={}
    with open(path) as f:
        for line in f:
            c=line.split("\t"); chr_,a,b,st,uniq=c[0],int(c[1]),int(c[2]),c[3],int(c[6])
            if uniq<=0: continue
            obs[(chr_,a,b,{"1":"+","2":"-"}.get(st,"."))]=uniq
    def inset(j,S):
        chr_,a,b,st=j
        if st==".": return (chr_,a,b,"+") in S or (chr_,a,b,"-") in S
        return j in S
    total_j=len(obs); total_r=sum(obs.values())
    o={"observed_junctions_unique_reads_ge1":total_j,"unique_read_support":total_r}
    for thr in (1,3,10):
        sel={j:n for j,n in obs.items() if n>=thr}
        tj=len(sel); tr=sum(sel.values()); row={"junctions":tj,"reads":tr}
        for v,d in A.items():
            inv=[j for j in sel if inset(j,d["junctions"])]
            row[f"in_{v}_junction_frac"]=round(len(inv)/tj,4); row[f"in_{v}_read_frac"]=round(sum(sel[j] for j in inv)/tr,4)
        # in v49 but not v32 (what a v32 junction annotation would miss among annotated-in-v49)
        gained=[j for j in sel if inset(j,A["v49"]["junctions"]) and not inset(j,A["v32"]["junctions"])]
        lost=[j for j in sel if inset(j,A["v32"]["junctions"]) and not inset(j,A["v49"]["junctions"])]
        none=[j for j in sel if not any(inset(j,A[v]["junctions"]) for v in A)]
        row["in_v49_not_v32_junction_frac"]=round(len(gained)/tj,4); row["in_v49_not_v32_read_frac"]=round(sum(sel[j] for j in gained)/tr,4)
        row["in_v32_not_v49_junction_frac"]=round(len(lost)/tj,4); row["in_v32_not_v49_read_frac"]=round(sum(sel[j] for j in lost)/tr,4)
        row["in_no_version_junction_frac"]=round(len(none)/tj,4); row["in_no_version_read_frac"]=round(sum(sel[j] for j in none)/tr,4)
        o[f"min_unique_reads_{thr}"]=row
    res["observed"][ds]=o
json.dump(res,open(OUT/"summary.json","w"),indent=1)
# markdown summary
L=["# Annotation stability: transcripts versus junctions","","## Per version",""]
L.append("| version | genes | transcripts | intron chains | junctions | junctions in one transcript |"); L.append("|---|---:|---:|---:|---:|---:|")
for v,r in res["versions"].items(): L.append(f"| {v} | {r['genes']:,} | {r['transcripts']:,} | {r['distinct_intron_chains']:,} | {r['junctions']:,} | {r['junctions_in_one_transcript']:,} |")
L+=["","## Pairwise change","","| pair | set | old | new | shared | removed | added | Jaccard | new present in old |","|---|---|---:|---:|---:|---:|---:|---:|---:|"]
for pr,r in res["pairs"].items():
    for name,s in r.items(): L.append(f"| {pr} | {name} | {s['old']:,} | {s['new']:,} | {s['shared']:,} | {s['removed']:,} | {s['added']:,} | {s['jaccard']:.3f} | {s['frac_of_new_present_in_old']:.3f} |")
L+=["","## Observed junctions (annotation-free two-pass STAR, uniquely mapped reads)",""]
for ds,o in res["observed"].items():
    L.append(f"### {ds}: {o['observed_junctions_unique_reads_ge1']:,} junctions, {o['unique_read_support']:,} unique junction reads"); L.append("")
    L.append("| min reads | junctions | reads | in v32 (j / reads) | in v48 | in v49 | v49-only vs v32 | v32-only vs v49 | in none |"); L.append("|---:|---:|---:|---|---|---|---|---|---|")
    for thr in (1,3,10):
        r=o[f"min_unique_reads_{thr}"]
        L.append(f"| {thr} | {r['junctions']:,} | {r['reads']:,} | {r['in_v32_junction_frac']:.3f} / {r['in_v32_read_frac']:.3f} | {r['in_v48_junction_frac']:.3f} / {r['in_v48_read_frac']:.3f} | {r['in_v49_junction_frac']:.3f} / {r['in_v49_read_frac']:.3f} | {r['in_v49_not_v32_junction_frac']:.4f} / {r['in_v49_not_v32_read_frac']:.4f} | {r['in_v32_not_v49_junction_frac']:.4f} / {r['in_v32_not_v49_read_frac']:.4f} | {r['in_no_version_junction_frac']:.3f} / {r['in_no_version_read_frac']:.3f} |")
    L.append("")
open(OUT/"summary.md","w").write("\n".join(L)); print("\n".join(L))
