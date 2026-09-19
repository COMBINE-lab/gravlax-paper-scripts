#!/usr/bin/env python
"""Merge per-archive cassette discovery into cohort-recurrent events with the same rule as
collection find-events: >=8 pooled UMI classes, >=2 include-only and >=2 exclude-only pooled,
>=2 classes in each of astro_nsc and mature_neuron, in >=4 donors (=samples here)."""
# Archived from the analysis host as scripts/155_naive_merge.py (2026-09-11); renumbered here to follow the repository sequence.
import json,glob,os,sys,time,collections
O=sys.argv[1]; t0=time.time()
per=collections.defaultdict(dict)   # event -> sample -> dict(inc,exc,both,groups{g:(inc,exc)})
for f in sorted(glob.glob(f"{O}/naive/*-chr*.json")):
    s=os.path.basename(f).split("-")[0]
    d=json.load(open(f)); tabs={t["name"]:t for t in d["data"]["tables"]}
    ev={r[0]:r for r in tabs["events"]["rows"]}
    for r in tabs["counts"]["rows"]:
        eid,agg,ent,inc,exc,both=r[:6]
        if agg!="group": continue
        e=per[eid].setdefault(s,{"inc":0,"exc":0,"both":0,"groups":{}})
        e["inc"]+=inc; e["exc"]+=exc; e["both"]+=both; e["groups"][ent]=(inc,exc)
def ok(e):  # per-sample eligibility for the per-group requirement; pooled thresholds across samples
    g=e["groups"]; return sum(g.get("astro_nsc",(0,0)))>=2 and sum(g.get("mature_neuron",(0,0)))>=2
rec=[]
for eid,samples in per.items():
    supp={s:e for s,e in samples.items() if (e["inc"]+e["exc"])>0}
    if len(supp)<4: continue
    inc=sum(e["inc"] for e in supp.values()); exc=sum(e["exc"] for e in supp.values())
    if inc+exc<8 or inc<2 or exc<2: continue
    ga=sum(sum(e["groups"].get("astro_nsc",(0,0))) for e in supp.values()); gm=sum(sum(e["groups"].get("mature_neuron",(0,0))) for e in supp.values())
    if ga<2 or gm<2: continue
    rec.append((len(supp),inc+exc,eid))
rec.sort(reverse=True)
elapsed=time.time()-t0
fn="cassette:chr9:129908999-129915965:129919242-129923843:129908999-129923843"
out={"candidate_events_union":len(per),"recurrent_events":len(rec),"merge_seconds":round(elapsed,2),
     "fnbp1_present":any(e[2]==fn for e in rec),"top20":[{"donors":a,"informative":b,"event":c} for a,b,c in rec[:20]]}
json.dump(out,open(f"{O}/naive-merge.json","w"),indent=1); print(json.dumps(out,indent=1)[:2500])
