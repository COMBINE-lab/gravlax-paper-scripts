#!/usr/bin/env bash
# Genome-wide recurrent-event discovery over the eight-donor SEZ collection versus a naive
# per-archive discovery-and-merge baseline. No coordinates are supplied to the collection search.
# Archived from the analysis host as scripts/156_collection_discovery_v023.sh (2026-09-11); renumbered here to follow the repository sequence.
set -euo pipefail
P=${GRAVLAX_PROJECT_ROOT:?set GRAVLAX_PROJECT_ROOT to the analysis project root (see README)}
AIE=${AIE:?set AIE to the released 0.2.3 binary}
PY=${GRAVLAX_PYTHON:-python3}
O=$P/runs/paper-followup-20260911/collection-discovery-0.2.3
ARCH=$P/runs/post-v1/archive-root-gate-a-r3/archives
ROUTED=$P/runs/post-v1/jshape-routing-gate-b-r3/collections/canonical/route-root.aicollection
COORD=$P/runs/post-v1/jshape-routing-gate-b-r3/collections/canonical/fallback-root.aicollection
AIC=$P/runs/genefull-20260910/gencode-v32.aic
G=$P/gravlax-paper-scripts/results/post-v1-sez-transcript-end-feasibility/groups
CPUSET=0-23; export RAYON_NUM_THREADS=24
cd $O
# inputs: design and a combined sample-scoped group map (astro_nsc, mature_neuron, ...)
true
true
THR=(--kind cassette --design design.tsv --groups groups.tsv --require-group astro_nsc --require-group mature_neuron --min-group-umi-classes 2 --min-donors 4 --min-samples 4 --min-umi-classes 8 --min-side-umi-classes 2 --min-support 2)
ANN=(--annotation $AIC --assembly GRCh38 --annotation-label GENCODE-v32)
run() { # label collection extra...
  local label=$1 coll=$2; shift 2
  echo "[$(date)] find-events $label"
  /usr/bin/time -v -o $O/find-$label.time.txt taskset -c $CPUSET $AIE collection find-events $coll "${THR[@]}" "$@" \
    --max-candidates 10000000 --max-candidates-considered 100000000 --max-routed-entries 1000000000 --max-exact-match-attempts 2500000000 --max-annotation-comparisons 1000000000 \
    --format json -o $O/find-$label.json > $O/find-$label.stdout.txt 2> $O/find-$label.stderr.txt || echo "find-events $label FAILED (see stderr)"
}
run routed-novel-v32 $ROUTED "${ANN[@]}" --novel-only
run coord-novel-v32  $COORD  "${ANN[@]}" --novel-only
run routed-all       $ROUTED "${ANN[@]}"
$PY - $O <<'PYEOF'
import json,sys,glob,re,collections
O=sys.argv[1]
def wall(f):
    t=open(f).read(); m=re.search(r"Elapsed \(wall clock\).*?: (.*)",t); h=m.group(1).strip().split(":")
    s=float(h[-1])+60*float(h[-2])+(3600*float(h[-3]) if len(h)==3 else 0)
    r=re.search(r"Maximum resident set size \(kbytes\): (\d+)",t); return s,int(r.group(1))/1048576
out=["# Collection discovery benchmark",""]
for lab in ["routed-novel-v32","coord-novel-v32","routed-all"]:
    try:
        s,g=wall(f"{O}/find-{lab}.time.txt")
        d=json.load(open(f"{O}/find-{lab}.json")); summ=d.get("data",{}).get("summary",{})
        keys={k:summ[k] for k in summ if any(t in k for t in ("candidates","entities","bytes","exact_match","archives","routed","retained","reported","classified","novel"))}
        out.append(f"## find-events {lab}: {s:.1f} s wall, {g:.2f} GiB peak"); out.append("```"); out.append(json.dumps(keys,indent=1)[:3000]); out.append("```")
        tabs=d.get("data",{}).get("tables",{})
        out.append(f"tables: {list(tabs.keys())}")
    except Exception as e: out.append(f"## find-events {lab}: unavailable ({e})")
tot=0.0; per=collections.Counter()
for line in open(f"{O}/naive/times.tsv"):
    s,c,t=line.split("\t"); tot+=float(t); per[s]+=float(t)
out.append(f"## naive per-archive discovery: {tot:.1f} s total wall over {sum(1 for _ in open(f'{O}/naive/times.tsv'))} runs; per archive: "+", ".join(f"{k} {v:.1f}s" for k,v in sorted(per.items())))
open(f"{O}/summary.md","w").write("\n".join(out)); print("\n".join(out))
PYEOF
echo "COLLECTION DISCOVERY RERUN DONE"
