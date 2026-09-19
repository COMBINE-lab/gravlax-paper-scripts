#!/usr/bin/env bash
# Genome-wide recurrent-event discovery over the eight-donor SEZ collection versus a naive
# per-archive discovery-and-merge baseline. No coordinates are supplied to the collection search.
# Archived from the analysis host as scripts/154_collection_discovery.sh (2026-09-11); renumbered here to follow the repository sequence.
set -euo pipefail
P=${GRAVLAX_PROJECT_ROOT:?set GRAVLAX_PROJECT_ROOT to the analysis project root (see README)}
AIE=$P/runs/gravlax-release-0.2.2-20260910/aie-0.2.2
PY=${GRAVLAX_PYTHON:-python3}
O=$P/runs/paper-followup-20260911/collection-discovery
ARCH=$P/runs/post-v1/archive-root-gate-a-r3/archives
ROUTED=$P/runs/post-v1/jshape-routing-gate-b-r3/collections/canonical/route-root.aicollection
COORD=$P/runs/post-v1/jshape-routing-gate-b-r3/collections/canonical/fallback-root.aicollection
AIC=$P/runs/genefull-20260910/gencode-v32.aic
G=$P/gravlax-paper-scripts/results/post-v1-sez-transcript-end-feasibility/groups
CPUSET=0-23; export RAYON_NUM_THREADS=24
cd $O
# inputs: design and a combined sample-scoped group map (astro_nsc, mature_neuron, ...)
printf 'sample\tdonor\n' > design.tsv; printf 'sample\tbarcode\tgroup\n' > groups.tsv
for s in A B C D E F G H; do printf '%s\t%s\n' $s $s >> design.tsv; awk -v s=$s -F'\t' 'BEGIN{OFS="\t"}{print s,$1,$2}' $G/donor-$s.tsv >> groups.tsv; done
THR=(--kind cassette --design design.tsv --groups groups.tsv --require-group astro_nsc --require-group mature_neuron --min-group-umi-classes 2 --min-donors 4 --min-samples 4 --min-umi-classes 8 --min-side-umi-classes 2 --min-support 2)
ANN=(--annotation $AIC --assembly GRCh38 --annotation-label GENCODE-v32)
run() { # label collection extra...
  local label=$1 coll=$2; shift 2
  echo "[$(date)] find-events $label"
  /usr/bin/time -v -o $O/find-$label.time.txt taskset -c $CPUSET $AIE collection find-events $coll "${THR[@]}" "$@" \
    --max-candidates 1000000 --max-candidates-considered 10000000 --max-routed-entries 100000000 --max-exact-match-attempts 250000000 --max-annotation-comparisons 100000000 \
    --format json -o $O/find-$label.json > $O/find-$label.stdout.txt 2> $O/find-$label.stderr.txt || echo "find-events $label FAILED (see stderr)"
}
run routed-novel-v32 $ROUTED "${ANN[@]}" --novel-only
run coord-novel-v32  $COORD  "${ANN[@]}" --novel-only
run routed-all       $ROUTED "${ANN[@]}"
# naive baseline: per-archive, per-chromosome cassette discovery with the donor's own group map, then merge
echo "[$(date)] naive per-archive discovery"
mkdir -p $O/naive; : > $O/naive/times.tsv
CHR="chr1 248956422 chr2 242193529 chr3 198295559 chr4 190214555 chr5 181538259 chr6 170805979 chr7 159345973 chr8 145138636 chr9 138394717 chr10 133797422 chr11 135086622 chr12 133275309 chr13 114364328 chr14 107043718 chr15 101991189 chr16 90338345 chr17 83257441 chr18 80373285 chr19 58617616 chr20 64444167 chr21 46709983 chr22 50818468 chrX 156040895 chrY 57227415 chrM 16569"
for s in A B C D E F G H; do
  set -- $CHR
  while [ $# -gt 0 ]; do c=$1; l=$2; shift 2
    st=$(date +%s.%N)
    taskset -c $CPUSET $AIE query $ARCH/$s.v2.aie events $c:1-$l --event-type cassette --min-support 2 --min-informative 1 \
      --groups $G/donor-$s.tsv --agg group --max-events 1000000 --format json -o $O/naive/$s-$c.json > /dev/null 2> $O/naive/$s-$c.err || echo "naive $s $c FAILED"
    printf '%s\t%s\t%.3f\n' $s $c "$(echo "$(date +%s.%N) - $st" | bc)" >> $O/naive/times.tsv
  done
done
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
echo "COLLECTION DISCOVERY DONE"
