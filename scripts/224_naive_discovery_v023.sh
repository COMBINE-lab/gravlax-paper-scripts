#!/usr/bin/env bash
# Archived from the analysis host as scripts/156b_naive_v023.sh (2026-09-11); renumbered here to follow the repository sequence.
# Genome-wide recurrent-event discovery over the eight-donor SEZ collection versus a naive
# per-archive discovery-and-merge baseline. No coordinates are supplied to the collection search.
set -euo pipefail
P=${GRAVLAX_PROJECT_ROOT:?set GRAVLAX_PROJECT_ROOT to the analysis project root (see README)}
AIE=$P/runs/gravlax-release-0.2.3-20260912/aie-0.2.3
PY=${GRAVLAX_PYTHON:-python3}
O=$P/runs/paper-followup-20260911/collection-discovery-0.2.3
ARCH=$P/runs/post-v1/archive-root-gate-a-r3/archives
ROUTED=$P/runs/post-v1/jshape-routing-gate-b-r3/collections/canonical/route-root.aicollection
COORD=$P/runs/post-v1/jshape-routing-gate-b-r3/collections/canonical/fallback-root.aicollection
AIC=$P/runs/genefull-20260910/gencode-v32.aic
G=$P/gravlax-paper-scripts/results/post-v1-sez-transcript-end-feasibility/groups
CPUSET=0-23; export RAYON_NUM_THREADS=24
THR=(--kind cassette --design design.tsv --groups groups.tsv --require-group astro_nsc --require-group mature_neuron --min-group-umi-classes 2 --min-donors 4 --min-samples 4 --min-umi-classes 8 --min-side-umi-classes 2 --min-support 2)
ANN=(--annotation $AIC --assembly GRCh38 --annotation-label GENCODE-v32)
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
