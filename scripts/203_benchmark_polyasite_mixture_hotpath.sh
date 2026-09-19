#!/usr/bin/env bash
# Re-run the primary arm with the sparse-likelihood hot path and require identical tables.
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 PROJECT_ROOT FROZEN_PRIMARY_DIR OUT_DIR" >&2
  exit 2
fi

GV_ROOT=$(realpath "$1")
GV_FROZEN=$(realpath "$2")
GV_OUT=$(realpath -m "$3")
GV_AIE=${AIE:-$GV_ROOT/env/cargo-target-transcript-ends/release/aie}
GV_THREADS=${THREADS:-24}

if [[ ! -s "$GV_OUT/summary.json" ]]; then
  mkdir -p "$(dirname "$GV_OUT")"
  RAYON_NUM_THREADS=$GV_THREADS /usr/bin/time -v -o "$GV_OUT.time.txt" \
    "$GV_AIE" cohort polyasite-mixture \
      --design "$GV_ROOT/gravlax-paper-scripts/experiments/designs/sez-transcript-end-atlas-v1.tsv" \
      --gtf "$GV_ROOT/runs/post-v1/fnbp1-multidonor-sez-validation-r1/gencode.v49.aic" \
      --genome "$GV_ROOT/data/GRCh38.primary_assembly.genome.fa.gz" \
      --polyasite "$GV_ROOT/runs/apa2/polyasite.bed.gz" \
      --group-contrast astro_nsc:mature_neuron --site-gap 24 --out-dir "$GV_OUT" \
      > "$GV_OUT.stdout.json" 2> "$GV_OUT.stderr.txt"
fi

for table in sites.tsv genes.tsv fragment-kernel.tsv; do
  cmp "$GV_FROZEN/$table" "$GV_OUT/$table"
done
sha256sum "$GV_AIE" "$GV_OUT"/{sites.tsv,genes.tsv,fragment-kernel.tsv,summary.json} \
  > "$GV_OUT.execution.sha256"
