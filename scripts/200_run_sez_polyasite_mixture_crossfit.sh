#!/usr/bin/env bash
# Cross-fitted production, catalogue sensitivity, shuffled-label, and legacy timing arms.
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 PROJECT_ROOT OUT_ROOT" >&2
  exit 2
fi

GV_ROOT=$(realpath "$1")
GV_OUT=$(realpath -m "$2")
GV_AIE=${AIE:-$GV_ROOT/env/cargo-target-transcript-ends/release/aie}
GV_LEGACY_AIE=${LEGACY_AIE:-$GV_ROOT/env/cargo-target-features/release/aie}
GV_DESIGN=$GV_ROOT/gravlax-paper-scripts/experiments/designs/sez-transcript-end-atlas-v1.tsv
GV_AIC=$GV_ROOT/runs/post-v1/fnbp1-multidonor-sez-validation-r1/gencode.v49.aic
GV_GENOME=$GV_ROOT/data/GRCh38.primary_assembly.genome.fa.gz
GV_POLYASITE=$GV_ROOT/runs/apa2/polyasite.bed.gz
GV_GROUP_ROOT=$GV_ROOT/gravlax-paper-scripts/results/post-v1-sez-transcript-end-feasibility/groups
GV_ARCHIVE_ROOT=$GV_ROOT/runs/post-v1/sez-transcript-end-atlas-r1
GV_THREADS=${THREADS:-24}
mkdir -p "$GV_OUT" "$GV_OUT/legacy" "$GV_OUT/legacy-groups"

for donor in A B C D E F G H; do
  test -s "$GV_ARCHIVE_ROOT/donor-$donor/sez.called.aie"
  awk -F '\t' '$2 == "astro_nsc" || $2 == "mature_neuron"' \
    "$GV_GROUP_ROOT/donor-$donor.tsv" > "$GV_OUT/legacy-groups/donor-$donor.tsv"
done

run_arm() {
  local label=$1
  local gap=$2
  shift 2
  if [[ -s "$GV_OUT/$label/summary.json" ]]; then
    echo "$label already complete" >&2
    return
  fi
  RAYON_NUM_THREADS=$GV_THREADS /usr/bin/time -v -o "$GV_OUT/$label.time.txt" \
    "$GV_AIE" cohort polyasite-mixture \
      --design "$GV_DESIGN" --gtf "$GV_AIC" --genome "$GV_GENOME" \
      --polyasite "$GV_POLYASITE" --group-contrast astro_nsc:mature_neuron \
      --site-gap "$gap" --out-dir "$GV_OUT/$label" "$@" \
      > "$GV_OUT/$label.stdout.json" 2> "$GV_OUT/$label.stderr.txt"
}

run_arm primary-gap24 24
run_arm sensitivity-gap12 12
run_arm sensitivity-gap48 48
run_arm shuffled-gap24 24 --shuffle-seed 20260901

if [[ ! -s "$GV_OUT/legacy.complete" ]]; then
  export GV_LEGACY_AIE GV_AIC GV_GENOME GV_ARCHIVE_ROOT GV_OUT
  /usr/bin/time -v -o "$GV_OUT/legacy.time.txt" /usr/bin/bash -c '
set -euo pipefail
for donor in A B C D E F G H; do
  "$GV_LEGACY_AIE" query "$GV_ARCHIVE_ROOT/donor-$donor/sez.called.aie" apa-test \
    --gtf "$GV_AIC" --groups "$GV_OUT/legacy-groups/donor-$donor.tsv" \
    --genome "$GV_GENOME" > "$GV_OUT/legacy/donor-$donor.tsv" \
    2> "$GV_OUT/legacy/donor-$donor.stderr.txt"
done
'
  sha256sum "$GV_OUT"/legacy/donor-*.tsv > "$GV_OUT/legacy.complete"
fi

sha256sum "$GV_AIE" "$GV_LEGACY_AIE" "$GV_DESIGN" "$GV_AIC" \
  "$GV_GENOME" "$GV_POLYASITE" > "$GV_OUT/execution-inputs.sha256"
