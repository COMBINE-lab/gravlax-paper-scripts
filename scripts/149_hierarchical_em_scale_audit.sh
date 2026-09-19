#!/usr/bin/env bash
# Post-hoc scale audit for the registered hierarchical-EM STOP result.
#
# This deliberately does not define a new confirmatory gate: D1 and D2-prime were inspected after
# the registered grid failed.  It diagnoses the effective pseudo-count scale and records the
# exploratory result that motivates a separately validated, candidate-normalized model.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN}"
: "${PREPARED:?set PREPARED to output from 145_prepare_hierarchical_em_groups.sh}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
threads=${THREADS:-24}

[[ -x "$AIE_BIN" ]] || { echo "AIE_BIN is not executable: $AIE_BIN" >&2; exit 2; }
[[ -r "$PREPARED/groups/groups.json" ]] || { echo "prepared groups are incomplete" >&2; exit 2; }
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
mkdir -p "$OUT_ROOT/d0-grid" "$OUT_ROOT/confirm" "$OUT_ROOT/shuffled"
export RAYON_NUM_THREADS=$threads

run_em() {
  local archive=$1 groups=$2 group_alpha=$3 global_alpha=$4 seed=$5 out=$6
  mkdir -p "$out"
  "$AIE_BIN" em "$archive" --gtf "$gtf" --mask 0.2 --seed "$seed" --alpha 20 \
    --groups "$groups" --group-alpha "$group_alpha" --global-alpha "$global_alpha" \
    --metrics-json "$out/metrics.json" >"$out/stdout.txt" 2>"$out/stderr.txt"
}

# Total prior masses span the point at which a unit local count ceases to dominate ordinary
# candidate genes.  Fractions span calibrated global-only through increasingly group-heavy priors.
for pair in \
  0:6400 1280:5120 3200:3200 5120:1280 \
  0:25600 5120:20480 12800:12800 20480:5120 \
  0:102400 20480:81920 51200:51200 81920:20480 \
  327680:81920
do
  group_alpha=${pair%%:*}
  global_alpha=${pair##*:}
  for seed in 7 17 29; do
    run_em "$PROJ_ROOT/runs/archive/d0.aie" "$PREPARED/groups/D0.real.tsv" \
      "$group_alpha" "$global_alpha" "$seed" \
      "$OUT_ROOT/d0-grid/ag-$group_alpha.aa-$global_alpha/seed-$seed"
  done
done

# D0 alone selected 12,800/12,800 in the exploratory grid.  The 0/25,600 arm separates calibrated
# local/global shrinkage from the incremental value of real groups.
for item in D1:d1 D2_prime:d2p; do
  dataset=${item%%:*}
  stem=${item##*:}
  for seed in 7 17 29; do
    run_em "$PROJ_ROOT/runs/archive/$stem.aie" "$PREPARED/groups/$dataset.real.tsv" \
      12800 12800 "$seed" "$OUT_ROOT/confirm/$dataset/hierarchical/seed-$seed"
    run_em "$PROJ_ROOT/runs/archive/$stem.aie" "$PREPARED/groups/$dataset.real.tsv" \
      0 25600 "$seed" "$OUT_ROOT/confirm/$dataset/calibrated_global/seed-$seed"
  done
done

for item in D0:d0 D1:d1 D2_prime:d2p; do
  dataset=${item%%:*}
  stem=${item##*:}
  for seed in 7 17 29; do
    run_em "$PROJ_ROOT/runs/archive/$stem.aie" "$PREPARED/groups/$dataset.shuffled.tsv" \
      12800 12800 "$seed" "$OUT_ROOT/shuffled/$dataset/seed-$seed"
  done
done

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-31\n'
  printf 'status\tpost_hoc_exploratory_not_independent_confirmation\n'
  printf 'threads\t%s\n' "$threads"
  printf 'gravlax_commit\t%s\n' "$(git -C "$PROJ_ROOT/work/gravlax" rev-parse HEAD)"
  printf 'scripts_commit\t%s\n' "$(git -C "$scripts_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

"$scripts_repo/scripts/150_summarize_hierarchical_em_scale_audit.py" \
  "$OUT_ROOT" --project-root "$PROJ_ROOT" --out "$OUT_ROOT/summary.json"
printf 'exploratory hierarchical-EM scale audit complete: %s\n' "$OUT_ROOT"
