#!/usr/bin/env bash
# Locked D0/D4 promotion grid plus full-mode baselines for paired diagnostics.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN to Gravlax 0446fea or a documented descendant}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
manifest=${MANIFEST:-$scripts_repo/experiments/manifests/hybrid-em-promotion.yaml}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
threads=${THREADS:-24}
d0_groups=${D0_GROUPS:-$PROJ_ROOT/runs/post-v1/hierarchical-em-groups-r1/groups/D0.real.tsv}
d4_groups=${D4_GROUPS:-$PROJ_ROOT/runs/post-v1/em-confirmation-d4-prep-r1/groups/D4.real.tsv}

for input in "$AIE_BIN" "$manifest" "$gtf" "$d0_groups" "$d4_groups" \
  "$PROJ_ROOT/runs/archive/d0.aie" "$PROJ_ROOT/runs/archive/d4.aie"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
mkdir -p "$OUT_ROOT"
export RAYON_NUM_THREADS=$threads
printf 'dataset\tseed\trun_type\tdepth_scale\tdepth_power\tgroups\tout\n' >"$OUT_ROOT/runs.tsv"

run_one() {
  local dataset=$1 seed=$2 run_type=$3 scale=$4 power=$5 groups=$6 out=$7
  local archive=$PROJ_ROOT/runs/archive/${dataset,,}.aie
  local mode_args=()
  if [[ "$run_type" == grid ]]; then
    mode_args=(--hybrid-only --hybrid-depth-scale "$scale" --hybrid-depth-power "$power")
  else
    mode_args=(--hybrid-depth-scale "$scale" --hybrid-depth-power "$power")
  fi
  mkdir -p "$out"
  /usr/bin/time -v -o "$out/time.txt" \
    "$AIE_BIN" em "$archive" --gtf "$gtf" --mask 0.2 --seed "$seed" --alpha 20 \
      --groups "$groups" \
      --convex-cell-weight 0.20 --convex-group-weight 0.60 --convex-group-prior 80 \
      --dirichlet-cell-prior 64 --dirichlet-group-prior 80 \
      "${mode_args[@]}" --metrics-json "$out/metrics.json" \
      >"$out/stdout.txt" 2>"$out/stderr.txt"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$dataset" "$seed" "$run_type" "$scale" "$power" "$groups" "$out" >>"$OUT_ROOT/runs.tsv"
}

for dataset in D0 D4; do
  groups=$d0_groups
  [[ "$dataset" == D4 ]] && groups=$d4_groups
  for seed in 7 17 29; do
    run_one "$dataset" "$seed" baseline 64 2 "$groups" \
      "$OUT_ROOT/baseline/$dataset/seed-$seed"
  done
  for power in 0.5 1 2 4 8; do
    power_dir=${power/./p}
    for scale in 2 4 8 16 32 64 128 256; do
      for seed in 7 17 29; do
        run_one "$dataset" "$seed" grid "$scale" "$power" "$groups" \
          "$OUT_ROOT/grid/$dataset/power-$power_dir/scale-$scale/seed-$seed"
      done
    done
  done
done

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-31\n'
  printf 'scope\tD0_D4_development\n'
  printf 'threads\t%s\n' "$threads"
  printf 'gravlax_commit\t%s\n' "$(git -C "$PROJ_ROOT/work/gravlax" rev-parse HEAD)"
  printf 'scripts_commit\t%s\n' "$(git -C "$scripts_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$manifest" | cut -d' ' -f1)"
  printf 'archive_D0_sha256\t%s\n' "$(sha256sum "$PROJ_ROOT/runs/archive/d0.aie" | cut -d' ' -f1)"
  printf 'archive_D4_sha256\t%s\n' "$(sha256sum "$PROJ_ROOT/runs/archive/d4.aie" | cut -d' ' -f1)"
  printf 'gtf_sha256\t%s\n' "$(sha256sum "$gtf" | cut -d' ' -f1)"
  printf 'groups_D0_sha256\t%s\n' "$(sha256sum "$d0_groups" | cut -d' ' -f1)"
  printf 'groups_D4_sha256\t%s\n' "$(sha256sum "$d4_groups" | cut -d' ' -f1)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

python3 "$scripts_repo/scripts/165_select_hybrid_promotion.py" "$OUT_ROOT" \
  --out "$OUT_ROOT/selection.json"
printf 'hybrid promotion development complete: %s\n' "$OUT_ROOT"

