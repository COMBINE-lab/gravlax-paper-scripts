#!/usr/bin/env bash
# Run the locked hybrid/proxy factorial and fixed-convex depth guardrail.
set -euo pipefail

dataset=${1:-}
case "$dataset" in D4|D3) ;; *) echo "usage: $0 D4|D3" >&2; exit 2 ;; esac
: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN}"
: "${PYTHON_BIN:?set PYTHON_BIN to an environment with PyYAML}"
: "${PREPARED:?set PREPARED to output from script 160}"
: "${OUT_ROOT:?set OUT_ROOT to a new evaluation directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
manifest=${MANIFEST:-$scripts_repo/experiments/manifests/depth-hybrid-em-confirmation.yaml}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
archive=$PROJ_ROOT/runs/archive/${dataset,,}.aie
real_groups=$PREPARED/groups/$dataset.real.tsv
shuffled_groups=$PREPARED/groups/$dataset.shuffled.tsv
threads=${THREADS:-24}

for input in "$AIE_BIN" "$PYTHON_BIN" "$manifest" "$gtf" "$archive" "$real_groups" "$shuffled_groups"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
mkdir -p "$OUT_ROOT"
export RAYON_NUM_THREADS=$threads
printf 'dataset\tseed\tmodel\tcontrol\tgroups\tout\n' >"$OUT_ROOT/runs.tsv"

run_mode() {
  local model=$1 control=$2 seed=$3 groups=$4 out=$5
  local mode_args=()
  case "$model" in
    depth-hybrid) mode_args=(--hybrid-only --hybrid-depth-scale 64) ;;
    dirichlet-proxy) mode_args=(--dirichlet-only --dirichlet-cell-prior 64 --dirichlet-group-prior 80) ;;
    fixed-convex) mode_args=(--convex-only --convex-cell-weight 0.20 --convex-group-weight 0.60 --convex-group-prior 80) ;;
    *) echo "unknown model: $model" >&2; exit 2 ;;
  esac
  mkdir -p "$out"
  /usr/bin/time -v -o "$out/time.txt" \
    "$AIE_BIN" em "$archive" --gtf "$gtf" --mask 0.2 --seed "$seed" --alpha 20 \
      --groups "$groups" "${mode_args[@]}" --metrics-json "$out/metrics.json" \
      >"$out/stdout.txt" 2>"$out/stderr.txt"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$dataset" "$seed" "$model" "$control" "$groups" "$out" >>"$OUT_ROOT/runs.tsv"
}

for model in depth-hybrid dirichlet-proxy; do
  for control in real shuffled; do
    groups=$real_groups
    [[ "$control" == shuffled ]] && groups=$shuffled_groups
    for seed in 7 17 29; do
      run_mode "$model" "$control" "$seed" "$groups" \
        "$OUT_ROOT/$model/$control/seed-$seed"
    done
  done
done
for seed in 7 17 29; do
  run_mode fixed-convex real "$seed" "$real_groups" "$OUT_ROOT/fixed-convex/real/seed-$seed"
done

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-31\n'
  printf 'dataset\t%s\n' "$dataset"
  printf 'threads\t%s\n' "$threads"
  printf 'gravlax_commit\t%s\n' "$(git -C "$PROJ_ROOT/work/gravlax" rev-parse HEAD)"
  printf 'scripts_commit\t%s\n' "$(git -C "$scripts_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$manifest" | cut -d' ' -f1)"
  printf 'archive_sha256\t%s\n' "$(sha256sum "$archive" | cut -d' ' -f1)"
  printf 'gtf_sha256\t%s\n' "$(sha256sum "$gtf" | cut -d' ' -f1)"
  printf 'groups_sha256\t%s\n' "$(sha256sum "$real_groups" | cut -d' ' -f1)"
  printf 'shuffled_groups_sha256\t%s\n' "$(sha256sum "$shuffled_groups" | cut -d' ' -f1)"
  printf 'group_selection_sha256\t%s\n' "$(sha256sum "$PREPARED/groups/groups.json" | cut -d' ' -f1)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

"$PYTHON_BIN" "$scripts_repo/scripts/163_summarize_em_confirmation.py" "$OUT_ROOT" \
  --manifest "$manifest" --out "$OUT_ROOT/summary.json"
printf '%s EM confirmation complete: %s\n' "$dataset" "$OUT_ROOT"
