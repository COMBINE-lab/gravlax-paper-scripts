#!/usr/bin/env bash
# Locked D0 grid for the monotone evidence-depth EM hybrid.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN to Gravlax 9ccc3f6 or a documented descendant}"
: "${PREPARED:?set PREPARED to the leakage-controlled hierarchical group directory}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
manifest=${MANIFEST:-$scripts_repo/experiments/manifests/depth-gated-hybrid-em.yaml}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
archive=${ARCHIVE:-$PROJ_ROOT/runs/archive/d0.aie}
group_map=${EM_GROUPS:-$PREPARED/groups/D0.real.tsv}
shuffled_group_map=${SHUFFLED_EM_GROUPS:-$PREPARED/groups/D0.shuffled.tsv}
threads=${THREADS:-24}

[[ -x "$AIE_BIN" ]] || { echo "AIE_BIN is not executable: $AIE_BIN" >&2; exit 2; }
[[ -r "$manifest" ]] || { echo "manifest is not readable: $manifest" >&2; exit 2; }
[[ -r "$archive" ]] || { echo "archive is not readable: $archive" >&2; exit 2; }
[[ -r "$group_map" ]] || { echo "groups are not readable: $group_map" >&2; exit 2; }
[[ -r "$shuffled_group_map" ]] || { echo "shuffled groups are not readable" >&2; exit 2; }
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
mkdir -p "$OUT_ROOT/grid"
export RAYON_NUM_THREADS=$threads

printf 'dataset\tseed\tmodel\tcontrol\tdepth_scale\tgroups\tout\n' >"$OUT_ROOT/runs.tsv"

run_mode() {
  local model=$1 control=$2 seed=$3 group_file=$4 depth_scale=$5 out=$6
  mkdir -p "$out"
  local mode_args=()
  case "$model" in
    depth-hybrid)
      mode_args=(--hybrid-only --hybrid-depth-scale "$depth_scale")
      ;;
    fixed-convex)
      mode_args=(--convex-only --convex-cell-weight 0.20 --convex-group-weight 0.60 --convex-group-prior 80)
      ;;
    dirichlet-proxy)
      mode_args=(--dirichlet-only --dirichlet-cell-prior 64 --dirichlet-group-prior 80)
      ;;
    *) echo "unknown model: $model" >&2; exit 2 ;;
  esac
  /usr/bin/time -v -o "$out/time.txt" \
    "$AIE_BIN" em "$archive" --gtf "$gtf" --mask 0.2 --seed "$seed" --alpha 20 \
      --groups "$group_file" "${mode_args[@]}" --metrics-json "$out/metrics.json" \
      >"$out/stdout.txt" 2>"$out/stderr.txt"
  printf 'D0\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$seed" "$model" "$control" "$depth_scale" "$group_file" "$out" >>"$OUT_ROOT/runs.tsv"
}

for depth_scale in 2 4 8 16 32 64 128 256 512 1024; do
  for seed in 7 17 29; do
    run_mode depth-hybrid real "$seed" "$group_map" "$depth_scale" \
      "$OUT_ROOT/grid/scale-$depth_scale/seed-$seed"
  done
done

for model in fixed-convex dirichlet-proxy; do
  for control in real shuffled; do
    group_file=$group_map
    [[ "$control" == shuffled ]] && group_file=$shuffled_group_map
    for seed in 7 17 29; do
      run_mode "$model" "$control" "$seed" "$group_file" NA \
        "$OUT_ROOT/constituents/$model/$control/seed-$seed"
    done
  done
done

"$scripts_repo/scripts/158_select_depth_hybrid.py" "$OUT_ROOT" --out "$OUT_ROOT/selection.json"
selected_scale=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["selected"]["depth_scale"])' "$OUT_ROOT/selection.json")
for seed in 7 17 29; do
  run_mode depth-hybrid shuffled "$seed" "$shuffled_group_map" "$selected_scale" \
    "$OUT_ROOT/selected-shuffled/seed-$seed"
done

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-31\n'
  printf 'scope\tD0_development_only\n'
  printf 'threads\t%s\n' "$threads"
  printf 'gravlax_commit\t%s\n' "$(git -C "$PROJ_ROOT/work/gravlax" rev-parse HEAD)"
  printf 'scripts_commit\t%s\n' "$(git -C "$scripts_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$manifest" | cut -d' ' -f1)"
  printf 'groups_sha256\t%s\n' "$(sha256sum "$group_map" | cut -d' ' -f1)"
  printf 'shuffled_groups_sha256\t%s\n' "$(sha256sum "$shuffled_group_map" | cut -d' ' -f1)"
  printf 'archive_sha256\t%s\n' "$(sha256sum "$archive" | cut -d' ' -f1)"
  printf 'gtf_sha256\t%s\n' "$(sha256sum "$gtf" | cut -d' ' -f1)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

"$scripts_repo/scripts/159_summarize_depth_hybrid.py" "$OUT_ROOT" --out "$OUT_ROOT/summary.json"
printf 'depth-hybrid D0 gate complete: %s\n' "$OUT_ROOT"
