#!/usr/bin/env bash
# Run the single untouched D3 confirmation of the frozen robust hybrid candidate.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN}"
: "${PREPARED:?set PREPARED to the frozen D3 preparation from script 160}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
manifest=${MANIFEST:-$scripts_repo/experiments/manifests/hybrid-em-promotion.yaml}
selection=${SELECTION:-$PROJ_ROOT/runs/post-v1/hybrid-em-promotion-robust-selection-r1/selection.json}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
archive=${ARCHIVE:-$PROJ_ROOT/runs/archive/d3.aie}
real_groups=$PREPARED/groups/D3.real.tsv
shuffled_groups=$PREPARED/groups/D3.shuffled.tsv
threads=${THREADS:-24}

for input in "$AIE_BIN" "$manifest" "$selection" "$gtf" "$archive" \
  "$real_groups" "$shuffled_groups" "$PREPARED/groups/groups.json"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
read -r scale power < <(python3 -c '
import json, sys
x=json.load(open(sys.argv[1]))
if x["status"] != "ELIGIBLE_MODEL_SELECTED": raise SystemExit("selection is not eligible")
print(x["selected"]["depth_scale"], x["selected"]["depth_power"])
' "$selection")
[[ "$scale" == 8 && "$power" == 8 ]] || {
  echo "selection does not match frozen D3 candidate: scale=$scale power=$power" >&2; exit 2;
}

mkdir -p "$OUT_ROOT"
export RAYON_NUM_THREADS=$threads
printf 'dataset\tseed\tcontrol\tdepth_scale\tdepth_power\tgroups\tout\n' >"$OUT_ROOT/runs.tsv"
for control in real shuffled; do
  groups=$real_groups
  mode_args=()
  [[ "$control" == shuffled ]] && { groups=$shuffled_groups; mode_args=(--hybrid-only); }
  for seed in 7 17 29; do
    out=$OUT_ROOT/$control/seed-$seed
    mkdir -p "$out"
    /usr/bin/time -v -o "$out/time.txt" \
      "$AIE_BIN" em "$archive" --gtf "$gtf" --mask 0.2 --seed "$seed" --alpha 20 \
        --groups "$groups" \
        --convex-cell-weight 0.20 --convex-group-weight 0.60 --convex-group-prior 80 \
        --dirichlet-cell-prior 64 --dirichlet-group-prior 80 \
        --hybrid-depth-scale "$scale" --hybrid-depth-power "$power" \
        "${mode_args[@]}" --metrics-json "$out/metrics.json" \
        >"$out/stdout.txt" 2>"$out/stderr.txt"
    printf 'D3\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$seed" "$control" "$scale" "$power" "$groups" "$out" >>"$OUT_ROOT/runs.tsv"
  done
done

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-31\n'
  printf 'scope\tuntouched_D3_confirmation\n'
  printf 'depth_scale\t%s\n' "$scale"
  printf 'depth_power\t%s\n' "$power"
  printf 'threads\t%s\n' "$threads"
  printf 'gravlax_commit\t%s\n' "$(git -C "$PROJ_ROOT/work/gravlax" rev-parse HEAD)"
  printf 'scripts_commit\t%s\n' "$(git -C "$scripts_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$manifest" | cut -d' ' -f1)"
  printf 'selection_sha256\t%s\n' "$(sha256sum "$selection" | cut -d' ' -f1)"
  printf 'preparation_protocol_sha256\t%s\n' "$(sha256sum "$PREPARED/protocol.tsv" | cut -d' ' -f1)"
  printf 'group_selection_sha256\t%s\n' "$(sha256sum "$PREPARED/groups/groups.json" | cut -d' ' -f1)"
  printf 'real_groups_sha256\t%s\n' "$(sha256sum "$real_groups" | cut -d' ' -f1)"
  printf 'shuffled_groups_sha256\t%s\n' "$(sha256sum "$shuffled_groups" | cut -d' ' -f1)"
  printf 'archive_sha256\t%s\n' "$(sha256sum "$archive" | cut -d' ' -f1)"
  printf 'gtf_sha256\t%s\n' "$(sha256sum "$gtf" | cut -d' ' -f1)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

python3 "$scripts_repo/scripts/170_summarize_hybrid_promotion_d3.py" "$OUT_ROOT" \
  --out "$OUT_ROOT/summary.json"
printf 'D3 hybrid promotion confirmation complete: %s\n' "$OUT_ROOT"

