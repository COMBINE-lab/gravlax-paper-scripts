#!/usr/bin/env bash
# Rerun the selected D0/D4 candidate with all reference modes for paired diagnostics.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN}"
: "${DEVELOPMENT_ROOT:?set DEVELOPMENT_ROOT to output from script 164}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
manifest=${MANIFEST:-$scripts_repo/experiments/manifests/hybrid-em-promotion.yaml}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
threads=${THREADS:-24}
d0_groups=${D0_GROUPS:-$PROJ_ROOT/runs/post-v1/hierarchical-em-groups-r1/groups/D0.real.tsv}
d4_groups=${D4_GROUPS:-$PROJ_ROOT/runs/post-v1/em-confirmation-d4-prep-r1/groups/D4.real.tsv}
selection=$DEVELOPMENT_ROOT/selection.json

for input in "$AIE_BIN" "$manifest" "$gtf" "$selection" "$d0_groups" "$d4_groups"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
read -r scale power < <(python3 -c '
import json, sys
x=json.load(open(sys.argv[1]))
if x["status"] != "ELIGIBLE_MODEL_SELECTED":
    raise SystemExit("no eligible model; selected audit is forbidden")
print(x["selected"]["depth_scale"], x["selected"]["depth_power"])
' "$selection")

mkdir -p "$OUT_ROOT"
export RAYON_NUM_THREADS=$threads
printf 'dataset\tseed\tdepth_scale\tdepth_power\tgroups\tout\n' >"$OUT_ROOT/runs.tsv"
for dataset in D0 D4; do
  groups=$d0_groups
  [[ "$dataset" == D4 ]] && groups=$d4_groups
  archive=$PROJ_ROOT/runs/archive/${dataset,,}.aie
  for seed in 7 17 29; do
    out=$OUT_ROOT/$dataset/seed-$seed
    mkdir -p "$out"
    /usr/bin/time -v -o "$out/time.txt" \
      "$AIE_BIN" em "$archive" --gtf "$gtf" --mask 0.2 --seed "$seed" --alpha 20 \
        --groups "$groups" \
        --convex-cell-weight 0.20 --convex-group-weight 0.60 --convex-group-prior 80 \
        --dirichlet-cell-prior 64 --dirichlet-group-prior 80 \
        --hybrid-depth-scale "$scale" --hybrid-depth-power "$power" \
        --metrics-json "$out/metrics.json" >"$out/stdout.txt" 2>"$out/stderr.txt"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$dataset" "$seed" "$scale" "$power" "$groups" "$out" >>"$OUT_ROOT/runs.tsv"
  done
done

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-31\n'
  printf 'scope\tselected_D0_D4_development_audit\n'
  printf 'depth_scale\t%s\n' "$scale"
  printf 'depth_power\t%s\n' "$power"
  printf 'threads\t%s\n' "$threads"
  printf 'gravlax_commit\t%s\n' "$(git -C "$PROJ_ROOT/work/gravlax" rev-parse HEAD)"
  printf 'scripts_commit\t%s\n' "$(git -C "$scripts_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$manifest" | cut -d' ' -f1)"
  printf 'selection_sha256\t%s\n' "$(sha256sum "$selection" | cut -d' ' -f1)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

python3 "$scripts_repo/scripts/167_summarize_hybrid_promotion_selected.py" "$OUT_ROOT" \
  --selection "$selection" --out "$OUT_ROOT/summary.json"
printf 'selected hybrid paired audit complete: %s\n' "$OUT_ROOT"

