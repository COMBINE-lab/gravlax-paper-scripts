#!/usr/bin/env bash
# D0-only development grid for candidate-normalized convex partial pooling.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN to Gravlax f6fd759 or a documented descendant}"
: "${PREPARED:?set PREPARED to the leakage-controlled hierarchical group directory}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
manifest=${MANIFEST:-$scripts_repo/experiments/manifests/candidate-normalized-convex-em.yaml}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
archive=${ARCHIVE:-$PROJ_ROOT/runs/archive/d0.aie}
group_map=${EM_GROUPS:-$PREPARED/groups/D0.real.tsv}
shuffled_group_map=${SHUFFLED_EM_GROUPS:-$PREPARED/groups/D0.shuffled.tsv}
threads=${THREADS:-24}

[[ -x "$AIE_BIN" ]] || { echo "AIE_BIN is not executable: $AIE_BIN" >&2; exit 2; }
[[ -r "$manifest" ]] || { echo "manifest is not readable: $manifest" >&2; exit 2; }
[[ -r "$archive" ]] || { echo "archive is not readable: $archive" >&2; exit 2; }
[[ -r "$group_map" ]] || { echo "groups are not readable: $group_map" >&2; exit 2; }
[[ -r "$shuffled_group_map" ]] || {
  echo "shuffled groups are not readable: $shuffled_group_map" >&2
  exit 2
}
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
mkdir -p "$OUT_ROOT/grid"
export RAYON_NUM_THREADS=$threads

run_em() {
  local cell_weight=$1 group_weight=$2 group_prior=$3 seed=$4 group_file=$5 control=$6 out=$7
  mkdir -p "$out"
  /usr/bin/time -v -o "$out/time.txt" \
    "$AIE_BIN" em "$archive" --gtf "$gtf" --mask 0.2 --seed "$seed" --alpha 20 \
      --groups "$group_file" --convex-only \
      --convex-cell-weight "$cell_weight" \
      --convex-group-weight "$group_weight" \
      --convex-group-prior "$group_prior" \
      --metrics-json "$out/metrics.json" \
      >"$out/stdout.txt" 2>"$out/stderr.txt"
  printf 'D0\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$seed" "$cell_weight" "$group_weight" "$group_prior" "$control" "$group_file" "$out" \
    >>"$OUT_ROOT/runs.tsv"
}

printf 'dataset\tseed\tcell_weight\tgroup_weight\tgroup_prior\tcontrol\tgroups\tout\n' \
  >"$OUT_ROOT/runs.tsv"

# Explicit decimals implement group_weight=(1-cell_weight)*rho without depending on bc/awk
# floating-point behavior. Group-zero configurations are deduplicated because kappa is irrelevant.
for config in \
  0.00:0.0000:20 0.00:0.2500:0 0.00:0.2500:5 0.00:0.2500:20 0.00:0.2500:80 \
  0.00:0.5000:0 0.00:0.5000:5 0.00:0.5000:20 0.00:0.5000:80 \
  0.00:0.7500:0 0.00:0.7500:5 0.00:0.7500:20 0.00:0.7500:80 \
  0.00:1.0000:0 0.00:1.0000:5 0.00:1.0000:20 0.00:1.0000:80 \
  0.05:0.0000:20 0.05:0.2375:0 0.05:0.2375:5 0.05:0.2375:20 0.05:0.2375:80 \
  0.05:0.4750:0 0.05:0.4750:5 0.05:0.4750:20 0.05:0.4750:80 \
  0.05:0.7125:0 0.05:0.7125:5 0.05:0.7125:20 0.05:0.7125:80 \
  0.05:0.9500:0 0.05:0.9500:5 0.05:0.9500:20 0.05:0.9500:80 \
  0.10:0.0000:20 0.10:0.2250:0 0.10:0.2250:5 0.10:0.2250:20 0.10:0.2250:80 \
  0.10:0.4500:0 0.10:0.4500:5 0.10:0.4500:20 0.10:0.4500:80 \
  0.10:0.6750:0 0.10:0.6750:5 0.10:0.6750:20 0.10:0.6750:80 \
  0.10:0.9000:0 0.10:0.9000:5 0.10:0.9000:20 0.10:0.9000:80 \
  0.20:0.0000:20 0.20:0.2000:0 0.20:0.2000:5 0.20:0.2000:20 0.20:0.2000:80 \
  0.20:0.4000:0 0.20:0.4000:5 0.20:0.4000:20 0.20:0.4000:80 \
  0.20:0.6000:0 0.20:0.6000:5 0.20:0.6000:20 0.20:0.6000:80 \
  0.20:0.8000:0 0.20:0.8000:5 0.20:0.8000:20 0.20:0.8000:80
do
  cell_weight=${config%%:*}
  rest=${config#*:}
  group_weight=${rest%%:*}
  group_prior=${config##*:}
  label="cw-$cell_weight.gw-$group_weight.k-$group_prior"
  for seed in 7 17 29; do
    run_em "$cell_weight" "$group_weight" "$group_prior" "$seed" \
      "$group_map" real "$OUT_ROOT/grid/$label/seed-$seed"
  done
done

"$scripts_repo/scripts/152_select_candidate_normalized_convex.py" \
  "$OUT_ROOT" --out "$OUT_ROOT/selection.json"

read -r selected_cell_weight selected_group_weight selected_group_prior < <(
  python3 -c 'import json,sys; x=json.load(open(sys.argv[1]))["selected"]; print(x["cell_weight"], x["group_weight"], x["group_prior"])' \
    "$OUT_ROOT/selection.json"
)
for seed in 7 17 29; do
  run_em "$selected_cell_weight" "$selected_group_weight" "$selected_group_prior" "$seed" \
    "$shuffled_group_map" shuffled "$OUT_ROOT/selected-shuffled/seed-$seed"
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

printf 'candidate-normalized convex D0 grid complete: %s\n' "$OUT_ROOT"
