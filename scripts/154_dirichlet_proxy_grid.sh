#!/usr/bin/env bash
# Locked D0 posterior-mean Dirichlet-proxy screen and matched convex controls.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN to Gravlax 16c012d or a documented descendant}"
: "${PREPARED:?set PREPARED to the leakage-controlled hierarchical group directory}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
manifest=${MANIFEST:-$scripts_repo/experiments/manifests/dirichlet-proxy-screen.yaml}
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
mkdir -p "$OUT_ROOT/grid" "$OUT_ROOT/fixed-convex/real" "$OUT_ROOT/fixed-convex/shuffled"
export RAYON_NUM_THREADS=$threads

printf 'dataset\tseed\tmodel\tcontrol\tcell_prior\tgroup_prior\tgroups\tout\n' \
  >"$OUT_ROOT/runs.tsv"

run_proxy() {
  local cell_prior=$1 group_prior=$2 seed=$3 group_file=$4 control=$5 out=$6
  mkdir -p "$out"
  /usr/bin/time -v -o "$out/time.txt" \
    "$AIE_BIN" em "$archive" --gtf "$gtf" --mask 0.2 --seed "$seed" --alpha 20 \
      --groups "$group_file" --dirichlet-only \
      --dirichlet-cell-prior "$cell_prior" \
      --dirichlet-group-prior "$group_prior" \
      --metrics-json "$out/metrics.json" \
      >"$out/stdout.txt" 2>"$out/stderr.txt"
  printf 'D0\t%s\tdirichlet-proxy\t%s\t%s\t%s\t%s\t%s\n' \
    "$seed" "$control" "$cell_prior" "$group_prior" "$group_file" "$out" \
    >>"$OUT_ROOT/runs.tsv"
}

run_convex() {
  local seed=$1 group_file=$2 control=$3 out=$4
  mkdir -p "$out"
  /usr/bin/time -v -o "$out/time.txt" \
    "$AIE_BIN" em "$archive" --gtf "$gtf" --mask 0.2 --seed "$seed" --alpha 20 \
      --groups "$group_file" --convex-only \
      --convex-cell-weight 0.20 --convex-group-weight 0.60 --convex-group-prior 80 \
      --metrics-json "$out/metrics.json" \
      >"$out/stdout.txt" 2>"$out/stderr.txt"
  printf 'D0\t%s\tfixed-convex\t%s\tNA\tNA\t%s\t%s\n' \
    "$seed" "$control" "$group_file" "$out" >>"$OUT_ROOT/runs.tsv"
}

for cell_prior in 1 4 16 64 256; do
  for group_prior in 1 5 20 80 320; do
    label="kh-$cell_prior.k0-$group_prior"
    for seed in 7 17 29; do
      run_proxy "$cell_prior" "$group_prior" "$seed" "$group_map" real \
        "$OUT_ROOT/grid/$label/seed-$seed"
    done
  done
done

for seed in 7 17 29; do
  run_convex "$seed" "$group_map" real "$OUT_ROOT/fixed-convex/real/seed-$seed"
  run_convex "$seed" "$shuffled_group_map" shuffled \
    "$OUT_ROOT/fixed-convex/shuffled/seed-$seed"
done

"$scripts_repo/scripts/155_select_dirichlet_proxy.py" \
  "$OUT_ROOT" --out "$OUT_ROOT/selection.json"

read -r selected_cell_prior selected_group_prior < <(
  python3 -c 'import json,sys; x=json.load(open(sys.argv[1]))["selected"]; print(x["cell_prior"], x["group_prior"])' \
    "$OUT_ROOT/selection.json"
)
for seed in 7 17 29; do
  run_proxy "$selected_cell_prior" "$selected_group_prior" "$seed" \
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

"$scripts_repo/scripts/156_summarize_dirichlet_proxy.py" \
  "$OUT_ROOT" --out "$OUT_ROOT/summary.json"
printf 'Dirichlet-proxy D0 screen complete: %s\n' "$OUT_ROOT"
