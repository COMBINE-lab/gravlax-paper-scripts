#!/usr/bin/env bash
# Registered D0 selection and untouched D1/D2-prime hierarchical-EM confirmation.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN}"
: "${PREPARED:?set PREPARED to output from 145_prepare_hierarchical_em_groups.sh}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
manifest=${MANIFEST:-$scripts_repo/experiments/manifests/hierarchical-em.yaml}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
threads=${THREADS:-24}

[[ -x "$AIE_BIN" ]] || { echo "AIE_BIN is not executable: $AIE_BIN" >&2; exit 2; }
[[ -r "$PREPARED/groups/groups.json" ]] || { echo "prepared groups are incomplete" >&2; exit 2; }
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
mkdir -p "$OUT_ROOT/d0-grid" "$OUT_ROOT/selected"
export RAYON_NUM_THREADS=$threads

run_em() {
  local dataset=$1 archive=$2 groups=$3 group_alpha=$4 global_alpha=$5 seed=$6 out=$7 collapse=$8
  mkdir -p "$out"
  local collapse_arg=()
  if [[ "$collapse" == yes ]]; then
    collapse_arg+=(--collapse-groups)
  fi
  /usr/bin/time -v -o "$out/time.txt" \
    "$AIE_BIN" em "$archive" --gtf "$gtf" --mask 0.2 --seed "$seed" --alpha 20 \
      --groups "$groups" --group-alpha "$group_alpha" --global-alpha "$global_alpha" \
      --metrics-json "$out/metrics.json" "${collapse_arg[@]}" \
      >"$out/stdout.txt" 2>"$out/stderr.txt"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$dataset" "$seed" "$group_alpha" "$global_alpha" "$collapse" "$groups" "$out" \
    >>"$OUT_ROOT/runs.tsv"
}

printf 'dataset\tseed\tgroup_alpha\tglobal_alpha\tcollapse\tgroups\tout\n' >"$OUT_ROOT/runs.tsv"
d0_archive=$PROJ_ROOT/runs/archive/d0.aie
d0_groups=$PREPARED/groups/D0.real.tsv
for group_alpha in 5 20 80; do
  for global_alpha in 0 5 20; do
    for seed in 7 17 29; do
      run_em D0 "$d0_archive" "$d0_groups" "$group_alpha" "$global_alpha" "$seed" \
        "$OUT_ROOT/d0-grid/ag-$group_alpha.aa-$global_alpha/seed-$seed" no
    done
  done
done

"$scripts_repo/scripts/146_select_hierarchical_em.py" "$OUT_ROOT/d0-grid" \
  --out "$OUT_ROOT/selection.json"
read -r selected_group_alpha selected_global_alpha < <(
  python3 -c 'import json,sys; x=json.load(open(sys.argv[1]))["selected"]; print(x["group_alpha"], x["global_alpha"])' \
    "$OUT_ROOT/selection.json"
)

for item in D0:d0 D1:d1 D2_prime:d2p; do
  dataset=${item%%:*}
  stem=${item##*:}
  archive=$PROJ_ROOT/runs/archive/$stem.aie
  for control in real shuffled collapse; do
    groups=$PREPARED/groups/$dataset.real.tsv
    collapse=no
    if [[ "$control" == shuffled ]]; then
      groups=$PREPARED/groups/$dataset.shuffled.tsv
    elif [[ "$control" == collapse ]]; then
      collapse=yes
    fi
    for seed in 7 17 29; do
      run_em "$dataset" "$archive" "$groups" "$selected_group_alpha" \
        "$selected_global_alpha" "$seed" \
        "$OUT_ROOT/selected/$dataset/$control/seed-$seed" "$collapse"
    done
  done
done

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-31\n'
  printf 'threads\t%s\n' "$threads"
  printf 'selected_group_alpha\t%s\n' "$selected_group_alpha"
  printf 'selected_global_alpha\t%s\n' "$selected_global_alpha"
  printf 'gravlax_commit\t%s\n' "$(git -C "$PROJ_ROOT/work/gravlax" rev-parse HEAD)"
  printf 'scripts_commit\t%s\n' "$(git -C "$scripts_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$manifest" | cut -d' ' -f1)"
  printf 'prepared_groups_sha256\t%s\n' "$(sha256sum "$PREPARED/groups/groups.json" | cut -d' ' -f1)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

printf 'hierarchical EM gate complete: %s\n' "$OUT_ROOT"
