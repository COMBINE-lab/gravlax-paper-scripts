#!/usr/bin/env bash
# Reproduce the locked cell/group scope and exact junction-set gates.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT to the annotation-independent-evidence project root}"
: "${AIE_BIN:?set AIE_BIN to Gravlax commit f1f99ece33f3d89b6f54d160cb081b3e49e5e3a7}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
d0_archive=${D0_ARCHIVE:-$PROJ_ROOT/runs/archive/d0.aie}
d1_archive=${D1_ARCHIVE:-$PROJ_ROOT/runs/archive/d1.aie}
d0_group_path=${D0_GROUP_PATH:-$PROJ_ROOT/runs/post-v1/hierarchical-em-groups-r1/groups/D0.real.tsv}
d1_group_path=${D1_GROUP_PATH:-$PROJ_ROOT/runs/post-v1/hierarchical-em-groups-r1/groups/D1.real.tsv}
plan=${PLAN:-$script_dir/../experiments/plans/junction-set-d0.tsv}
threads=${THREADS:-24}
commit=${AIE_COMMIT:-f1f99ece33f3d89b6f54d160cb081b3e49e5e3a7}
[[ "$threads" == 24 ]] || {
  echo "the locked gate requires THREADS=24 (received $threads)" >&2
  exit 2
}

for input in "$AIE_BIN" "$d0_archive" "$d1_archive" \
  "$d0_group_path" "$d1_group_path" "$plan"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
mkdir -p "$OUT_ROOT"/{d0,d1,performance,compatibility}

include_locus=$(awk -F '\t' '$1 == "include" {print $2}' "$plan")
exclude_locus=$(awk -F '\t' '$1 == "exclude" {print $2}' "$plan")
[[ -n "$include_locus" && -n "$exclude_locus" ]] || {
  echo "plan must define one include and one exclude locus" >&2
  exit 2
}
export RAYON_NUM_THREADS=$threads
awk -F '\t' '{print $1}' "$d0_group_path" >"$OUT_ROOT/d0/cells.txt"

run_jset() {
  "$AIE_BIN" query "$1" jset \
    --include "$include_locus" --exclude "$exclude_locus" "${@:2}"
}

run_jset "$d0_archive" --groups "$d0_group_path" --json \
  >"$OUT_ROOT/d0/groups.json" 2>"$OUT_ROOT/d0/groups.stderr"
run_jset "$d0_archive" --groups "$d0_group_path" --agg cell --top 0 --json \
  >"$OUT_ROOT/d0/cells-from-groups.json" 2>"$OUT_ROOT/d0/cells-from-groups.stderr"
run_jset "$d0_archive" --cells "$OUT_ROOT/d0/cells.txt" --agg cell --top 0 --json \
  >"$OUT_ROOT/d0/cells-from-list.json" 2>"$OUT_ROOT/d0/cells-from-list.stderr"
run_jset "$d0_archive" --groups "$d0_group_path" --agg bulk --json \
  >"$OUT_ROOT/d0/bulk.json" 2>"$OUT_ROOT/d0/bulk.stderr"
run_jset "$d0_archive" --groups "$d0_group_path" --json \
  >"$OUT_ROOT/d0/groups-repeat.json" 2>"$OUT_ROOT/d0/groups-repeat.stderr"
"$AIE_BIN" query "$d0_archive" junction "$include_locus" \
  --groups "$d0_group_path" --agg cell --top 0 --json \
  >"$OUT_ROOT/d0/include-point.json" 2>"$OUT_ROOT/d0/include-point.stderr"
"$AIE_BIN" query "$d0_archive" junction "$exclude_locus" \
  --groups "$d0_group_path" --agg cell --top 0 --json \
  >"$OUT_ROOT/d0/exclude-point.json" 2>"$OUT_ROOT/d0/exclude-point.stderr"

run_jset "$d1_archive" --groups "$d1_group_path" --json \
  >"$OUT_ROOT/d1/groups.json" 2>"$OUT_ROOT/d1/groups.stderr"
run_jset "$d1_archive" --groups "$d1_group_path" --agg cell --top 0 --json \
  >"$OUT_ROOT/d1/cells.json" 2>"$OUT_ROOT/d1/cells.stderr"
run_jset "$d1_archive" --groups "$d1_group_path" --agg bulk --json \
  >"$OUT_ROOT/d1/bulk.json" 2>"$OUT_ROOT/d1/bulk.stderr"

run_jset "$d0_archive" --json >/dev/null 2>/dev/null
: >"$OUT_ROOT/performance/jset-unscoped.tsv"
: >"$OUT_ROOT/performance/jset-grouped.tsv"
: >"$OUT_ROOT/performance/two-point.tsv"
for run in 1 2 3 4 5 6 7; do
  /usr/bin/time -f "$run\t%e\t%M" -o "$OUT_ROOT/performance/jset-unscoped.tsv" -a \
    "$AIE_BIN" query "$d0_archive" jset \
      --include "$include_locus" --exclude "$exclude_locus" --json \
      >/dev/null 2>/dev/null
  /usr/bin/time -f "$run\t%e\t%M" -o "$OUT_ROOT/performance/jset-grouped.tsv" -a \
    "$AIE_BIN" query "$d0_archive" jset \
      --include "$include_locus" --exclude "$exclude_locus" \
      --groups "$d0_group_path" --json >/dev/null 2>/dev/null
  AIE_GATE_BIN="$AIE_BIN" AIE_GATE_ARCHIVE="$d0_archive" \
    AIE_GATE_INCLUDE="$include_locus" AIE_GATE_EXCLUDE="$exclude_locus" \
    /usr/bin/time -f "$run\t%e\t%M" -o "$OUT_ROOT/performance/two-point.tsv" -a \
    /usr/bin/bash -c '
      "$AIE_GATE_BIN" query "$AIE_GATE_ARCHIVE" junction "$AIE_GATE_INCLUDE" --json >/dev/null 2>/dev/null
      "$AIE_GATE_BIN" query "$AIE_GATE_ARCHIVE" junction "$AIE_GATE_EXCLUDE" --json >/dev/null 2>/dev/null
    '
done

# Locked unscoped outputs ensure adding the scope substrate did not change existing interfaces.
batch_plan=$script_dir/../experiments/plans/batched-query-d0.tsv
"$AIE_BIN" query "$d0_archive" batch --plan "$batch_plan" \
  >"$OUT_ROOT/compatibility/batch.json" 2>"$OUT_ROOT/compatibility/batch.stderr"
"$AIE_BIN" query "$d0_archive" junction "$include_locus" --top 0 --json \
  >"$OUT_ROOT/compatibility/junction.json" 2>"$OUT_ROOT/compatibility/junction.stderr"
"$AIE_BIN" query "$d0_archive" junctions chr11:34050000-34075000 --with-cells --json \
  >"$OUT_ROOT/compatibility/junctions.json" 2>"$OUT_ROOT/compatibility/junctions.stderr"
"$AIE_BIN" query "$d0_archive" region chr11:34050000-34075000 --top 0 \
  >"$OUT_ROOT/compatibility/region.txt" 2>"$OUT_ROOT/compatibility/region.stderr"

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-31\n'
  printf 'gravlax_commit\t%s\n' "$commit"
  printf 'threads\t%s\n' "$threads"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'd0_archive_sha256\t%s\n' "$(sha256sum "$d0_archive" | cut -d' ' -f1)"
  printf 'd1_archive_sha256\t%s\n' "$(sha256sum "$d1_archive" | cut -d' ' -f1)"
  printf 'd0_groups_sha256\t%s\n' "$(sha256sum "$d0_group_path" | cut -d' ' -f1)"
  printf 'd1_groups_sha256\t%s\n' "$(sha256sum "$d1_group_path" | cut -d' ' -f1)"
  printf 'plan_sha256\t%s\n' "$(sha256sum "$plan" | cut -d' ' -f1)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'cpu\t%s\n' "$(lscpu | awk -F: '$1 == "Model name" {sub(/^[[:space:]]+/, "", $2); print $2}')"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

python3 "$script_dir/176_validate_cell_group_junction_set_query.py" \
  "$OUT_ROOT" --d0-groups "$d0_group_path" --d1-groups "$d1_group_path" \
  --commit "$commit" --out "$OUT_ROOT/summary.json"
printf 'cell/group junction-set gates passed: %s\n' "$OUT_ROOT/summary.json"
