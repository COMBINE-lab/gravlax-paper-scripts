#!/usr/bin/env bash
# Randomized same-binary comparison of bounded streaming and eager .aie Gene replay.
#
# Required environment:
#   PROJ_ROOT, AIE_BIN, DATASET, ARCHIVE, BARCODES, REFERENCE_MATRIX, OUT_ROOT
# Optional environment:
#   GTF, THREADS (24), CPUSET (0-23), WARM_BLOCKS (5), GRAVLAX_REPO
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN}"
: "${DATASET:?set DATASET}"
: "${ARCHIVE:?set ARCHIVE}"
: "${BARCODES:?set BARCODES}"
: "${REFERENCE_MATRIX:?set REFERENCE_MATRIX}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
threads=${THREADS:-24}
cpuset=${CPUSET:-0-23}
warm_blocks=${WARM_BLOCKS:-5}
gravlax_repo=${GRAVLAX_REPO:-$PROJ_ROOT/work/gravlax}

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be positive" >&2; exit 2; }
[[ "$warm_blocks" =~ ^[1-9][0-9]*$ ]] || { echo "WARM_BLOCKS must be positive" >&2; exit 2; }
[[ -x "$AIE_BIN" ]] || { echo "AIE_BIN is not executable: $AIE_BIN" >&2; exit 2; }
command -v taskset >/dev/null || { echo "taskset not found" >&2; exit 2; }
git -C "$gravlax_repo" rev-parse --verify HEAD >/dev/null 2>&1 || {
  echo "GRAVLAX_REPO is not a Git checkout: $gravlax_repo" >&2
  exit 2
}
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }

for input in "$ARCHIVE" "$gtf" "$BARCODES" \
  "$REFERENCE_MATRIX/matrix.mtx" "$REFERENCE_MATRIX/features.tsv" \
  "$REFERENCE_MATRIX/barcodes.tsv"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done

mkdir -p "$OUT_ROOT"
export RAYON_NUM_THREADS=$threads

run_one() {
  local block=$1
  local position=$2
  local arm=$3
  local label="warm-b${block}-p${position}-${arm}"
  local run_dir=$OUT_ROOT/$label
  local -a mode=()
  if [[ "$arm" == eager ]]; then
    mode=(--eager)
  fi
  mkdir -p "$run_dir/replay"
  /usr/bin/time -v -o "$run_dir/time.txt" \
    taskset -c "$cpuset" "$AIE_BIN" replay-rows "$ARCHIVE" "${mode[@]}" \
      --gtf "$gtf" --barcodes "$BARCODES" --out-dir "$run_dir/replay" \
      >"$run_dir/command.stdout.txt" 2>"$run_dir/command.stderr.txt"
  for artifact in matrix.mtx features.tsv barcodes.tsv; do
    cmp "$run_dir/replay/$artifact" "$REFERENCE_MATRIX/$artifact"
  done
  printf 'warm\t%s\t%s\t%s\t%s\n' "$block" "$position" "$arm" "$label" \
    >>"$OUT_ROOT/schedule.tsv"
}

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-30\n'
  printf 'dataset\t%s\n' "$DATASET"
  printf 'threads\t%s\n' "$threads"
  printf 'cpuset\t%s\n' "$cpuset"
  printf 'warm_blocks\t%s\n' "$warm_blocks"
  printf 'schedule_seed\t20260830\n'
  printf 'schedule_design\tbalanced alternating two-arm crossover; eager first in odd blocks\n'
  printf 'streaming_mode\tdefault bounded archive chunk decode and global cell-shard reducer\n'
  printf 'eager_mode\t--eager full archive materialization; same binary and reducer semantics\n'
  printf 'gravlax_commit\t%s\n' "$(git -C "$gravlax_repo" rev-parse HEAD)"
  printf 'dirty_worktree\t%s\n' "$(git -C "$gravlax_repo" status --porcelain | wc -l)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'archive_sha256\t%s\n' "$(sha256sum "$ARCHIVE" | cut -d' ' -f1)"
  printf 'gtf_sha256\t%s\n' "$(sha256sum "$gtf" | cut -d' ' -f1)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"
printf 'phase\tblock\tposition\tarm\tlabel\n' >"$OUT_ROOT/schedule.tsv"

for block in $(seq 1 "$warm_blocks"); do
  if (( block % 2 == 1 )); then
    arms=(eager streaming)
  else
    arms=(streaming eager)
  fi
  position=0
  for arm in "${arms[@]}"; do
    position=$((position + 1))
    run_one "$block" "$position" "$arm"
  done
done

printf 'streaming replay gate passed exactness: %s\n' "$OUT_ROOT"
