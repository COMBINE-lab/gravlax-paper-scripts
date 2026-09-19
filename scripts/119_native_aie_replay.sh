#!/usr/bin/env bash
# Dataset-generic native .aie replay timing and exactness gate.
#
# One best-effort advisory-cold observation is followed by warm replicates. The command is pinned
# to the same CPU set used by scripts/117_stream_molecule_cram_replay.sh so the native and streamed
# CRAM replay paths have the same hardware-thread ceiling.
#
# Required environment:
#   PROJ_ROOT, AIE_BIN, ARCHIVE, DATASET, BARCODES, REFERENCE_MATRIX, OUT_ROOT
# Optional environment:
#   GTF (default gencode v49), THREADS (24), CPUSET (0-23), WARM_REPS (5), GRAVLAX_REPO
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN}"
: "${ARCHIVE:?set ARCHIVE}"
: "${DATASET:?set DATASET}"
: "${BARCODES:?set BARCODES}"
: "${REFERENCE_MATRIX:?set REFERENCE_MATRIX}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
threads=${THREADS:-24}
cpuset=${CPUSET:-0-23}
warm_reps=${WARM_REPS:-5}
gravlax_repo=${GRAVLAX_REPO:-$PROJ_ROOT/work/gravlax}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be positive" >&2; exit 2; }
[[ "$warm_reps" =~ ^[1-9][0-9]*$ ]] || { echo "WARM_REPS must be positive" >&2; exit 2; }
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
  local label=$1
  local run_dir=$OUT_ROOT/$label
  local start_ns end_ns
  mkdir -p "$run_dir/replay"

  start_ns=$(date +%s%N)
  /usr/bin/time -v -o "$run_dir/replay.time.txt" \
    taskset -c "$cpuset" "$AIE_BIN" replay-rows "$ARCHIVE" \
      --gtf "$gtf" --barcodes "$BARCODES" --out-dir "$run_dir/replay" \
      >"$run_dir/replay.stdout.txt" 2>"$run_dir/replay.stderr.txt"
  end_ns=$(date +%s%N)

  python3 -c 'import sys; print((int(sys.argv[2])-int(sys.argv[1]))/1e9)' \
    "$start_ns" "$end_ns" >"$run_dir/pipeline-wall-seconds.txt"
  for artifact in matrix.mtx barcodes.tsv features.tsv; do
    cmp "$run_dir/replay/$artifact" "$REFERENCE_MATRIX/$artifact"
  done
}

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-30\n'
  printf 'dataset\t%s\n' "$DATASET"
  printf 'threads\t%s\n' "$threads"
  printf 'cpuset\t%s\n' "$cpuset"
  printf 'hardware_thread_budget\t%s\n' "$threads"
  printf 'warm_replicates\t%s\n' "$warm_reps"
  printf 'gravlax_commit\t%s\n' "$(git -C "$gravlax_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'archive_sha256\t%s\n' "$(sha256sum "$ARCHIVE" | cut -d' ' -f1)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
  printf 'cold_definition\tPOSIX_FADV_DONTNEED on archive, GTF, and barcodes; advisory\n'
} >"$OUT_ROOT/protocol.tsv"

python3 "$script_dir/115_fadvise_dontneed.py" --log "$OUT_ROOT/cold-cache.tsv" \
  "$ARCHIVE" "$gtf" "$BARCODES"
run_one cold-r1

for rep in $(seq 1 "$warm_reps"); do
  run_one "warm-r$rep"
done

printf 'native .aie replay passed: %s\n' "$OUT_ROOT"
