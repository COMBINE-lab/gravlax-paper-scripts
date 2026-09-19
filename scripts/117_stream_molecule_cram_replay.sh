#!/usr/bin/env bash
# Quantify directly from a post-correction molecule CRAM without materializing a BAM.
#
# samtools decodes CRAM to BAM through a named pipe while Gravlax reconstructs the molecule
# abstraction and replays a GTF. Both processes share one CPU affinity set, so their combined
# execution cannot exceed the declared hardware-thread budget.
#
# Required environment:
#   PROJ_ROOT, AIE_BIN, MOLECULE_CRAM, DATASET, BARCODES, REFERENCE_MATRIX, OUT_ROOT
# Optional environment:
#   GTF (default gencode v49), GENOME (default GRCh38 primary), SOLO_STRAND (forward),
#   SAMTOOLS, THREADS (24),
#   CPUSET (0-23), WARM_REPS (5), GRAVLAX_REPO
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN}"
: "${MOLECULE_CRAM:?set MOLECULE_CRAM}"
: "${DATASET:?set DATASET}"
: "${BARCODES:?set BARCODES}"
: "${REFERENCE_MATRIX:?set REFERENCE_MATRIX}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
genome=${GENOME:-$PROJ_ROOT/ref/GRCh38.primary_assembly.genome.fa}
samtools_bin=${SAMTOOLS:-samtools}
threads=${THREADS:-24}
cpuset=${CPUSET:-0-23}
warm_reps=${WARM_REPS:-5}
gravlax_repo=${GRAVLAX_REPO:-$PROJ_ROOT/work/gravlax}
solo_strand=${SOLO_STRAND:-forward}
run_date=${RUN_DATE:-$(date +%F)}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be positive" >&2; exit 2; }
[[ "$warm_reps" =~ ^[1-9][0-9]*$ ]] || { echo "WARM_REPS must be positive" >&2; exit 2; }
[[ -x "$AIE_BIN" ]] || { echo "AIE_BIN is not executable: $AIE_BIN" >&2; exit 2; }
command -v "$samtools_bin" >/dev/null || { echo "samtools not found: $samtools_bin" >&2; exit 2; }
command -v taskset >/dev/null || { echo "taskset not found" >&2; exit 2; }
git -C "$gravlax_repo" rev-parse --verify HEAD >/dev/null 2>&1 || {
  echo "GRAVLAX_REPO is not a Git checkout: $gravlax_repo" >&2
  exit 2
}
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }

for input in "$MOLECULE_CRAM" "$genome" "$genome.fai" "$gtf" "$BARCODES" \
  "$REFERENCE_MATRIX/matrix.mtx" "$REFERENCE_MATRIX/features.tsv" \
  "$REFERENCE_MATRIX/barcodes.tsv"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done

mkdir -p "$OUT_ROOT"
"$samtools_bin" quickcheck -v "$MOLECULE_CRAM"
export RAYON_NUM_THREADS=$threads

run_one() {
  local label=$1
  local run_dir=$OUT_ROOT/$label
  local fifo=$run_dir/molecules.bam.fifo
  mkdir -p "$run_dir/replay"
  mkfifo "$fifo"

  local start_ns end_ns decoder_pid decoder_status replay_status
  start_ns=$(date +%s%N)
  /usr/bin/time -v -o "$run_dir/decoder.time.txt" \
    taskset -c "$cpuset" "$samtools_bin" view -@ "$threads" -T "$genome" -b \
      -o "$fifo" "$MOLECULE_CRAM" \
      >"$run_dir/decoder.stdout.txt" 2>"$run_dir/decoder.stderr.txt" &
  decoder_pid=$!

  set +e
  /usr/bin/time -v -o "$run_dir/replay.time.txt" \
    taskset -c "$cpuset" "$AIE_BIN" replay-rows "$fifo" --from-molecule-bam \
      --gtf "$gtf" --barcodes "$BARCODES" --out-dir "$run_dir/replay" \
      --solo-strand "$solo_strand" \
      >"$run_dir/replay.stdout.txt" 2>"$run_dir/replay.stderr.txt"
  replay_status=$?
  wait "$decoder_pid"
  decoder_status=$?
  set -e
  end_ns=$(date +%s%N)
  rm -f "$fifo"

  python3 -c 'import sys; print((int(sys.argv[2])-int(sys.argv[1]))/1e9)' \
    "$start_ns" "$end_ns" >"$run_dir/pipeline-wall-seconds.txt"
  printf 'decoder_exit\t%s\nreplay_exit\t%s\n' "$decoder_status" "$replay_status" \
    >"$run_dir/exit-status.tsv"
  if [[ "$decoder_status" != 0 || "$replay_status" != 0 ]]; then
    echo "$label failed: decoder=$decoder_status replay=$replay_status" >&2
    exit 1
  fi

  for artifact in matrix.mtx barcodes.tsv features.tsv; do
    cmp "$run_dir/replay/$artifact" "$REFERENCE_MATRIX/$artifact"
  done
}

{
  printf 'key\tvalue\n'
  printf 'date\t%s\n' "$run_date"
  printf 'dataset\t%s\n' "$DATASET"
  printf 'threads_per_process\t%s\n' "$threads"
  printf 'shared_cpuset\t%s\n' "$cpuset"
  printf 'combined_hardware_thread_budget\t%s\n' "$threads"
  printf 'warm_replicates\t%s\n' "$warm_reps"
  printf 'solo_strand\t%s\n' "$solo_strand"
  printf 'gravlax_commit\t%s\n' "$(git -C "$gravlax_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'samtools\t%s\n' "$("$samtools_bin" --version | head -1)"
  printf 'cram_sha256\t%s\n' "$(sha256sum "$MOLECULE_CRAM" | cut -d' ' -f1)"
  printf 'reference_sha256\t%s\n' "$(sha256sum "$genome" | cut -d' ' -f1)"
  printf 'transport\tnamed pipe; no materialized BAM\n'
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

python3 "$script_dir/115_fadvise_dontneed.py" --log "$OUT_ROOT/cold-cache.tsv" \
  "$MOLECULE_CRAM" "$genome" "$genome.fai" "$gtf" "$BARCODES"
run_one cold-r1

for rep in $(seq 1 "$warm_reps"); do
  run_one "warm-r$rep"
done

printf 'streamed molecule-CRAM replay passed: %s\n' "$OUT_ROOT"
