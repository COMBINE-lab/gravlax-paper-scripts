#!/usr/bin/env bash
# Dataset-generic randomized Gene-only STARsolo versus native .aie replay benchmark.
#
# One advisory-cold observation per arm is followed by five warm paired blocks in an alternating
# two-arm crossover. Every run is pinned to one CPU set and compared with its own frozen reference.
# STARsolo emits counts only; alignment output is forbidden.
#
# Required environment:
#   PROJ_ROOT, AIE_BIN, STAR_BIN, DATASET, ARCHIVE, BARCODES, AIE_REFERENCE_MATRIX,
#   R1_CSV, R2_CSV, STAR_ORACLE_SOLO, OUT_ROOT
# Optional environment:
#   GTF, GENOME_DIR, WHITELIST, THREADS (24), CPUSET (0-23), WARM_BLOCKS (5), GRAVLAX_REPO
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN}"
: "${STAR_BIN:?set STAR_BIN}"
: "${DATASET:?set DATASET}"
: "${ARCHIVE:?set ARCHIVE}"
: "${BARCODES:?set BARCODES}"
: "${AIE_REFERENCE_MATRIX:?set AIE_REFERENCE_MATRIX}"
: "${R1_CSV:?set R1_CSV to comma-separated FASTQs}"
: "${R2_CSV:?set R2_CSV to comma-separated FASTQs}"
: "${STAR_ORACLE_SOLO:?set STAR_ORACLE_SOLO}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
genome_dir=${GENOME_DIR:-$PROJ_ROOT/ref/star-49}
whitelist=${WHITELIST:-$PROJ_ROOT/data/3M-february-2018.txt}
threads=${THREADS:-24}
cpuset=${CPUSET:-0-23}
warm_blocks=${WARM_BLOCKS:-5}
gravlax_repo=${GRAVLAX_REPO:-$PROJ_ROOT/work/gravlax}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be positive" >&2; exit 2; }
[[ "$warm_blocks" =~ ^[1-9][0-9]*$ ]] || { echo "WARM_BLOCKS must be positive" >&2; exit 2; }
[[ -x "$AIE_BIN" ]] || { echo "AIE_BIN is not executable: $AIE_BIN" >&2; exit 2; }
[[ -x "$STAR_BIN" ]] || { echo "STAR_BIN is not executable: $STAR_BIN" >&2; exit 2; }
command -v taskset >/dev/null || { echo "taskset not found" >&2; exit 2; }
git -C "$gravlax_repo" rev-parse --verify HEAD >/dev/null 2>&1 || {
  echo "GRAVLAX_REPO is not a Git checkout: $gravlax_repo" >&2
  exit 2
}
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }

IFS=',' read -r -a r1_files <<<"$R1_CSV"
IFS=',' read -r -a r2_files <<<"$R2_CSV"
[[ ${#r1_files[@]} -eq ${#r2_files[@]} ]] || {
  echo "R1_CSV and R2_CSV must have the same number of files" >&2
  exit 2
}
for input in "$ARCHIVE" "$BARCODES" "$gtf" "$whitelist" \
  "$AIE_REFERENCE_MATRIX/matrix.mtx" "$AIE_REFERENCE_MATRIX/features.tsv" \
  "$AIE_REFERENCE_MATRIX/barcodes.tsv" "$genome_dir/Genome" "$genome_dir/SA" \
  "$genome_dir/SAindex" "${r1_files[@]}" "${r2_files[@]}" \
  "$STAR_ORACLE_SOLO/Gene/raw/matrix.mtx" "$STAR_ORACLE_SOLO/Gene/raw/features.tsv" \
  "$STAR_ORACLE_SOLO/Gene/raw/barcodes.tsv"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done

mkdir -p "$OUT_ROOT"
export RAYON_NUM_THREADS=$threads

run_aie() {
  local label=$1
  local run_dir=$OUT_ROOT/$label
  mkdir -p "$run_dir/replay"
  /usr/bin/time -v -o "$run_dir/time.txt" \
    taskset -c "$cpuset" "$AIE_BIN" replay-rows "$ARCHIVE" \
      --gtf "$gtf" --barcodes "$BARCODES" --out-dir "$run_dir/replay" \
      >"$run_dir/command.stdout.txt" 2>"$run_dir/command.stderr.txt"
  for artifact in matrix.mtx barcodes.tsv features.tsv; do
    cmp "$run_dir/replay/$artifact" "$AIE_REFERENCE_MATRIX/$artifact"
  done
  printf 'completed exact %s\n' "$label"
}

run_star() {
  local label=$1
  local run_dir=$OUT_ROOT/$label
  mkdir -p "$run_dir"
  /usr/bin/time -v -o "$run_dir/time.txt" \
    taskset -c "$cpuset" "$STAR_BIN" \
      --runThreadN "$threads" \
      --genomeDir "$genome_dir" \
      --readFilesIn "$R2_CSV" "$R1_CSV" \
      --readFilesCommand zcat \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist "$whitelist" \
      --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
      --soloBarcodeReadLength 0 \
      --soloFeatures Gene \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
      --soloUMIdedup 1MM_CR \
      --soloUMIfiltering MultiGeneUMI_CR \
      --soloCellFilter EmptyDrops_CR \
      --outSAMtype None \
      --outSAMunmapped Within \
      --outFileNamePrefix "$run_dir/" \
      >"$run_dir/command.stdout.txt" 2>"$run_dir/command.stderr.txt"

  while IFS= read -r artifact; do
    cmp "$run_dir/Solo.out/$artifact" "$STAR_ORACLE_SOLO/$artifact"
  done < <(cd "$run_dir/Solo.out" && find . -type f -printf '%P\n' | sort)
  if find "$run_dir" -maxdepth 1 -type f \( -name '*.bam' -o -name '*.sam' \) \
      -print -quit | grep -q .; then
    echo "$label unexpectedly emitted an alignment file" >&2
    exit 1
  fi
  printf 'completed exact %s\n' "$label"
}

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-30\n'
  printf 'dataset\t%s\n' "$DATASET"
  printf 'threads\t%s\n' "$threads"
  printf 'cpuset\t%s\n' "$cpuset"
  printf 'warm_blocks\t%s\n' "$warm_blocks"
  printf 'schedule_seed\t20260830\n'
  printf 'schedule_design\tbalanced alternating two-arm crossover; STAR first in odd blocks\n'
  printf 'gravlax_commit\t%s\n' "$(git -C "$gravlax_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'archive_sha256\t%s\n' "$(sha256sum "$ARCHIVE" | cut -d' ' -f1)"
  printf 'star_version\t%s\n' "$("$STAR_BIN" --version)"
  printf 'star_features\tGene\n'
  printf 'star_alignment_output\tNone\n'
  printf 'r1_csv\t%s\n' "$R1_CSV"
  printf 'r2_csv\t%s\n' "$R2_CSV"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
  printf 'cold_definition\tPOSIX_FADV_DONTNEED on named arm inputs; advisory and non-destructive\n'
} >"$OUT_ROOT/protocol.tsv"
printf 'phase\tblock\tposition\tarm\tlabel\n' >"$OUT_ROOT/schedule.tsv"

python3 "$script_dir/115_fadvise_dontneed.py" --log "$OUT_ROOT/cold-aie-cache.tsv" \
  "$ARCHIVE" "$gtf" "$BARCODES"
printf 'cold\t0\t1\taie\tcold-aie\n' >>"$OUT_ROOT/schedule.tsv"
run_aie cold-aie

python3 "$script_dir/115_fadvise_dontneed.py" --log "$OUT_ROOT/cold-star-cache.tsv" \
  "${r1_files[@]}" "${r2_files[@]}" "$genome_dir" "$whitelist"
printf 'cold\t0\t2\tstar\tcold-star\n' >>"$OUT_ROOT/schedule.tsv"
run_star cold-star

for block in $(seq 1 "$warm_blocks"); do
  if (( block % 2 == 1 )); then
    arms=(star aie)
  else
    arms=(aie star)
  fi
  position=0
  for arm in "${arms[@]}"; do
    position=$((position + 1))
    label=warm-b${block}-p${position}-${arm}
    printf 'warm\t%s\t%s\t%s\t%s\n' "$block" "$position" "$arm" "$label" \
      >>"$OUT_ROOT/schedule.tsv"
    if [[ "$arm" == aie ]]; then
      run_aie "$label"
    else
      run_star "$label"
    fi
  done
done

printf 'scale randomized counts comparison passed: %s\n' "$OUT_ROOT"
