#!/usr/bin/env bash
# Prepare one untouched dataset without observing any EM score.
set -euo pipefail

dataset=${1:-}
case "$dataset" in D4|D3) ;; *) echo "usage: $0 D4|D3" >&2; exit 2 ;; esac
: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN}"
: "${STAR_BIN:?set STAR_BIN}"
: "${PYTHON_BIN:?set PYTHON_BIN to the frozen Scanpy environment}"
: "${OUT_ROOT:?set OUT_ROOT to a new dataset preparation directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
manifest=${MANIFEST:-$scripts_repo/experiments/manifests/depth-hybrid-em-confirmation.yaml}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
genome_dir=${GENOME_DIR:-$PROJ_ROOT/ref/star-49}
whitelist=${WHITELIST:-$PROJ_ROOT/data/3M-february-2018.txt}
threads=${THREADS:-24}
resume=${RESUME:-0}
archive=$PROJ_ROOT/runs/archive/${dataset,,}.aie
if [[ "$dataset" == D4 ]]; then
  fastq_dir=$PROJ_ROOT/data/d4/500_PBMC_3p_LT_Chromium_X_fastqs
  prefix=500_PBMC_3p_LT_Chromium_X_S4
else
  fastq_dir=$PROJ_ROOT/data/d3/pbmc_10k_v3_fastqs
  prefix=pbmc_10k_v3_S1
fi
r1_csv=$(find "$fastq_dir" -maxdepth 1 -type f -name "${prefix}_*_R1_*.fastq.gz" -print | sort | paste -sd,)
r2_csv=$(find "$fastq_dir" -maxdepth 1 -type f -name "${prefix}_*_R2_*.fastq.gz" -print | sort | paste -sd,)

for input in "$AIE_BIN" "$STAR_BIN" "$PYTHON_BIN" "$manifest" "$gtf" "$archive" \
  "$genome_dir/Genome" "$genome_dir/SA" "$genome_dir/SAindex" "$whitelist"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done
[[ -n "$r1_csv" && -n "$r2_csv" ]] || { echo "FASTQ inputs not found" >&2; exit 2; }
if [[ -e "$OUT_ROOT" && "$resume" != 1 ]]; then
  echo "refusing to overwrite OUT_ROOT: $OUT_ROOT (set RESUME=1 for a checked resume)" >&2
  exit 2
fi
mkdir -p "$OUT_ROOT/starsolo" "$OUT_ROOT/candidates" "$OUT_ROOT/logs"
export RAYON_NUM_THREADS=$threads
export MPLCONFIGDIR=$OUT_ROOT/matplotlib
export NUMBA_CACHE_DIR=$OUT_ROOT/numba-cache

called=$OUT_ROOT/starsolo/Solo.out/Gene/filtered/barcodes.tsv
if [[ ! -s "$called" ]]; then
  /usr/bin/time -v -o "$OUT_ROOT/starsolo/time.txt" "$STAR_BIN" \
    --runThreadN "$threads" --genomeDir "$genome_dir" \
    --readFilesIn "$r2_csv" "$r1_csv" --readFilesCommand zcat \
    --soloType CB_UMI_Simple --soloCBwhitelist "$whitelist" \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 --soloFeatures Gene \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIdedup 1MM_CR --soloUMIfiltering MultiGeneUMI_CR \
    --soloCellFilter EmptyDrops_CR --outSAMtype None --outSAMunmapped Within \
    --outFileNamePrefix "$OUT_ROOT/starsolo/" \
    >"$OUT_ROOT/starsolo/command.stdout.txt" 2>"$OUT_ROOT/starsolo/command.stderr.txt"
fi
[[ -s "$called" ]] || { echo "fresh STARsolo emitted no called-cell barcodes" >&2; exit 1; }
for seed in 7 17 29; do
  candidate=$OUT_ROOT/candidates/$dataset.seed-$seed.genes.txt
  if [[ ! -s "$candidate" ]]; then
    /usr/bin/time -v -o "$OUT_ROOT/logs/candidates.seed-$seed.time.txt" \
      "$AIE_BIN" em "$archive" --gtf "$gtf" --mask 0.2 --seed "$seed" \
        --candidate-genes-out "$candidate" --candidate-genes-only \
        >"$OUT_ROOT/logs/candidates.seed-$seed.stdout.txt" \
        2>"$OUT_ROOT/logs/candidates.seed-$seed.stderr.txt"
  fi
done

if [[ ! -f "$OUT_ROOT/replay.complete" ]]; then
  [[ ! -e "$OUT_ROOT/replay" ]] || {
    echo "incomplete replay directory exists; preserve or move it before resume" >&2
    exit 2
  }
  /usr/bin/time -v -o "$OUT_ROOT/replay.time.txt" \
    "$AIE_BIN" replay-rows "$archive" --gtf "$gtf" --barcodes "$whitelist" \
      --out-dir "$OUT_ROOT/replay" \
      >"$OUT_ROOT/logs/replay.stdout.txt" 2>"$OUT_ROOT/logs/replay.stderr.txt"
  touch "$OUT_ROOT/replay.complete"
fi

if [[ ! -s "$OUT_ROOT/groups/groups.json" ]]; then
  [[ ! -e "$OUT_ROOT/groups" ]] || {
    echo "incomplete group directory exists; preserve or move it before resume" >&2
    exit 2
  }
  "$PYTHON_BIN" "$scripts_repo/scripts/161_build_em_confirmation_groups.py" \
    --dataset "$dataset" --manifest "$manifest" --replay-dir "$OUT_ROOT/replay" \
    --called-barcodes "$called" --candidate-dir "$OUT_ROOT/candidates" \
    --out-dir "$OUT_ROOT/groups"
fi

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-31\n'
  printf 'dataset\t%s\n' "$dataset"
  printf 'threads\t%s\n' "$threads"
  printf 'gravlax_commit\t%s\n' "$(git -C "$PROJ_ROOT/work/gravlax" rev-parse HEAD)"
  printf 'scripts_commit\t%s\n' "$(git -C "$scripts_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$manifest" | cut -d' ' -f1)"
  printf 'archive_sha256\t%s\n' "$(sha256sum "$archive" | cut -d' ' -f1)"
  printf 'gtf_sha256\t%s\n' "$(sha256sum "$gtf" | cut -d' ' -f1)"
  printf 'whitelist_sha256\t%s\n' "$(sha256sum "$whitelist" | cut -d' ' -f1)"
  printf 'STAR\t%s\n' "$($STAR_BIN --version)"
  printf 'python\t%s\n' "$PYTHON_BIN"
  printf 'python_version\t%s\n' "$($PYTHON_BIN --version 2>&1)"
  printf 'r1_csv\t%s\n' "$r1_csv"
  printf 'r2_csv\t%s\n' "$r2_csv"
  printf 'called_cells\t%s\n' "$(wc -l < "$called")"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

printf '%s confirmation inputs prepared: %s\n' "$dataset" "$OUT_ROOT"
