#!/usr/bin/env bash
# Execute the prospectively locked FNBP1 confirmation on the public 10x normal-brain dataset.
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 ROOT OUT" >&2
  exit 2
fi

ROOT=$(realpath "$1")
OUT=$(realpath -m "$2")
mkdir -p "$OUT"

AIE=${AIE:-$ROOT/env/cargo-target-features/release/aie}
SAMTOOLS=${SAMTOOLS:-$ROOT/env/sc/bin/samtools}
STAR=${STAR:-$ROOT/env/sc/bin/STAR}
PYTHON=${PYTHON:-python3}
THREADS=${THREADS:-8}
BASE=https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.2.0/10k_Human_Brain_MO_gemx

fetch() {
  local name=$1
  local expected=$2
  local url=$3
  if [[ ! -f "$OUT/$name" ]] || ! echo "$expected  $OUT/$name" | sha256sum -c - >/dev/null 2>&1; then
    aria2c -c -x 8 -s 8 --file-allocation=none --dir="$OUT" --out="$name" "$url"
  fi
  echo "$expected  $OUT/$name" | sha256sum -c -
}

fetch filtered_feature_bc_matrix.h5 \
  5dc9d83c2ce7b0f6fbc5f6b9ee75b2fa3cdb7c470b027b7c6284a8381fe2f99a \
  "$BASE/10k_Human_Brain_MO_gemx_filtered_feature_bc_matrix.h5"
fetch gex_possorted_bam.bam \
  e49697023eb7f5cca2bac7330a6c032d0ce40c0f3e9dcc703b11e4e72fb1dbbb \
  "$BASE/10k_Human_Brain_MO_gemx_gex_possorted_bam.bam"
fetch gex_possorted_bam.bam.bai \
  61fe1a7d796fce3ed1388743f80e7ed0107b5b4055a986d1e32e86311793f046 \
  "$BASE/10k_Human_Brain_MO_gemx_gex_possorted_bam.bam.bai"
"$SAMTOOLS" quickcheck -v "$OUT/gex_possorted_bam.bam"

PYTHONDONTWRITEBYTECODE=1 "$PYTHON" \
  "$ROOT/gravlax-paper-scripts/scripts/182_prepare_fnbp1_confirmation.py" \
  --h5 "$OUT/filtered_feature_bc_matrix.h5" \
  --bam "$OUT/gex_possorted_bam.bam" \
  --samtools "$SAMTOOLS" --out "$OUT" --threads "$THREADS"

AIC="$OUT/gencode.v49.aic"
if [[ ! -s "$AIC" ]]; then
  "$AIE" compile-annotation "$ROOT/annotations/gencode.v49.annotation.gtf" --out "$AIC"
fi

if [[ ! -s "$OUT/fnbp1.called.aie" ]]; then
  "$AIE" ingest-archive "$OUT/fnbp1.called.bam" \
    --whitelist "$OUT/called-whitelist.txt" \
    --genome "$ROOT/data/GRCh38.primary_assembly.genome.fa.gz" \
    --out "$OUT/fnbp1.called.aie"
fi
"$AIE" query "$OUT/fnbp1.called.aie" events chr9:129907500-129926000 \
  --min-support 2 --min-informative 1 --max-events 100000 --gtf "$AIC" --json \
  > "$OUT/events.json"

STAR_OUT="$OUT/star-targeted-t8"
if [[ ! -s "$STAR_OUT/Aligned.sortedByCoord.out.bam" ]]; then
  mkdir -p "$STAR_OUT"
  "$STAR" --runThreadN "$THREADS" --genomeDir "$ROOT/ref/star-base" --twopassMode Basic \
    --readFilesIn "$OUT/fnbp1.R2.fastq" "$OUT/fnbp1.R1.fastq" \
    --soloType CB_UMI_Simple --soloCBwhitelist "$OUT/called-whitelist.txt" \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 --soloFeatures SJ \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --outSAMattributes NH HI AS nM CR UR CY UY --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 32000000000 --outSAMunmapped Within --outFileNamePrefix "$STAR_OUT/"
fi
if [[ ! -s "$OUT/fnbp1.star-targeted.aie" ]]; then
  "$AIE" ingest-archive "$STAR_OUT/Aligned.sortedByCoord.out.bam" \
    --whitelist "$OUT/called-whitelist.txt" \
    --genome "$ROOT/data/GRCh38.primary_assembly.genome.fa.gz" \
    --out "$OUT/fnbp1.star-targeted.aie"
fi
"$AIE" query "$OUT/fnbp1.star-targeted.aie" events chr9:129907500-129926000 \
  --min-support 2 --min-informative 1 --max-events 100000 --gtf "$AIC" --json \
  > "$OUT/star-events.json"

PYTHONDONTWRITEBYTECODE=1 "$PYTHON" \
  "$ROOT/gravlax-paper-scripts/scripts/183_summarize_fnbp1_confirmation.py" \
  --run "$OUT" --samtools "$SAMTOOLS" --out "$OUT/summary.json"
