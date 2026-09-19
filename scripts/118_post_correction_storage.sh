#!/usr/bin/env bash
# Dataset-generic post-correction molecule BAM/CRAM storage and exactness gate.
#
# Required environment:
#   PROJ_ROOT, AIE_BIN, ARCHIVE, DATASET, BARCODES, REFERENCE_MATRIX, OUT_DIR
# Optional: GENOME, GTF, SOLO_STRAND (forward), SAMTOOLS, THREADS, GRAVLAX_REPO
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN}"
: "${ARCHIVE:?set ARCHIVE}"
: "${DATASET:?set DATASET}"
: "${BARCODES:?set BARCODES}"
: "${REFERENCE_MATRIX:?set REFERENCE_MATRIX}"
: "${OUT_DIR:?set OUT_DIR to a new directory}"

genome=${GENOME:-$PROJ_ROOT/ref/GRCh38.primary_assembly.genome.fa}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
samtools_bin=${SAMTOOLS:-samtools}
threads=${THREADS:-24}
gravlax_repo=${GRAVLAX_REPO:-$PROJ_ROOT/work/gravlax}
solo_strand=${SOLO_STRAND:-forward}
run_date=${RUN_DATE:-$(date +%F)}

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be positive" >&2; exit 2; }
[[ -x "$AIE_BIN" ]] || { echo "AIE_BIN is not executable: $AIE_BIN" >&2; exit 2; }
command -v "$samtools_bin" >/dev/null || { echo "samtools not found: $samtools_bin" >&2; exit 2; }
git -C "$gravlax_repo" rev-parse --verify HEAD >/dev/null 2>&1 || {
  echo "GRAVLAX_REPO is not a Git checkout: $gravlax_repo" >&2
  exit 2
}
[[ ! -e "$OUT_DIR" ]] || { echo "refusing to overwrite OUT_DIR: $OUT_DIR" >&2; exit 2; }

for input in "$ARCHIVE" "$genome" "$genome.fai" "$gtf" "$BARCODES" \
  "$REFERENCE_MATRIX/matrix.mtx" "$REFERENCE_MATRIX/features.tsv" \
  "$REFERENCE_MATRIX/barcodes.tsv"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done

mkdir -p "$OUT_DIR/logs" "$OUT_DIR/replay-bam"
molecule_bam=$OUT_DIR/post-correction.molecules.bam
molecule_cram=$OUT_DIR/post-correction.molecules.cram
export RAYON_NUM_THREADS=$threads

run_timed() {
  local label=$1
  shift
  /usr/bin/time -v -o "$OUT_DIR/logs/$label.time.txt" \
    "$@" >"$OUT_DIR/logs/$label.stdout.txt" 2>"$OUT_DIR/logs/$label.stderr.txt"
}

run_timed export-molecule-bam "$AIE_BIN" export-molecule-bam "$ARCHIVE" \
  --fai "$genome.fai" --out "$molecule_bam"
"$samtools_bin" quickcheck -v "$molecule_bam"

run_timed make-molecule-cram "$samtools_bin" view -@ "$threads" -T "$genome" -C \
  --output-fmt-option version=3.1 --output-fmt-option archive \
  --output-fmt-option use_lzma=1 -o "$molecule_cram" "$molecule_bam"
"$samtools_bin" quickcheck -v "$molecule_cram"

run_timed replay-molecule-bam "$AIE_BIN" replay-rows "$molecule_bam" \
  --from-molecule-bam --gtf "$gtf" --barcodes "$BARCODES" \
  --out-dir "$OUT_DIR/replay-bam" --solo-strand "$solo_strand"
for artifact in matrix.mtx barcodes.tsv features.tsv; do
  cmp "$OUT_DIR/replay-bam/$artifact" "$REFERENCE_MATRIX/$artifact"
done

# A textual record-stream hash proves that CRAM preserved every local tag, flag, CIGAR, and record
# in order. The separate streaming benchmark proves that the decoded stream is quantifiable.
bam_sam_sha=$("$samtools_bin" view -@ "$threads" "$molecule_bam" | sha256sum | cut -d' ' -f1)
cram_sam_sha=$("$samtools_bin" view -@ "$threads" -T "$genome" "$molecule_cram" \
  | sha256sum | cut -d' ' -f1)
[[ "$bam_sam_sha" == "$cram_sam_sha" ]] || {
  echo "CRAM record stream differs from molecule BAM" >&2
  exit 1
}

{
  printf 'artifact\tbytes\tsha256\n'
  for artifact in "$ARCHIVE" "$molecule_bam" "$molecule_cram"; do
    printf '%s\t%s\t%s\n' "$artifact" "$(stat -c%s "$artifact")" \
      "$(sha256sum "$artifact" | cut -d' ' -f1)"
  done
  printf 'normalized_sam_stream\tNA\t%s\n' "$bam_sam_sha"
} >"$OUT_DIR/artifact-metrics.tsv"

{
  printf 'key\tvalue\n'
  printf 'date\t%s\n' "$run_date"
  printf 'dataset\t%s\n' "$DATASET"
  printf 'gravlax_commit\t%s\n' "$(git -C "$gravlax_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'samtools\t%s\n' "$("$samtools_bin" --version | head -1)"
  printf 'threads\t%s\n' "$threads"
  printf 'solo_strand\t%s\n' "$solo_strand"
  printf 'cram_version\t3.1\n'
  printf 'cram_profile\tarchive,use_lzma=1\n'
  printf 'reference_sha256\t%s\n' "$(sha256sum "$genome" | cut -d' ' -f1)"
} >"$OUT_DIR/protocol.tsv"

printf 'post-correction storage gate passed: %s\n' "$OUT_DIR"
