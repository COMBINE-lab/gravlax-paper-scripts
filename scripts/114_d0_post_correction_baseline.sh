#!/usr/bin/env bash
# D0 post-correction molecule BAM/CRAM storage baseline.
#
# This measures a standards-compliant alignment container carrying the exact Gravlax
# post-correction abstraction. It is not a generic read BAM: local tags encode opaque UMI classes,
# molecule/group structure, and explicit 1MM edges. See Gravlax docs-notes/molecule-bam.md.
#
# Required environment:
#   PROJ_ROOT   project root containing ref/, annotations/, and runs/
#   AIE_BIN     Gravlax binary containing export-molecule-bam
# Optional environment:
#   GRAVLAX_REPO source checkout used to build AIE_BIN (default: $PROJ_ROOT/work/gravlax)
#   SAMTOOLS    samtools binary (default: samtools on PATH)
#   OUT_DIR     new output directory
#   THREADS     samtools and Rayon threads (default: 24)
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT to the project root}"
: "${AIE_BIN:?set AIE_BIN to the post-v1 Gravlax binary}"
samtools_bin=${SAMTOOLS:-samtools}
threads=${THREADS:-24}
out_dir=${OUT_DIR:-$PROJ_ROOT/runs/post-v1/fair-baselines/d0-post-correction}
gravlax_repo=${GRAVLAX_REPO:-$PROJ_ROOT/work/gravlax}

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 2; }
[[ -x "$AIE_BIN" ]] || { echo "AIE_BIN is not executable: $AIE_BIN" >&2; exit 2; }
command -v "$samtools_bin" >/dev/null || { echo "samtools not found: $samtools_bin" >&2; exit 2; }
git -C "$gravlax_repo" rev-parse --verify HEAD >/dev/null 2>&1 || {
  echo "GRAVLAX_REPO is not a Git checkout: $gravlax_repo" >&2
  exit 2
}
[[ ! -e "$out_dir" ]] || { echo "refusing to overwrite existing OUT_DIR: $out_dir" >&2; exit 2; }

archive=$PROJ_ROOT/runs/archive/d0.aie
genome=$PROJ_ROOT/ref/GRCh38.primary_assembly.genome.fa
fai=$genome.fai
gtf=$PROJ_ROOT/annotations/gencode.v49.annotation.gtf
barcodes=$PROJ_ROOT/runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv
reference_matrix=$PROJ_ROOT/runs/replay/full-frombam

for input in "$archive" "$genome" "$fai" "$gtf" "$barcodes" \
  "$reference_matrix/matrix.mtx" "$reference_matrix/barcodes.tsv" \
  "$reference_matrix/features.tsv"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done

mkdir -p "$out_dir/logs" "$out_dir/replay-bam" "$out_dir/replay-cram-decoded"
molecule_bam=$out_dir/d0.post-correction.molecules.bam
molecule_cram=$out_dir/d0.post-correction.molecules.cram
decoded_bam=$out_dir/d0.post-correction.cram-decoded.bam
export RAYON_NUM_THREADS=$threads

run_timed() {
  local label=$1
  shift
  /usr/bin/time -v -o "$out_dir/logs/$label.time.txt" \
    "$@" >"$out_dir/logs/$label.stdout.txt" 2>"$out_dir/logs/$label.stderr.txt"
}

run_timed export-molecule-bam "$AIE_BIN" export-molecule-bam "$archive" \
  --fai "$fai" --out "$molecule_bam"
"$samtools_bin" quickcheck -v "$molecule_bam"

run_timed make-molecule-cram "$samtools_bin" view -@ "$threads" -T "$genome" -C \
  --output-fmt-option version=3.1 --output-fmt-option archive \
  --output-fmt-option use_lzma=1 -o "$molecule_cram" "$molecule_bam"
"$samtools_bin" quickcheck -v "$molecule_cram"

run_timed decode-molecule-cram "$samtools_bin" view -@ "$threads" -T "$genome" -b \
  -o "$decoded_bam" "$molecule_cram"
"$samtools_bin" quickcheck -v "$decoded_bam"

run_timed replay-molecule-bam "$AIE_BIN" replay-rows "$molecule_bam" \
  --from-molecule-bam --gtf "$gtf" --barcodes "$barcodes" \
  --out-dir "$out_dir/replay-bam"
run_timed replay-cram-decoded "$AIE_BIN" replay-rows "$decoded_bam" \
  --from-molecule-bam --gtf "$gtf" --barcodes "$barcodes" \
  --out-dir "$out_dir/replay-cram-decoded"

for replay_dir in replay-bam replay-cram-decoded; do
  for artifact in matrix.mtx barcodes.tsv features.tsv; do
    cmp "$out_dir/$replay_dir/$artifact" "$reference_matrix/$artifact"
  done
done

bam_sam_sha=$("$samtools_bin" view -@ "$threads" "$molecule_bam" | sha256sum | cut -d' ' -f1)
cram_sam_sha=$("$samtools_bin" view -@ "$threads" -T "$genome" "$molecule_cram" \
  | sha256sum | cut -d' ' -f1)
[[ "$bam_sam_sha" == "$cram_sam_sha" ]] || {
  echo "CRAM record stream differs from molecule BAM" >&2
  exit 1
}

{
  printf 'artifact\tbytes\tsha256\n'
  for artifact in "$archive" "$molecule_bam" "$molecule_cram" "$decoded_bam"; do
    printf '%s\t%s\t%s\n' "$artifact" "$(stat -c%s "$artifact")" \
      "$(sha256sum "$artifact" | cut -d' ' -f1)"
  done
  printf 'normalized_sam_stream\tNA\t%s\n' "$bam_sam_sha"
} >"$out_dir/artifact-metrics.tsv"

{
  printf 'gravlax_commit\t%s\n' "$(git -C "$gravlax_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'samtools\t%s\n' "$("$samtools_bin" --version | head -1)"
  printf 'threads\t%s\n' "$threads"
  printf 'cram_version\t3.1\n'
  printf 'cram_profile\tarchive,use_lzma=1\n'
  printf 'reference_sha256\t%s\n' "$(sha256sum "$genome" | cut -d' ' -f1)"
} >"$out_dir/tool-and-input-identities.tsv"

printf 'D0 post-correction baseline passed: %s\n' "$out_dir"
