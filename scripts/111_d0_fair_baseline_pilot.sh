#!/usr/bin/env bash
# Reproduce the post-v1 D0 ingest-equivalent storage pilot.
#
# This is a correctness and protocol pilot, not the randomized five-replicate
# performance experiment.  It fails rather than overwriting an existing run.
#
# Required environment:
#   PROJ_ROOT   project root containing data/, ref/, annotations/, and runs/
# Optional environment:
#   AIE_BIN     Gravlax binary (default: $PROJ_ROOT/src/target/release/aie)
#   SAMTOOLS    samtools binary (default: samtools on PATH)
#   PILOT_OUT   new output directory
#   THREADS     samtools and Rayon threads (default: 24)
set -euo pipefail

: "${PROJ_ROOT:?source environment/setup.sh or set PROJ_ROOT explicitly}"
aie_bin=${AIE_BIN:-$PROJ_ROOT/src/target/release/aie}
samtools_bin=${SAMTOOLS:-samtools}
threads=${THREADS:-24}
out_dir=${PILOT_OUT:-$PROJ_ROOT/runs/post-v1/fair-baselines/d0-reproduction}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 2; }
[[ -x "$aie_bin" ]] || { echo "AIE_BIN is not executable: $aie_bin" >&2; exit 2; }
command -v "$samtools_bin" >/dev/null || { echo "samtools not found: $samtools_bin" >&2; exit 2; }
[[ ! -e "$out_dir" ]] || { echo "refusing to overwrite existing PILOT_OUT: $out_dir" >&2; exit 2; }

full_bam=$PROJ_ROOT/runs/ingest/full/Aligned.sortedByCoord.out.bam
full_cram=$PROJ_ROOT/runs/gatee/ingest_annotation_free.cram
canonical_aie=$PROJ_ROOT/runs/archive/d0.aie
whitelist=$PROJ_ROOT/data/3M-february-2018.txt
genome=$PROJ_ROOT/ref/GRCh38.primary_assembly.genome.fa
gtf=$PROJ_ROOT/annotations/gencode.v49.annotation.gtf
barcodes=$PROJ_ROOT/runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv
reference_matrix=$PROJ_ROOT/runs/replay/full-frombam

for input in "$full_bam" "$full_cram" "$canonical_aie" "$whitelist" "$genome" "$gtf" \
  "$barcodes" "$reference_matrix/matrix.mtx" "$reference_matrix/barcodes.tsv" \
  "$reference_matrix/features.tsv"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done

mkdir -p "$out_dir/logs" "$out_dir/replay-archive" "$out_dir/replay-full-bam" \
  "$out_dir/replay-minimal-bam"
minimal_bam=$out_dir/ingest.replay-minimal.bam
minimal_cram=$out_dir/ingest.replay-minimal.cram
full_aie=$out_dir/d0.full-reingest.aie
minimal_aie=$out_dir/d0.replay-minimal.aie
stamped_aie=$out_dir/d0.replay-minimal.stamped.aie

export RAYON_NUM_THREADS=$threads

run_timed() {
  local label=$1
  shift
  /usr/bin/time -v -o "$out_dir/logs/$label.time.txt" \
    "$@" >"$out_dir/logs/$label.stdout.txt" 2>"$out_dir/logs/$label.stderr.txt"
}

run_timed make-minimal env SAMTOOLS="$samtools_bin" \
  "$script_dir/110_make_replay_minimal_bam.sh" "$full_bam" "$minimal_bam" "$threads"

run_timed make-minimal-cram "$samtools_bin" view -@ "$threads" -T "$genome" -C \
  --output-fmt-option version=3.1 --output-fmt-option archive \
  --output-fmt-option use_lzma=1 -o "$minimal_cram" "$minimal_bam"
"$samtools_bin" quickcheck -v "$minimal_bam" "$minimal_cram"

run_timed ingest-full "$aie_bin" ingest-archive "$full_bam" --whitelist "$whitelist" \
  --out "$full_aie" --zstd-level 19
run_timed ingest-minimal "$aie_bin" ingest-archive "$minimal_bam" --whitelist "$whitelist" \
  --out "$minimal_aie" --zstd-level 19
cmp "$full_aie" "$minimal_aie"

run_timed stamp-minimal "$aie_bin" stamp-genome "$minimal_aie" --genome "$genome" \
  --out "$stamped_aie"
cmp "$canonical_aie" "$stamped_aie"

run_timed replay-archive "$aie_bin" replay-rows "$minimal_aie" --gtf "$gtf" \
  --barcodes "$barcodes" --out-dir "$out_dir/replay-archive"
run_timed replay-minimal-bam "$aie_bin" replay-rows "$minimal_bam" --from-bam \
  --whitelist "$whitelist" --gtf "$gtf" --barcodes "$barcodes" \
  --out-dir "$out_dir/replay-minimal-bam"
run_timed replay-full-bam "$aie_bin" replay-rows "$full_bam" --from-bam \
  --whitelist "$whitelist" --gtf "$gtf" --barcodes "$barcodes" \
  --out-dir "$out_dir/replay-full-bam"

for replay_dir in replay-archive replay-minimal-bam replay-full-bam; do
  for artifact in matrix.mtx barcodes.tsv features.tsv; do
    cmp "$out_dir/$replay_dir/$artifact" "$reference_matrix/$artifact"
  done
done

bam_sam_sha=$("$samtools_bin" view -@ "$threads" "$minimal_bam" | sha256sum | cut -d' ' -f1)
cram_sam_sha=$("$samtools_bin" view -@ "$threads" -T "$genome" "$minimal_cram" \
  | sha256sum | cut -d' ' -f1)
[[ "$bam_sam_sha" == "$cram_sam_sha" ]] || { echo "CRAM record stream differs from BAM" >&2; exit 1; }

{
  printf 'artifact\tbytes\tsha256\n'
  for artifact in "$full_bam" "$full_cram" "$minimal_bam" "$minimal_cram" \
    "$canonical_aie" "$minimal_aie" "$stamped_aie"; do
    printf '%s\t%s\t%s\n' "$artifact" "$(stat -c%s "$artifact")" \
      "$(sha256sum "$artifact" | cut -d' ' -f1)"
  done
  printf 'normalized_sam_stream\tNA\t%s\n' "$bam_sam_sha"
} >"$out_dir/artifact-metrics.tsv"

printf 'D0 fair-baseline pilot passed: %s\n' "$out_dir"

