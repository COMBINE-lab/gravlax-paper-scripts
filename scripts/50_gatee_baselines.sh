#!/usr/bin/env bash
# Gate G-E baselines: the sizes the archive must beat.
#
# The competitor that matters is NOT the count matrix -- it is the barcoded BAM, which already
# stores intervals, blocks, junctions, edits and corrected CB/UB, and which people already
# re-quantify against new GTFs today (STARsolo --readFilesType SAM SE, featureCounts --byReadGroup,
# cellCounts). See docs/prior-art.md §4. So CRAM 3.1 of that BAM, with tags retained, is the
# primary baseline, and `--capture-all-tags` is not optional: a CRAM that silently dropped CB/UB
# would be an unfairly small number to compare against.
#
# Usage: scripts/50_gatee_baselines.sh [reads_tag]
set -euo pipefail
READS_TAG="${1:-full}"
: "${PROJ_ROOT:?source env/setup.sh first}"
cd "$PROJ_ROOT"

FASTA=ref/GRCh38.primary_assembly.genome.fa
OUT=results/gatee-baselines-${READS_TAG}.tsv
mkdir -p results runs/gatee
THREADS="${AIE_THREADS:-16}"

# The archive ingests from the annotation-free alignment, so that is the BAM whose CRAM is the
# like-for-like baseline. The oracle BAM is also measured: it is what a user would actually have
# kept, and it is the number a reviewer will quote.
declare -A BAMS=(
  [ingest_annotation_free]="runs/ingest/${READS_TAG}/Aligned.sortedByCoord.out.bam"
  [oracle_v49]="runs/oracle/${READS_TAG}/v49/Aligned.sortedByCoord.out.bam"
)

{
  printf "artifact\tbytes\tnote\n"
  for f in data/fastq/*.fastq.gz; do
    [[ "$f" == *dev_* && "$READS_TAG" == full ]] && continue
    [[ "$f" != *dev_* && "$READS_TAG" == dev  ]] && continue
    printf "fastq_gz:%s\t%s\traw reads\n" "$(basename "$f")" "$(stat -c%s "$f")"
  done

  for name in "${!BAMS[@]}"; do
    bam="${BAMS[$name]}"
    [[ -f "$bam" ]] || { echo "missing $bam, skipping" >&2; continue; }
    printf "bam:%s\t%s\tuncompressed-tag BAM\n" "$name" "$(stat -c%s "$bam")"

    cram="runs/gatee/${name}.cram"
    if [[ ! -f "$cram" ]]; then
      samtools view -@ "$THREADS" -T "$FASTA" -C \
        --output-fmt-option version=3.1 \
        --output-fmt-option archive \
        --output-fmt-option use_lzma=1 \
        -o "$cram" "$bam"
    fi
    printf "cram31:%s\t%s\tCRAM 3.1 archive mode, all tags retained\n" "$name" "$(stat -c%s "$cram")"
  done
} | tee "$OUT"

echo
echo "wrote $OUT"
echo "Reference points measured earlier on this same dataset:"
echo "  Malva index (pbmc_1k_v3): 393,216,000 bytes approx (375 MB) -- ${MALVA_ROOT}/indices/pbmc_1k_v3"
