#!/usr/bin/env bash
# The archive's own ingest alignment: annotation-free, but barcode-aware.
#
# This is the run the archive is built from, so it must touch no GTF. It still needs cell barcodes
# and UMIs, which plain STAR does not extract -- hence STARsolo with --soloFeatures SJ, the one
# feature that is defined without a gene model. Two-pass mode gives annotation-free alignment its
# best shot at the novel junctions that sjdb insertion would otherwise supply.
#
# Emits ONLY the raw CR/UR barcode and UMI. STARsolo refuses to write corrected CB/UB unless a
# gene feature (Gene/GeneFull/...) is requested, and those require a GTF -- so an annotation-free
# ingest cannot obtain STARsolo's corrected barcodes, and must correct them itself against the
# whitelist. That is the honest arrangement anyway: barcode correction is whitelist-based and
# genuinely annotation-independent, whereas STARsolo's UMI correction is performed per-gene and
# must not leak into our grouping (docs/decision_log.md D9, D11).
#
# Usage: scripts/40_ingest_align.sh [reads_tag]
set -euo pipefail
READS_TAG="${1:-full}"
: "${PROJ_ROOT:?source env/setup.sh first}"
cd "$PROJ_ROOT"

case "$READS_TAG" in
  full) R1=data/fastq/pbmc_1k_v3_S1_R1_001.fastq.gz; R2=data/fastq/pbmc_1k_v3_S1_R2_001.fastq.gz ;;
  dev)  R1=data/fastq/dev_R1.fastq.gz;               R2=data/fastq/dev_R2.fastq.gz ;;
  *) echo "unknown reads_tag: $READS_TAG" >&2; exit 2 ;;
esac

OUT="runs/ingest/${READS_TAG}"
mkdir -p "$OUT"

STAR \
  --runThreadN "${AIE_THREADS:-32}" \
  --genomeDir ref/star-base \
  --twopassMode Basic \
  --readFilesIn "$R2" "$R1" \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist data/3M-february-2018.txt \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloFeatures SJ \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --outSAMattributes NH HI AS nM CR UR CY UY \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 32000000000 \
  --outSAMunmapped Within \
  --outFileNamePrefix "$OUT/" \
  2>&1 | tee "logs/ingest-${READS_TAG}.log"

echo "annotation-free ingest alignment done -> $OUT"
