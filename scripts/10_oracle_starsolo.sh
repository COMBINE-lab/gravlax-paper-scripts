#!/usr/bin/env bash
# Fresh-processing oracle Q(D, A): STARsolo under one annotation.
#
# Emits both the count matrices AND a BAM carrying CB/UB/GX/GN per read, so the oracle is
# available at MOLECULE level, not just as a count matrix. Gate -1 compares molecule by molecule,
# and Gate G-D needs STARsolo's own per-(cell,gene) UMI collapse (the UB tag) to measure the error
# of annotation-independent grouping against it.
#
# Usage: scripts/10_oracle_starsolo.sh <anno_id> <gtf.gz> [reads_tag]
#   anno_id   e.g. v49 | v32 | v48
#   gtf.gz    path to the GENCODE GTF for this annotation
#   reads_tag "full" (default) or "dev" (subsampled, for minute-scale iteration)
set -euo pipefail

ANNO_ID="${1:?usage: $0 <anno_id> <gtf.gz> [reads_tag]}"
GTF="${2:?}"
READS_TAG="${3:-full}"

: "${PROJ_ROOT:?source env/setup.sh first}"
cd "$PROJ_ROOT"

case "$READS_TAG" in
  full) R1=data/fastq/pbmc_1k_v3_S1_R1_001.fastq.gz; R2=data/fastq/pbmc_1k_v3_S1_R2_001.fastq.gz ;;
  dev)  R1=data/fastq/dev_R1.fastq.gz;               R2=data/fastq/dev_R2.fastq.gz ;;
  *) echo "unknown reads_tag: $READS_TAG" >&2; exit 2 ;;
esac

OUT="runs/oracle/${READS_TAG}/${ANNO_ID}"
mkdir -p "$OUT"

# The GTF must be uncompressed for STAR's --sjdbGTFfile.
GTF_PLAIN="annotations/$(basename "${GTF%.gz}")"
[[ -f "$GTF_PLAIN" ]] || unpigz -p 8 -c "$GTF" > "$GTF_PLAIN"

# Prefer the pre-built annotated index (scripts/05). Inserting v49's 507k transcripts at mapping
# time costs ~30 min per run, and there are many runs. Fall back to on-the-fly insertion off the
# annotation-free base index if the annotated one has not been built.
GENOME_DIR="ref/star-${ANNO_ID#v}"
if [[ -f "$GENOME_DIR/SAindex" ]]; then
  SJDB_ARGS=()
  echo "using pre-built annotated index $GENOME_DIR"
else
  GENOME_DIR=ref/star-base
  SJDB_ARGS=(--sjdbGTFfile "$GTF_PLAIN")
  echo "annotated index not found; inserting junctions on the fly from $GTF_PLAIN"
fi

# 10x 3' v3 geometry: R1 = 16 bp CB + 12 bp UMI. cDNA read is given FIRST to --readFilesIn.
# soloCBmatchWLtype/soloUMIdedup/soloUMIfiltering/soloCellFilter are the CellRanger-equivalent
# settings; they are frozen here and must not be varied between annotations, or the A1-vs-A2
# comparison would confound annotation change with policy change.
STAR \
  --runThreadN "${AIE_THREADS:-32}" \
  --genomeDir "$GENOME_DIR" \
  "${SJDB_ARGS[@]}" \
  --readFilesIn "$R2" "$R1" \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist data/3M-february-2018.txt \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloFeatures Gene GeneFull SJ Velocyto \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIdedup 1MM_CR \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloCellFilter EmptyDrops_CR \
  --outSAMattributes NH HI AS nM CR UR CB UB GX GN \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 80000000000 \
  --outSAMunmapped Within \
  --outFileNamePrefix "$OUT/" \
  2>&1 | tee "logs/oracle-${READS_TAG}-${ANNO_ID}.log"

echo "oracle $ANNO_ID ($READS_TAG) done -> $OUT"
