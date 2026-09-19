#!/usr/bin/env bash
# Build one STAR index per annotation, so splice-junction insertion is paid once rather than on
# every oracle run. Inserting GENCODE v49's 507k transcripts at mapping time costs ~30 min per run;
# there are many runs.
#
# The annotation-free base index (ref/star-base) is built separately by scripts/00 and is the one
# the archive itself ingests from. These annotated indices exist only to produce the fresh-processing
# oracle Q(D,A) that the archive is measured against.
set -euo pipefail
: "${PROJ_ROOT:?source env/setup.sh first}"
cd "$PROJ_ROOT"

FASTA=ref/GRCh38.primary_assembly.genome.fa
# sjdbOverhang = read length - 1. R2 is 91 bp for pbmc_1k_v3 10x 3' v3.
OVERHANG=90

for v in 49 32 48; do
  out="ref/star-${v}"
  if [[ -f "$out/SAindex" ]]; then echo "ref/star-${v} exists, skipping"; continue; fi
  mkdir -p "$out" "runs/star-index-${v}"
  echo "=== building annotated index for GENCODE v${v} ==="
  STAR --runMode genomeGenerate \
    --runThreadN "${AIE_THREADS:-32}" \
    --genomeDir "$out" \
    --genomeFastaFiles "$FASTA" \
    --sjdbGTFfile "annotations/gencode.v${v}.annotation.gtf" \
    --sjdbOverhang "$OVERHANG" \
    --outFileNamePrefix "runs/star-index-${v}/" \
    > "logs/star-index-${v}.log" 2>&1 &
done
wait
echo "annotated indices complete"
du -sh ref/star-*
