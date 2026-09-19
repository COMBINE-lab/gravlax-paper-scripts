#!/usr/bin/env bash
# Gate G-C: alignment invariance under annotation change.
#
# The archive is built from an annotation-free alignment, but the oracle it must reproduce uses an
# annotation-aware one (STAR inserts the GTF's junctions into the index at mapping time). If those
# two disagree on where reads land, the archive has a fidelity ceiling that no amount of stored
# evidence can repair. This measures that gap directly, and it is the cheapest way to kill the
# project if it deserves killing.
#
# Barcodes are irrelevant to where a read aligns, so these are plain STAR runs on the cDNA read
# only. Config (i) uses two-pass mode: that is the annotation-free ingest's best shot at
# discovering the junctions sjdb insertion would otherwise hand it.
#
# Usage: scripts/20_gatec_align_configs.sh [reads_tag]
set -euo pipefail

READS_TAG="${1:-full}"
: "${PROJ_ROOT:?source env/setup.sh first}"
cd "$PROJ_ROOT"

case "$READS_TAG" in
  full) R2=data/fastq/pbmc_1k_v3_S1_R2_001.fastq.gz ;;
  dev)  R2=data/fastq/dev_R2.fastq.gz ;;
  *) echo "unknown reads_tag: $READS_TAG" >&2; exit 2 ;;
esac

# Read names are the join key across configs, so all three must see the same reads in the same
# order. --outSAMattributes includes jM/jI: STAR's own junction annotation-status and coordinates,
# which lets the comparator distinguish "novel junction found" from "annotated junction used".
run_cfg () {
  local tag="$1"; local genome="$2"; shift 2
  local out="runs/gatec/${READS_TAG}/${tag}"
  mkdir -p "$out"
  echo "=== G-C config ${tag} (genome ${genome}) ==="
  STAR \
    --runThreadN "${AIE_THREADS:-32}" \
    --genomeDir "$genome" \
    --readFilesIn "$R2" \
    --readFilesCommand zcat \
    --outSAMattributes NH HI AS nM jM jI \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --outFileNamePrefix "$out/" \
    "$@" \
    2>&1 | tee "logs/gatec-${READS_TAG}-${tag}.log"
}

# (i) annotation-free ingest, two-pass novel-junction discovery. This is the archive's own
#     ingest path, so it gets the best shot the annotation-free setting allows.
run_cfg noanno_2pass ref/star-base --twopassMode Basic

# (i-b) annotation-free, single pass. Isolates how much of any gap two-pass discovery closes.
run_cfg noanno_1pass ref/star-base

# (ii) annotation-aware, A1 = GENCODE v49
run_cfg sjdb_v49 ref/star-49

# (iii) annotation-aware, A2-far = GENCODE v32
run_cfg sjdb_v32 ref/star-32

echo "G-C alignment configs complete for reads_tag=${READS_TAG}"
