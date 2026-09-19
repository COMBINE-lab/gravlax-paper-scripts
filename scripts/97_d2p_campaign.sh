#!/usr/bin/env bash
# D2 campaign: Parent_SC3v3_Human_Glioblastoma — gates in docs/decisions/2026-08-30-d2p-gates.md.
set -uo pipefail
: "${PROJ_ROOT:?source env/setup.sh first}"
cd "$PROJ_ROOT"
AIE=src/target/release/aie
PY=${GRAVLAX_PYTHON:-python3}
D2=$(dirname "$(find data/d2p -name '*_R2_*.fastq.gz' | head -1)")
R1=$(ls $D2/*_R1_*.fastq.gz | paste -sd,)
R2=$(ls $D2/*_R2_*.fastq.gz | paste -sd,)
mkdir -p runs/ingest/d2p runs/oracle/d2p/v49 runs/oracle/d2p/v32

# Stage 1: three STAR runs — annotation-free ingest, v49 oracle, v32 oracle (D2-SIGNAL).
(
  STAR --runThreadN 20 --genomeDir ref/star-base --twopassMode Basic \
    --readFilesIn "$R2" "$R1" --readFilesCommand zcat \
    --soloType CB_UMI_Simple --soloCBwhitelist data/3M-february-2018.txt \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 --soloFeatures SJ \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --outSAMattributes NH HI AS nM CR UR CY UY \
    --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 32000000000 \
    --outSAMunmapped Within --outFileNamePrefix runs/ingest/d2p/ \
    > logs/d2p-ingest.log 2>&1; echo $? > runs/ingest/d2p/EXITCODE
) &
(
  STAR --runThreadN 20 --genomeDir ref/star-49 \
    --readFilesIn "$R2" "$R1" --readFilesCommand zcat \
    --soloType CB_UMI_Simple --soloCBwhitelist data/3M-february-2018.txt \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 --soloFeatures Gene GeneFull SJ Velocyto \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIdedup 1MM_CR --soloUMIfiltering MultiGeneUMI_CR \
    --soloCellFilter EmptyDrops_CR \
    --outSAMattributes NH HI AS nM CR UR CB UB GX GN \
    --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 80000000000 \
    --outSAMunmapped Within --outFileNamePrefix runs/oracle/d2p/v49/ \
    > logs/d2p-oracle-v49.log 2>&1; echo $? > runs/oracle/d2p/v49/EXITCODE
) &
wait
for f in runs/ingest/d2p/EXITCODE runs/oracle/d2p/v49/EXITCODE; do
  [[ "$(cat $f)" = 0 ]] || { echo "D2P STAGE1 FAILED: $f"; exit 1; }
done
# v32 oracle (matrix-only) after the pair, to keep the thread budget polite.
STAR --runThreadN 24 --genomeDir ref/star-32 \
  --readFilesIn "$R2" "$R1" --readFilesCommand zcat \
  --soloType CB_UMI_Simple --soloCBwhitelist data/3M-february-2018.txt \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 --soloFeatures Gene \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIdedup 1MM_CR --soloUMIfiltering MultiGeneUMI_CR \
  --soloCellFilter EmptyDrops_CR \
  --outSAMtype None --outFileNamePrefix runs/oracle/d2p/v32/ \
  > logs/d2p-oracle-v32.log 2>&1
echo "D2P_STAGE1_DONE"

BAM=runs/ingest/d2p/Aligned.sortedByCoord.out.bam
samtools view -@ 24 -T ref/GRCh38.primary_assembly.genome.fa -C \
  --output-fmt-option version=3.1 --output-fmt-option archive --output-fmt-option use_lzma=1 \
  -o runs/gatee/d2p_ingest_annotation_free.cram "$BAM"
/usr/bin/time -f "d2p_aie_ingest_wall=%e" $AIE ingest-archive "$BAM" \
  --whitelist data/3M-february-2018.txt --out runs/archive/d2p.aie --zstd-level 19 \
  2>&1 | grep -E 'extracted|wall'
echo "D2P_STAGE2_DONE"

BARC=runs/oracle/d2p/v49/Solo.out/Gene/raw/barcodes.tsv
CELLS=runs/oracle/d2p/v49/Solo.out/Gene/filtered/barcodes.tsv
/usr/bin/time -f "d2p_replay_wall=%e" $AIE replay-rows runs/archive/d2p.aie \
  --gtf annotations/gencode.v49.annotation.gtf --barcodes "$BARC" \
  --out-dir runs/replay/d2p-aie 2>&1 | grep -E 'molecules|wall'
$AIE replay-rows "$BAM" --from-bam --whitelist data/3M-february-2018.txt \
  --gtf annotations/gencode.v49.annotation.gtf --barcodes "$BARC" \
  --out-dir runs/replay/d2p-frombam > /dev/null 2>&1
ok=1
for f in matrix.mtx barcodes.tsv features.tsv; do
  cmp -s runs/replay/d2p-aie/$f runs/replay/d2p-frombam/$f || ok=0
done
[[ $ok = 1 ]] && echo "D2P-REG: BYTE-IDENTICAL" || echo "D2P-REG: FAILED"
$PY scripts/30_gateb_compare.py runs/oracle/d2p/v49/Solo.out/Gene/raw runs/replay/d2p-aie \
  oracle replay --cells "$CELLS" --out results/d2p-fid.json | grep -E 'moved|cells'
$PY scripts/30_gateb_compare.py runs/oracle/d2p/v49/Solo.out/Gene/raw runs/oracle/d2p/v32/Solo.out/Gene/raw \
  v49 v32 --cells "$CELLS" --out results/d2p-signal.json | grep -E 'moved|genes chang'
$AIE replay-rows runs/archive/d2p.aie --gtf annotations/gencode.v49.annotation.gtf \
  --barcodes runs/oracle/d2p/v49/Solo.out/Velocyto/raw/barcodes.tsv \
  --out-dir runs/replay/d2p-velocity --velocity > /dev/null 2>&1
$PY scripts/91_velo_gate.py runs/oracle/d2p/v49/Solo.out/Velocyto/raw runs/replay/d2p-velocity "$CELLS" | grep -E 'mtx|VELO'
$AIE em runs/archive/d2p.aie --gtf annotations/gencode.v49.annotation.gtf | grep -E 'classes|mode'
echo "---- D2 storage ----"
for f in $D2/*_R[12]_*.fastq.gz; do stat -c%s "$f"; done | awk '{s+=$1} END {print "fastq\t" s}'
stat -c 'bam	%s' "$BAM"; stat -c 'cram31	%s' runs/gatee/d2p_ingest_annotation_free.cram
stat -c 'aie	%s' runs/archive/d2p.aie
echo "D2P_CAMPAIGN_DONE"
