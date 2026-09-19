#!/usr/bin/env bash
# D1 re-verification campaign: every headline number, re-measured on pbmc_5k_v3 (project rule:
# single-dataset numbers do not extrapolate — five dev!=full instances say so, D26/D27 numbers
# included). Stage 1 runs the two STAR alignments in parallel (annotation-free ingest + v49
# oracle); stage 2 builds the CRAM 3.1 baseline and the .aie v1.2 archive; stage 3 replays,
# checks byte-identity, measures fidelity vs the fresh oracle, and prints the storage ladder.
set -uo pipefail
: "${PROJ_ROOT:?source env/setup.sh first}"
cd "$PROJ_ROOT"
AIE=src/target/release/aie
D1=data/d1/5k_pbmc_v3_fastqs
R1=$(ls $D1/*_R1_*.fastq.gz | paste -sd,)
R2=$(ls $D1/*_R2_*.fastq.gz | paste -sd,)

mkdir -p runs/ingest/d1 runs/oracle/d1/v49 runs/gatee logs

# ---- Stage 1: two STAR runs, parallel, 24 threads each ----
(
  STAR --runThreadN 24 --genomeDir ref/star-base --twopassMode Basic \
    --readFilesIn "$R2" "$R1" --readFilesCommand zcat \
    --soloType CB_UMI_Simple --soloCBwhitelist data/3M-february-2018.txt \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 --soloFeatures SJ \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --outSAMattributes NH HI AS nM CR UR CY UY \
    --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 32000000000 \
    --outSAMunmapped Within --outFileNamePrefix runs/ingest/d1/ \
    > logs/d1-ingest.log 2>&1
  echo $? > runs/ingest/d1/EXITCODE
) &
ING_PID=$!
(
  /usr/bin/time -f "oracle_wall=%e" STAR --runThreadN 24 --genomeDir ref/star-49 \
    --readFilesIn "$R2" "$R1" --readFilesCommand zcat \
    --soloType CB_UMI_Simple --soloCBwhitelist data/3M-february-2018.txt \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 --soloFeatures Gene GeneFull SJ Velocyto \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIdedup 1MM_CR --soloUMIfiltering MultiGeneUMI_CR \
    --soloCellFilter EmptyDrops_CR \
    --outSAMattributes NH HI AS nM CR UR CB UB GX GN \
    --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 80000000000 \
    --outSAMunmapped Within --outFileNamePrefix runs/oracle/d1/v49/ \
    > logs/d1-oracle-v49.log 2>&1
  echo $? > runs/oracle/d1/v49/EXITCODE
) &
ORA_PID=$!
wait $ING_PID; wait $ORA_PID
for f in runs/ingest/d1/EXITCODE runs/oracle/d1/v49/EXITCODE; do
  [[ "$(cat $f)" = 0 ]] || { echo "STAGE1 FAILED: $f = $(cat $f)"; exit 1; }
done
echo "STAGE1 DONE"

BAM=runs/ingest/d1/Aligned.sortedByCoord.out.bam

# ---- Stage 2: CRAM 3.1 baseline + archive build ----
samtools view -@ 24 -T ref/GRCh38.primary_assembly.genome.fa -C \
  --output-fmt-option version=3.1 --output-fmt-option archive \
  --output-fmt-option use_lzma=1 \
  -o runs/gatee/d1_ingest_annotation_free.cram "$BAM" || { echo "CRAM FAILED"; exit 1; }
/usr/bin/time -f "aie_ingest_wall=%e" $AIE ingest-archive "$BAM" \
  --whitelist data/3M-february-2018.txt --out runs/archive/d1.aie --zstd-level 19 \
  2>&1 | tee logs/d1-aie-ingest.log
[[ -s runs/archive/d1.aie ]] || { echo "AIE INGEST FAILED"; exit 1; }
echo "STAGE2 DONE"

# ---- Stage 3: replay, regression, fidelity, storage ladder ----
BARC=runs/oracle/d1/v49/Solo.out/Gene/raw/barcodes.tsv
/usr/bin/time -f "replay_wall=%e" $AIE replay-rows runs/archive/d1.aie \
  --gtf annotations/gencode.v49.annotation.gtf --barcodes "$BARC" \
  --out-dir runs/replay/d1-aie 2>&1 | tee logs/d1-replay.log
$AIE replay-rows "$BAM" --from-bam --whitelist data/3M-february-2018.txt \
  --gtf annotations/gencode.v49.annotation.gtf --barcodes "$BARC" \
  --out-dir runs/replay/d1-frombam > logs/d1-replay-frombam.log 2>&1
ok=1
for f in matrix.mtx barcodes.tsv features.tsv; do
  cmp -s runs/replay/d1-aie/$f runs/replay/d1-frombam/$f || { echo "MISMATCH $f"; ok=0; }
done
[[ $ok = 1 ]] && echo "REGRESSION: BYTE-IDENTICAL" || echo "REGRESSION: FAILED"

PY=${PY:-${GRAVLAX_PYTHON:-python3}}
"$PY" scripts/30_gateb_compare.py runs/oracle/d1/v49/Solo.out/Gene/raw runs/replay/d1-aie \
  oracle replay --cells runs/oracle/d1/v49/Solo.out/Gene/filtered/barcodes.tsv \
  --out results/gate1a-d1-v49.json
grep -E 'umi_mass_moved_frac|cells_changed|cells_live|total_umi' results/gate1a-d1-v49.json

echo "---- storage ladder (bytes) ----"
for f in $D1/*_R[12]_*.fastq.gz; do stat -c%s "$f"; done | awk '{s+=$1} END {print "fastq_R1R2\t" s}'
stat -c 'bam	%s' "$BAM"
stat -c 'cram31	%s' runs/gatee/d1_ingest_annotation_free.cram
stat -c 'aie	%s' runs/archive/d1.aie
$AIE debug runs/archive/d1.aie 2>/dev/null | sed -n '/sections/,/^== chunk streams/p'
echo "D1_REVERIFY_DONE"
