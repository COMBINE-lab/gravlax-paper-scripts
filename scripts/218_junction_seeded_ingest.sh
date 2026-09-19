#!/usr/bin/env bash
# D0: does seeding the annotation-free ingest alignment with a junction list close the gap to
# annotation-aware STARsolo, and does a five-year-old junction list (v32) do as well as the current one (v49)?
# Arms: unseeded two-pass (existing archive), v32-seeded two-pass, v49-seeded two-pass, v32-seeded one-pass.
# Archived from the analysis host as scripts/153_junction_seeded_ingest.sh (2026-09-11); renumbered here to follow the repository sequence.
set -euo pipefail
P=${GRAVLAX_PROJECT_ROOT:?set GRAVLAX_PROJECT_ROOT to the analysis project root (see README)}
AIE=$P/runs/gravlax-release-0.2.2-20260910/aie-0.2.2
STAR=$P/env/sc/bin/STAR
PY=${GRAVLAX_PYTHON:-python3}
O=$P/runs/paper-followup-20260911/junction-seeded-ingest
R1=$P/data/fastq/pbmc_1k_v3_S1_R1_001.fastq.gz; R2=$P/data/fastq/pbmc_1k_v3_S1_R2_001.fastq.gz
WL=$P/data/3M-february-2018.txt
GTF49=$P/annotations/gencode.v49.annotation.gtf; GTF32=$P/annotations/gencode.v32.annotation.gtf
ORA49=$P/runs/oracle/full/v49/Solo.out; ORA32=$P/runs/oracle/full/v32/Solo.out
BARC=$ORA49/Gene/raw/barcodes.tsv; CELLS=$ORA49/Gene/filtered/barcodes.tsv
THREADS=24; CPUSET=0-23; export RAYON_NUM_THREADS=$THREADS
cd $O
# 1. junction lists in STAR sjdbFileChrStartEnd format (1-based intron start/end, strand)
for v in 32 49; do
  [[ -s v$v.sjdb.tab ]] || $PY - $P/annotations/gencode.v$v.annotation.gtf v$v.sjdb.tab <<'PYEOF'
import sys,collections
gtf,out=sys.argv[1],sys.argv[2]
tx=collections.defaultdict(list); meta={}
for line in open(gtf):
    if line[0]=="#": continue
    c=line.split("\t",9)
    if c[2]!="exon": continue
    a=c[8]; i=a.find('transcript_id "'); tid=a[i+15:a.find('"',i+15)]
    tx[tid].append((int(c[3]),int(c[4]))); meta[tid]=(c[0],c[6])
J=set()
for tid,ex in tx.items():
    ex=sorted(set(ex)); ch,st=meta[tid]
    for k in range(len(ex)-1):
        if ex[k+1][0]-1>=ex[k][1]+1: J.add((ch,ex[k][1]+1,ex[k+1][0]-1,st))
with open(out,"w") as f:
    for ch,a,b,st in sorted(J): f.write(f"{ch}\t{a}\t{b}\t{st}\n")
print(out,len(J))
PYEOF
done
align() {  # name seedfile twopass(yes/no)
  local name=$1 seed=$2 tp=$3; local d=$O/align-$name; mkdir -p $d
  local tpargs=(); [[ $tp == yes ]] && tpargs=(--twopassMode Basic)
  /usr/bin/time -v -o $d/time.txt taskset -c $CPUSET $STAR --runThreadN $THREADS --genomeDir $P/ref/star-base \
    --sjdbFileChrStartEnd $seed --sjdbOverhang 90 "${tpargs[@]}" \
    --readFilesIn $R2 $R1 --readFilesCommand zcat --soloType CB_UMI_Simple --soloCBwhitelist $WL \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
    --soloFeatures SJ --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --outSAMattributes NH HI AS nM CR UR CY UY --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 32000000000 \
    --outSAMunmapped Within --outFileNamePrefix $d/ > $d/stdout.txt 2> $d/stderr.txt
  rm -rf $d/_STARtmp
}
ingest() {  # name seedfile discovery catalogue
  local name=$1 seed=$2 disc=$3 cat=$4; local d=$O/align-$name
  /usr/bin/time -v -o $O/ingest-$name.time.txt taskset -c $CPUSET $AIE ingest-archive $d/Aligned.sortedByCoord.out.bam \
    --whitelist $WL --out $O/d0-$name.aie --zstd-level 19 --junction-discovery $disc --junction-catalogue $cat \
    --alignment-annotation $seed --alignment-index-identity "$P/ref/star-base (annotation-free STAR 2.7.11b index; seed inserted at mapping, sjdbOverhang 90)" \
    --alignment-chemistry 10x-3p-v3 --alignment-log $d/Log.out --report-format json --report-output $O/ingest-$name.report.json \
    > $O/ingest-$name.stdout.txt 2> $O/ingest-$name.stderr.txt
}
echo "[$(date)] align v32 two-pass"; align v32-2pass $O/v32.sjdb.tab yes
ingest v32-2pass $O/v32.sjdb.tab per-library-two-pass $O/align-v32-2pass/_STARpass1/SJ.out.tab
echo "[$(date)] align v49 two-pass"; align v49-2pass $O/v49.sjdb.tab yes
ingest v49-2pass $O/v49.sjdb.tab per-library-two-pass $O/align-v49-2pass/_STARpass1/SJ.out.tab
echo "[$(date)] align v32 one-pass"; align v32-1pass $O/v32.sjdb.tab no
ingest v32-1pass $O/v32.sjdb.tab frozen-catalogue $O/v32.sjdb.tab
cp $P/runs/archive/d0.aie $O/d0-unseeded.aie
for arm in unseeded v32-2pass v49-2pass v32-1pass; do
  A=$O/d0-$arm.aie
  echo "[$(date)] replay $arm"
  taskset -c $CPUSET $AIE replay-rows $A --gtf $GTF49 --barcodes $BARC --out-dir $O/replay-$arm-v49 > $O/replay-$arm-v49.log 2>&1
  taskset -c $CPUSET $AIE replay-rows $A --gtf $GTF32 --barcodes $ORA32/Gene/raw/barcodes.tsv --out-dir $O/replay-$arm-v32 > $O/replay-$arm-v32.log 2>&1
  taskset -c $CPUSET $AIE replay-rows $A --gtf $GTF49 --barcodes $BARC --out-dir $O/replay-$arm-velo --velocity > $O/replay-$arm-velo.log 2>&1
done
{
  for arm in unseeded v32-2pass v49-2pass v32-1pass; do
    echo "===== $arm  ($(stat -c %s $O/d0-$arm.aie) bytes)"
    [[ -f $O/align-$arm/Log.final.out ]] && grep -E "Uniquely mapped reads %|Number of splices: Total|Number of splices: Annotated" $O/align-$arm/Log.final.out
    echo "--- Gene v49 vs STARsolo v49"; $PY $P/scripts/30_gateb_compare.py $ORA49/Gene/raw $O/replay-$arm-v49 oracle replay --cells $CELLS 2>/dev/null | grep -E "UMI mass moved \(shared|expressed genes|cells changing"
    echo "--- Gene v32 vs STARsolo v32"; $PY $P/scripts/30_gateb_compare.py $ORA32/Gene/raw $O/replay-$arm-v32 oracle replay --cells $CELLS 2>/dev/null | grep -E "UMI mass moved \(shared|expressed genes"
    echo "--- velocity v49"; $PY $P/scripts/91_velo_gate.py $ORA49/Velocyto/raw $O/replay-$arm-velo $CELLS 2>/dev/null | grep mtx
  done
} > $O/summary.txt 2>&1
cat $O/summary.txt; echo "JUNCTION SEEDED INGEST DONE"
