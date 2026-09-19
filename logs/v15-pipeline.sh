set -x
P=${GRAVLAX_PROJECT_ROOT}
AIE=$P/src/target/release/aie
GTF=$P/annotations/gencode.v49.annotation.gtf
declare -A BAM BARC REF
BAM[d0]=$P/runs/ingest/full/Aligned.sortedByCoord.out.bam;  BARC[d0]=$P/runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv; REF[d0]=$P/runs/replay/full-frombam
BAM[d1]=$P/runs/ingest/d1/Aligned.sortedByCoord.out.bam;    BARC[d1]=$P/runs/oracle/d1/v49/Solo.out/Gene/raw/barcodes.tsv;   REF[d1]=$P/runs/replay/d1-frombam
BAM[d2]=$P/runs/ingest/d2/Aligned.sortedByCoord.out.bam;    BARC[d2]=$P/runs/oracle/d2/v49/Solo.out/Gene/raw/barcodes.tsv;   REF[d2]=$P/runs/replay/d2-frombam
BAM[d2p]=$P/runs/ingest/d2p/Aligned.sortedByCoord.out.bam;  BARC[d2p]=$P/runs/oracle/d2p/v49/Solo.out/Gene/raw/barcodes.tsv; REF[d2p]=$P/runs/replay/d2p-frombam
for D in d0 d1 d2 d2p; do
  $AIE ingest-archive ${BAM[$D]} --whitelist $P/data/3M-february-2018.txt \
    --out $P/runs/archive/$D.aie --zstd-level 19 > /dev/null 2>&1 || { echo "$D INGEST FAILED"; exit 1; }
  stat -c "${D}_SIZE %s" $P/runs/archive/$D.aie
  $AIE replay-rows $P/runs/archive/$D.aie --gtf $GTF --barcodes ${BARC[$D]} \
    --out-dir $P/runs/replay/$D-v15 > /dev/null 2>&1
  ok=1
  for f in matrix.mtx barcodes.tsv features.tsv; do
    cmp -s $P/runs/replay/$D-v15/$f ${REF[$D]}/$f || ok=0
  done
  [ $ok = 1 ] && echo "$D REGRESSION: BYTE-IDENTICAL" || echo "$D REGRESSION: FAILED"
done
$AIE query $P/runs/archive/d2p.aie region chr6:73489308-73525587 | head -1
echo V15_PIPELINE_DONE
