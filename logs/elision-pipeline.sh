set -x
P=${GRAVLAX_PROJECT_ROOT}
AIE=$P/src/target/release/aie
$AIE ingest-archive $P/runs/ingest/full/Aligned.sortedByCoord.out.bam \
  --whitelist $P/data/3M-february-2018.txt --out $P/runs/archive/d0.aie --zstd-level 19 || exit 1
stat -c 'NEW SIZE %s' $P/runs/archive/d0.aie
$AIE replay-rows $P/runs/archive/d0.aie --gtf $P/annotations/gencode.v49.annotation.gtf \
  --barcodes $P/runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv \
  --out-dir $P/runs/replay/full-aie-v12 || exit 1
ok=1
for f in matrix.mtx barcodes.tsv features.tsv; do
  cmp -s $P/runs/replay/full-aie-v12/$f $P/runs/replay/full-frombam/$f || { echo "MISMATCH $f"; ok=0; }
done
[ $ok = 1 ] && echo "REGRESSION: BYTE-IDENTICAL" || echo "REGRESSION: FAILED"
$AIE query $P/runs/archive/d0.aie region chr6:73489308-73525587 | head -2
$AIE query $P/runs/archive/d0.aie junction chr16:89562391-89562883 | head -2
$AIE query $P/runs/archive/d0.aie apa chr6:73489308-73525587 --strand + --tsv | head -3
$AIE debug $P/runs/archive/d0.aie 2>&1 | sed -n '1,30p'
echo ELISION_PIPELINE_DONE
