set -x
P=${GRAVLAX_PROJECT_ROOT}
AIE=$P/src/target/release/aie
GTF=$P/annotations/gencode.v49.annotation.gtf
for D in d0 d1; do
  if [ $D = d0 ]; then
    BAM=$P/runs/ingest/full/Aligned.sortedByCoord.out.bam
    BARC=$P/runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv
    REF=$P/runs/replay/full-frombam
  else
    BAM=$P/runs/ingest/d1/Aligned.sortedByCoord.out.bam
    BARC=$P/runs/oracle/d1/v49/Solo.out/Gene/raw/barcodes.tsv
    REF=$P/runs/replay/d1-frombam
  fi
  /usr/bin/time -f "${D}_ingest_wall=%e" $AIE ingest-archive $BAM \
    --whitelist $P/data/3M-february-2018.txt --out $P/runs/archive/$D.aie --zstd-level 19 \
    2>&1 | grep -E 'extracted|ingest_wall' || { echo "${D} INGEST FAILED"; exit 1; }
  stat -c "${D}_SIZE %s" $P/runs/archive/$D.aie
  /usr/bin/time -f "${D}_replay_wall=%e" $AIE replay-rows $P/runs/archive/$D.aie \
    --gtf $GTF --barcodes $BARC --out-dir $P/runs/replay/$D-v13 2>&1 | grep -E 'molecules|replay_wall'
  ok=1
  for f in matrix.mtx barcodes.tsv features.tsv; do
    cmp -s $P/runs/replay/$D-v13/$f $REF/$f || { echo "${D} MISMATCH $f"; ok=0; }
  done
  [ $ok = 1 ] && echo "${D} REGRESSION: BYTE-IDENTICAL" || echo "${D} REGRESSION: FAILED"
  $AIE query $P/runs/archive/$D.aie region chr6:73489308-73525587 | head -1
  $AIE query $P/runs/archive/$D.aie junction chr16:89562391-89562883 | head -1
done
echo V14_PIPELINE_DONE
