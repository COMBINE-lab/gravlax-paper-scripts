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
  /usr/bin/time -f "${D}_aie_replay_wall=%e" $AIE replay-rows $P/runs/archive/$D.aie \
    --gtf $GTF --barcodes $BARC --out-dir $P/runs/replay/$D-mt 2>&1 | grep -E 'molecules|wall'
  /usr/bin/time -f "${D}_frombam_wall=%e" $AIE replay-rows $BAM --from-bam \
    --whitelist $P/data/3M-february-2018.txt --gtf $GTF --barcodes $BARC \
    --out-dir $P/runs/replay/$D-mt-frombam 2>&1 | grep -E 'molecules|wall'
  ok=1
  for f in matrix.mtx barcodes.tsv features.tsv; do
    cmp -s $P/runs/replay/$D-mt/$f $REF/$f || ok=0
    cmp -s $P/runs/replay/$D-mt-frombam/$f $REF/$f || ok=0
  done
  [ $ok = 1 ] && echo "${D} REGRESSION: BYTE-IDENTICAL" || echo "${D} REGRESSION: FAILED"
done
echo MTBGZF_DONE
