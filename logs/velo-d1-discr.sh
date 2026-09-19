set -x
P=${GRAVLAX_PROJECT_ROOT}
AIE=$P/src/target/release/aie
PY=${GRAVLAX_PYTHON}
# D1 velocity gate
$AIE replay-rows $P/runs/archive/d1.aie --gtf $P/annotations/gencode.v49.annotation.gtf \
  --barcodes $P/runs/oracle/d1/v49/Solo.out/Velocyto/raw/barcodes.tsv \
  --out-dir $P/runs/replay/d1-velocity --velocity
$PY $P/scripts/91_velo_gate.py $P/runs/oracle/d1/v49/Solo.out/Velocyto/raw \
  $P/runs/replay/d1-velocity $P/runs/oracle/d1/v49/Solo.out/Gene/filtered/barcodes.tsv
echo D1_VELO_DONE
# DISCR loop: discover on w1, emit GTF, replay under it
$AIE query $P/runs/archive/d0.aie discover --gtf $P/annotations/withheld/w1.gtf --tsv \
  --emit-gtf $P/runs/replay/d0-novel.gtf > $P/results/discover-w1-v14.tsv
$AIE replay-rows $P/runs/archive/d0.aie --gtf $P/runs/replay/d0-novel.gtf \
  --barcodes $P/runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv \
  --out-dir $P/runs/replay/d0-novel
echo DISCR_REPLAY_DONE
