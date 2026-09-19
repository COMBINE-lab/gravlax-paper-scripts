#!/usr/bin/env bash
# D1 geometry-fidelity experiment: build the D1 archive with --geometry-fidelity, then compare
# archive size, Gene replay deviation, and velocity-component deviation against the compact archive
# and the STARsolo v49 reference. Both archives are replayed with the same 0.2.2 binary.
# Archived from the analysis host as scripts/151_d1_geometry_fidelity.sh (2026-09-11); renumbered here to follow the repository sequence.
set -euo pipefail
P=${GRAVLAX_PROJECT_ROOT:?set GRAVLAX_PROJECT_ROOT to the analysis project root (see README)}
AIE=$P/runs/gravlax-release-0.2.2-20260910/aie-0.2.2
PY=${GRAVLAX_PYTHON:-python3}
OUT=$P/runs/paper-followup-20260911/d1-geometry-fidelity
BAM=$P/runs/ingest/d1/Aligned.sortedByCoord.out.bam
WL=$P/data/3M-february-2018.txt
GTF=$P/annotations/gencode.v49.annotation.gtf
COMPACT=$P/runs/archive/d1.aie
ORACLE=$P/runs/oracle/d1/v49/Solo.out
BARC=$ORACLE/Gene/raw/barcodes.tsv
CELLS=$ORACLE/Gene/filtered/barcodes.tsv
READS=383900000
THREADS=24; CPUSET=24-47
[[ ! -e $OUT ]] || { echo "refusing to overwrite $OUT" >&2; exit 2; }
mkdir -p $OUT/logs; export RAYON_NUM_THREADS=$THREADS
FID=$OUT/d1-fidelity.aie
echo "[$(date)] ingest with --geometry-fidelity"
/usr/bin/time -v -o $OUT/ingest.time.txt taskset -c $CPUSET $AIE ingest-archive $BAM --whitelist $WL \
  --out $FID --zstd-level 19 --geometry-fidelity --report-format json --report-output $OUT/ingest-report.json \
  > $OUT/logs/ingest.stdout.txt 2> $OUT/logs/ingest.stderr.txt
for kind in compact fidelity; do
  arc=$COMPACT; [[ $kind == fidelity ]] && arc=$FID
  echo "[$(date)] replay $kind: Gene"
  /usr/bin/time -v -o $OUT/replay-$kind-gene.time.txt taskset -c $CPUSET $AIE replay-rows $arc --gtf $GTF --barcodes $BARC \
    --out-dir $OUT/replay-$kind-gene > $OUT/logs/replay-$kind-gene.txt 2>&1
  echo "[$(date)] replay $kind: velocity"
  /usr/bin/time -v -o $OUT/replay-$kind-velocity.time.txt taskset -c $CPUSET $AIE replay-rows $arc --gtf $GTF --barcodes $BARC \
    --out-dir $OUT/replay-$kind-velocity --velocity > $OUT/logs/replay-$kind-velocity.txt 2>&1
done
# regression: compact replay under 0.2.2 must reproduce the archived D1 results
st=exact
for a in matrix.mtx barcodes.tsv features.tsv; do cmp -s $OUT/replay-compact-gene/$a $P/runs/replay/d1-aie/$a || st=DIFFERS; done
echo "compact_gene_vs_archived_replay $st" | tee $OUT/regression.txt
st=exact
for a in spliced.mtx unspliced.mtx ambiguous.mtx; do cmp -s $OUT/replay-compact-velocity/$a $P/runs/replay/d1-velocity/$a || st=DIFFERS; done
echo "compact_velocity_vs_archived_replay $st" | tee -a $OUT/regression.txt
{
  echo "== sizes"
  for f in $COMPACT $FID; do b=$(stat -c %s $f); printf '%s\t%d bytes\t%.2f bits/read\n' "$(basename $f)" $b "$(echo "$b*8/$READS" | bc -l)"; done
  for kind in compact fidelity; do
    echo "== Gene deviation ($kind) vs STARsolo v49, filtered cells"
    $PY $P/scripts/30_gateb_compare.py $ORACLE/Gene/raw $OUT/replay-$kind-gene oracle replay --cells $CELLS
    echo "== velocity components ($kind)"
    $PY $P/scripts/91_velo_gate.py $ORACLE/Velocyto/raw $OUT/replay-$kind-velocity $CELLS
    echo "== velocity completeness strata ($kind)"
    $PY $P/scripts/94_velo2_strata.py $ORACLE/Velocyto/raw $OUT/replay-$kind-velocity $CELLS
  done
} > $OUT/summary.txt 2>&1
cat $OUT/summary.txt
echo "D1 GEOMETRY FIDELITY DONE"
