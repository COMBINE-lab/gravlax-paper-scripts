#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
export RAYON_NUM_THREADS=24

SCRIPT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ROOT=$(dirname "$SCRIPT_REPO")
AIE_BIN=${AIE_BIN:-$ROOT/work/gravlax/target/release/aie}
PYTHON=${PYTHON:-python3}
RUN_DIR=$ROOT/runs/post-v1/query-junctions
ARCHIVE=$ROOT/runs/archive/d0.aie
LOCUS=chr11:35138870-35232402
MIN_SUPPORT=20

mkdir -p "$RUN_DIR"

for rep in 1 2 3 4 5; do
  /usr/bin/time -f '%e\t%M' -o "$RUN_DIR/index-$rep.time" \
    "$AIE_BIN" query "$ARCHIVE" junctions "$LOCUS" --min-support "$MIN_SUPPORT" --tsv \
    > "$RUN_DIR/index-$rep.tsv" 2> "$RUN_DIR/index-$rep.stderr"
  /usr/bin/time -f '%e\t%M' -o "$RUN_DIR/cells-$rep.time" \
    "$AIE_BIN" query "$ARCHIVE" junctions "$LOCUS" --min-support "$MIN_SUPPORT" \
      --with-cells --json > "$RUN_DIR/cells-$rep.json" 2> "$RUN_DIR/cells-$rep.stderr"
done

for rep in 2 3 4 5; do
  cmp "$RUN_DIR/index-1.tsv" "$RUN_DIR/index-$rep.tsv"
  cmp "$RUN_DIR/cells-1.json" "$RUN_DIR/cells-$rep.json"
done

"$AIE_BIN" query "$ARCHIVE" junctions "$LOCUS" --min-support "$MIN_SUPPORT" --json \
  > "$RUN_DIR/index.json" 2> "$RUN_DIR/index-json.stderr"
"$AIE_BIN" query "$ARCHIVE" junctions "$LOCUS" --min-support "$MIN_SUPPORT" \
  --gtf "$ROOT/annotations/gencode.v49.annotation.gtf" --tsv \
  > "$RUN_DIR/annotated.tsv" 2> "$RUN_DIR/annotated.stderr"
"$AIE_BIN" query "$ARCHIVE" junctions "$LOCUS" --min-support "$MIN_SUPPORT" \
  --min-cells 20 --tsv > "$RUN_DIR/min-cells.tsv" 2> "$RUN_DIR/min-cells.stderr"
"$AIE_BIN" query "$ARCHIVE" junctions chr11:35189834-35221733 \
  --min-support "$MIN_SUPPORT" --either --tsv \
  > "$RUN_DIR/either.tsv" 2> "$RUN_DIR/either.stderr"
"$AIE_BIN" query "$ARCHIVE" junctions chr11:35189834-35221733 \
  --min-support "$MIN_SUPPORT" --json \
  > "$RUN_DIR/contained-boundary.json" 2> "$RUN_DIR/contained-boundary.stderr"

CODE_COMMIT=$(git -C "$ROOT/work/gravlax" rev-parse HEAD)
"$PYTHON" "$SCRIPT_REPO/scripts/139_summarize_query_junctions.py" \
  --root "$ROOT" --run-dir "$RUN_DIR" --aie-bin "$AIE_BIN" \
  --code-commit "$CODE_COMMIT" \
  --output "$SCRIPT_REPO/results/post-v1-query-junctions.json"

