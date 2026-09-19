#!/usr/bin/env bash
set -euo pipefail

SCRIPT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ROOT=$(dirname "$SCRIPT_REPO")
AIE_BIN=${AIE_BIN:-$ROOT/work/gravlax/target/release/aie}
BASELINE_AIE_BIN=${BASELINE_AIE_BIN:?set BASELINE_AIE_BIN to a release build of 0efabed0d3e92a72c432b4c815c2e0b70c7155e4}
CPUSET=${CPUSET:-0-23}
RUN_DIR=${RUN_DIR:-$ROOT/runs/post-v1/hotpath-optimization}
PYTHON=${PYTHON:-python3}

mkdir -p "$RUN_DIR"

run_timed() {
  local time_file=$1
  local stdout_file=$2
  local stderr_file=$3
  shift 3
  /usr/bin/time -v -o "$time_file" taskset -c "$CPUSET" "$@" \
    > "$stdout_file" 2> "$stderr_file"
}

for build in baseline candidate; do
  if [[ $build == baseline ]]; then bin=$BASELINE_AIE_BIN; else bin=$AIE_BIN; fi
  for arm in d0 d1; do
    for rep in 1 2 3; do
      run_timed "$RUN_DIR/$build.discovery.$arm.$rep.time" \
        "$RUN_DIR/$build.discovery.$arm.$rep.tsv" \
        "$RUN_DIR/$build.discovery.$arm.$rep.stderr" \
        "$bin" query "$ROOT/runs/archive/$arm.aie" discover \
          --gtf "$ROOT/annotations/gencode.v32.annotation.gtf" --tsv

      if [[ $arm == d0 ]]; then
        barcodes=$ROOT/runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv
      else
        barcodes=$ROOT/runs/oracle/d1/v49/Solo.out/Gene/raw/barcodes.tsv
      fi
      out_dir=$RUN_DIR/$build.replay.$arm.$rep
      mkdir -p "$out_dir"
      run_timed "$RUN_DIR/$build.replay.$arm.$rep.time" \
        "$RUN_DIR/$build.replay.$arm.$rep.stdout" \
        "$RUN_DIR/$build.replay.$arm.$rep.stderr" \
        "$bin" replay-rows "$ROOT/runs/archive/$arm.aie" \
          --gtf "$ROOT/annotations/gencode.v49.annotation.gtf" \
          --barcodes "$barcodes" --out-dir "$out_dir"
    done
  done

  for rep in 1 2 3 4 5; do
    run_timed "$RUN_DIR/$build.junctions.d1.$rep.time" \
      "$RUN_DIR/$build.junctions.d1.$rep.tsv" \
      "$RUN_DIR/$build.junctions.d1.$rep.stderr" \
      "$bin" query "$ROOT/runs/archive/d1.aie" junctions chr1:0-249250000 \
        --min-support 20 --with-cells --tsv
  done
done

for arm in d0 d1; do
  for rep in 1 2 3; do
    cmp "$RUN_DIR/baseline.discovery.$arm.$rep.tsv" \
      "$RUN_DIR/candidate.discovery.$arm.$rep.tsv"
    for file in features.tsv barcodes.tsv matrix.mtx; do
      cmp "$RUN_DIR/baseline.replay.$arm.$rep/$file" \
        "$RUN_DIR/candidate.replay.$arm.$rep/$file"
    done
  done
done
for rep in 1 2 3 4 5; do
  cmp "$RUN_DIR/baseline.junctions.d1.$rep.tsv" \
    "$RUN_DIR/candidate.junctions.d1.$rep.tsv"
done

"$PYTHON" "$SCRIPT_REPO/scripts/143_summarize_hotpath_optimization.py" \
  --run-dir "$RUN_DIR" \
  --baseline-commit 0efabed0d3e92a72c432b4c815c2e0b70c7155e4 \
  --candidate-commit "$(git -C "$ROOT/work/gravlax" rev-parse HEAD)" \
  --output "$SCRIPT_REPO/results/post-v1-hotpath-optimization.json"

