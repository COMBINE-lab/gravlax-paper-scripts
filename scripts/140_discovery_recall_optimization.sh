#!/usr/bin/env bash
set -euo pipefail

SCRIPT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ROOT=$(dirname "$SCRIPT_REPO")
AIE_BIN=${AIE_BIN:-$ROOT/work/gravlax/target/release/aie}
PYTHON=${PYTHON:-python3}
CPUSET=${CPUSET:-0-23}
RUN_DIR=${RUN_DIR:-$ROOT/runs/post-v1/discovery-recall-optimization}
RESULT_BASE=$SCRIPT_REPO/results/post-v1-discovery-recall-optimization

mkdir -p "$RUN_DIR"
CODE_COMMIT=$(git -C "$ROOT/work/gravlax" rev-parse HEAD)

run_arm() {
  local label=$1
  local arm=$2
  shift 2
  /usr/bin/time -v -o "$RUN_DIR/$label.$arm.time" \
    taskset -c "$CPUSET" "$AIE_BIN" query "$ROOT/runs/archive/$arm.aie" discover \
      --gtf "$ROOT/annotations/gencode.v32.annotation.gtf" "$@" --tsv \
      > "$RUN_DIR/$label.$arm.tsv" 2> "$RUN_DIR/$label.$arm.stderr"
}

evaluate() {
  local label=$1
  "$PYTHON" "$SCRIPT_REPO/scripts/137_complete_denominator_discovery.py" \
    --root "$ROOT" \
    --d0-discovery "$RUN_DIR/$label.d0.tsv" \
    --d1-discovery "$RUN_DIR/$label.d1.tsv" \
    --aie-bin "$AIE_BIN" \
    --code-commit "$CODE_COMMIT" \
    --summary "$RUN_DIR/$label.summary.json" \
    --truth-tsv "$RUN_DIR/$label.truth.tsv" \
    --candidate-tsv "$RUN_DIR/$label.candidates.tsv"
}

for arm in d0 d1; do
  run_arm span "$arm"
  run_arm strand-span "$arm" --claim-mode strand-span
  run_arm compatible "$arm" --claim-mode compatible
  for support in 25 50 75 100; do
    run_arm "residual-sites-$support" "$arm" \
      --claim-mode residual-sites --residual-min-umis "$support"
  done
done

for label in span strand-span compatible \
  residual-sites-25 residual-sites-50 residual-sites-75 residual-sites-100; do
  evaluate "$label"
done

"$PYTHON" "$SCRIPT_REPO/scripts/141_summarize_discovery_recall_optimization.py" \
  --run-dir "$RUN_DIR" --output "$RESULT_BASE.json"

cp "$RUN_DIR/residual-sites-75.summary.json" "$RESULT_BASE-selected.json"
cp "$RUN_DIR/residual-sites-75.truth.tsv" "$RESULT_BASE-truth.tsv"
cp "$RUN_DIR/residual-sites-75.candidates.tsv" "$RESULT_BASE-candidates.tsv"

