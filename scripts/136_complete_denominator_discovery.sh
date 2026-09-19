#!/usr/bin/env bash
set -euo pipefail

SCRIPT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ROOT=$(dirname "$SCRIPT_REPO")
AIE_BIN=${AIE_BIN:-$ROOT/work/gravlax/target/release/aie}
PYTHON=${PYTHON:-python3}
RUN_DIR=$ROOT/runs/post-v1/complete-denominator-discovery
RESULT_BASE=$SCRIPT_REPO/results/post-v1-complete-denominator
TMP_DIR=$(mktemp -d /tmp/gravlax-discovery-determinism.XXXXXX)
trap 'rm -rf "$TMP_DIR"' EXIT

mkdir -p "$RUN_DIR"

"$AIE_BIN" query "$ROOT/runs/archive/d0.aie" discover \
  --gtf "$ROOT/annotations/gencode.v32.annotation.gtf" --tsv \
  > "$RUN_DIR/discover-v32-d0.tsv"
"$AIE_BIN" query "$ROOT/runs/archive/d1.aie" discover \
  --gtf "$ROOT/annotations/gencode.v32.annotation.gtf" --tsv \
  > "$RUN_DIR/discover-v32-d1.tsv"

CODE_COMMIT=$(git -C "$ROOT/work/gravlax" rev-parse HEAD)
COMMON=(
  --root "$ROOT"
  --d0-discovery "$RUN_DIR/discover-v32-d0.tsv"
  --d1-discovery "$RUN_DIR/discover-v32-d1.tsv"
  --aie-bin "$AIE_BIN"
  --code-commit "$CODE_COMMIT"
)

"$PYTHON" "$SCRIPT_REPO/scripts/137_complete_denominator_discovery.py" "${COMMON[@]}" \
  --summary "$RESULT_BASE-discovery.json" \
  --truth-tsv "$RESULT_BASE-truth.tsv" \
  --candidate-tsv "$RESULT_BASE-candidates.tsv"

"$PYTHON" "$SCRIPT_REPO/scripts/137_complete_denominator_discovery.py" "${COMMON[@]}" \
  --summary "$TMP_DIR/summary.json" \
  --truth-tsv "$TMP_DIR/truth.tsv" \
  --candidate-tsv "$TMP_DIR/candidates.tsv"

cmp "$RESULT_BASE-discovery.json" "$TMP_DIR/summary.json"
cmp "$RESULT_BASE-truth.tsv" "$TMP_DIR/truth.tsv"
cmp "$RESULT_BASE-candidates.tsv" "$TMP_DIR/candidates.tsv"
echo "complete-denominator discovery outputs are deterministic"

