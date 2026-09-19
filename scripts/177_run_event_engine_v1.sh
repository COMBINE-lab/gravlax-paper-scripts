#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:-${GRAVLAX_PROJECT_ROOT:?set GRAVLAX_PROJECT_ROOT to the analysis project root (see README)}}
BIN=${AIE_BIN:-$ROOT/env/cargo-target-features/release/aie}
OUT=${EVENT_ENGINE_RUN_DIR:-$ROOT/runs/post-v1/event-engine-v1}
REPETITIONS=${REPETITIONS:-7}

mkdir -p "$OUT"

args=(
  --root "$ROOT"
  --bin "$BIN"
  --out "$OUT"
  --repetitions "$REPETITIONS"
)
if [[ -n "${REFERENCE_AIE_BIN:-}" ]]; then
  args+=(--reference-bin "$REFERENCE_AIE_BIN")
fi

exec /tmp/gravlax-downstream-env/bin/python \
  "$ROOT/gravlax-paper-scripts/scripts/178_validate_event_engine_v1.py" \
  "${args[@]}"

