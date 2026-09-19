#!/usr/bin/env bash
# Run the locked oracle-only multiresolution selector. This stage must not open replay matrices.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${PYTHON_BIN:?set PYTHON_BIN to the frozen Scanpy environment Python}"
: "${OUT_DIR:?set OUT_DIR to a new directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
manifest=${MANIFEST:-$scripts_repo/experiments/manifests/downstream-multiresolution.yaml}
threads=${THREADS:-24}
run_date=${RUN_DATE:-$(date +%F)}

[[ -x "$PYTHON_BIN" ]] || { echo "PYTHON_BIN is not executable: $PYTHON_BIN" >&2; exit 2; }
[[ -r "$manifest" ]] || { echo "manifest is not readable: $manifest" >&2; exit 2; }
[[ ! -e "$OUT_DIR" ]] || { echo "refusing to overwrite OUT_DIR: $OUT_DIR" >&2; exit 2; }
mkdir -p "$OUT_DIR/matplotlib" "$OUT_DIR/numba-cache"

export OMP_NUM_THREADS=$threads
export OPENBLAS_NUM_THREADS=$threads
export MKL_NUM_THREADS=$threads
export NUMBA_NUM_THREADS=$threads
export VECLIB_MAXIMUM_THREADS=$threads
export NUMEXPR_NUM_THREADS=$threads
export MPLCONFIGDIR=$OUT_DIR/matplotlib
export NUMBA_CACHE_DIR=$OUT_DIR/numba-cache

{
  printf 'key\tvalue\n'
  printf 'stage\toracle_only_selection\n'
  printf 'date\t%s\n' "$run_date"
  printf 'threads\t%s\n' "$threads"
  printf 'python\t%s\n' "$PYTHON_BIN"
  printf 'python_version\t%s\n' "$($PYTHON_BIN --version 2>&1)"
  printf 'manifest\t%s\n' "$manifest"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$manifest" | cut -d' ' -f1)"
  printf 'script_sha256\t%s\n' "$(sha256sum "$scripts_repo/scripts/129_multiresolution_select.py" | cut -d' ' -f1)"
  printf 'scripts_commit\t%s\n' "$(git -C "$scripts_repo" rev-parse HEAD)"
  printf 'scripts_dirty_entries\t%s\n' "$(git -C "$scripts_repo" status --porcelain | wc -l)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_DIR/protocol.tsv"

/usr/bin/time -v -o "$OUT_DIR/time.txt" \
  "$PYTHON_BIN" "$scripts_repo/scripts/129_multiresolution_select.py" \
    --project-root "$PROJ_ROOT" --manifest "$manifest" \
    --out "$OUT_DIR/oracle-selection.json" \
    >"$OUT_DIR/stdout.txt" 2>"$OUT_DIR/stderr.txt"

"$PYTHON_BIN" -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d["status"] == "complete"; print("\n".join("{}\t{}".format(x["id"], x["locked_selection"]) for x in d["datasets"]))' \
  "$OUT_DIR/oracle-selection.json" | tee "$OUT_DIR/summary.txt"

printf 'oracle-only multiresolution selection complete: %s\n' "$OUT_DIR"
