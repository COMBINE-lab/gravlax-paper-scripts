#!/usr/bin/env bash
# Export masked candidate-gene unions, then construct replay-only candidate-excluded groups.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN}"
: "${PYTHON_BIN:?set PYTHON_BIN to the frozen Scanpy environment Python}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
manifest=${MANIFEST:-$scripts_repo/experiments/manifests/hierarchical-em.yaml}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
threads=${THREADS:-24}

[[ -x "$AIE_BIN" ]] || { echo "AIE_BIN is not executable: $AIE_BIN" >&2; exit 2; }
[[ -x "$PYTHON_BIN" ]] || { echo "PYTHON_BIN is not executable: $PYTHON_BIN" >&2; exit 2; }
[[ -r "$manifest" ]] || { echo "manifest is not readable: $manifest" >&2; exit 2; }
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
mkdir -p "$OUT_ROOT/candidates" "$OUT_ROOT/logs"
export RAYON_NUM_THREADS=$threads

printf 'dataset\tseed\tarchive\tcandidate_file\n' >"$OUT_ROOT/candidate-runs.tsv"
for item in D0:d0 D1:d1 D2_prime:d2p; do
  dataset=${item%%:*}
  stem=${item##*:}
  archive=$PROJ_ROOT/runs/archive/$stem.aie
  for seed in 7 17 29; do
    output=$OUT_ROOT/candidates/$dataset.seed-$seed.genes.txt
    /usr/bin/time -v -o "$OUT_ROOT/logs/$dataset.seed-$seed.time.txt" \
      "$AIE_BIN" em "$archive" --gtf "$gtf" --mask 0.2 --seed "$seed" \
        --candidate-genes-out "$output" --candidate-genes-only \
        >"$OUT_ROOT/logs/$dataset.seed-$seed.stdout.txt" \
        2>"$OUT_ROOT/logs/$dataset.seed-$seed.stderr.txt"
    printf '%s\t%s\t%s\t%s\n' "$dataset" "$seed" "$archive" "$output" \
      >>"$OUT_ROOT/candidate-runs.tsv"
  done
done

"$PYTHON_BIN" "$scripts_repo/scripts/144_build_hierarchical_em_groups.py" \
  --project-root "$PROJ_ROOT" --manifest "$manifest" \
  --candidate-dir "$OUT_ROOT/candidates" --out-dir "$OUT_ROOT/groups"

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-31\n'
  printf 'threads\t%s\n' "$threads"
  printf 'gravlax_commit\t%s\n' "$(git -C "$PROJ_ROOT/work/gravlax" rev-parse HEAD)"
  printf 'scripts_commit\t%s\n' "$(git -C "$scripts_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'manifest_sha256\t%s\n' "$(sha256sum "$manifest" | cut -d' ' -f1)"
  printf 'group_script_sha256\t%s\n' "$(sha256sum "$scripts_repo/scripts/144_build_hierarchical_em_groups.py" | cut -d' ' -f1)"
  printf 'python\t%s\n' "$PYTHON_BIN"
  printf 'python_version\t%s\n' "$($PYTHON_BIN --version 2>&1)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

printf 'hierarchical EM groups prepared: %s\n' "$OUT_ROOT"
