#!/usr/bin/env bash
# Reproduce the five-warm-run D0 Gene-only counts comparison reported post-v1.
#
# This intentionally reproduces the executed pilot order (five STARsolo runs,
# then five Gravlax runs). It is not the later randomized warm/cold protocol.
set -euo pipefail

: "${PROJ_ROOT:?source environment/setup.sh or set PROJ_ROOT explicitly}"
aie_bin=${AIE_BIN:-$PROJ_ROOT/src/target/release/aie}
star_bin=${STAR_BIN:-STAR}
threads=${THREADS:-24}
bench_root=${BENCH_ROOT:-$PROJ_ROOT/runs/post-v1/fair-baselines/d0/warm-comparison-reproduction}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

[[ -x "$aie_bin" ]] || { echo "AIE_BIN is not executable: $aie_bin" >&2; exit 2; }
command -v "$star_bin" >/dev/null || { echo "STAR not found: $star_bin" >&2; exit 2; }
[[ ! -e "$bench_root" ]] || { echo "refusing to overwrite BENCH_ROOT: $bench_root" >&2; exit 2; }
mkdir -p "$bench_root"

for rep in 1 2 3 4 5; do
  PROJ_ROOT="$PROJ_ROOT" STAR_BIN="$star_bin" THREADS="$threads" \
    STAR_COUNTS_OUT="$bench_root/star-r$rep" \
    "$script_dir/112_d0_starsolo_counts_only.sh" gene-only
done

archive=$PROJ_ROOT/runs/archive/d0.aie
gtf=$PROJ_ROOT/annotations/gencode.v49.annotation.gtf
barcodes=$PROJ_ROOT/runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv
reference=$PROJ_ROOT/runs/replay/full-frombam
export RAYON_NUM_THREADS=$threads

for rep in 1 2 3 4 5; do
  out=$bench_root/aie-r$rep
  mkdir -p "$out"
  /usr/bin/time -v -o "$out/time.txt" "$aie_bin" replay-rows "$archive" \
    --gtf "$gtf" --barcodes "$barcodes" --out-dir "$out" \
    >"$out/command.stdout.txt" 2>"$out/command.stderr.txt"
  for artifact in matrix.mtx barcodes.tsv features.tsv; do
    cmp "$out/$artifact" "$reference/$artifact"
  done
done

printf 'D0 five-run warm comparison passed: %s\n' "$bench_root"

