#!/usr/bin/env bash
# Locked packed-EM benchmark. Run once per dataset; raw runs remain excluded from Git.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN}"
: "${DATASET:?set DATASET}"
: "${ARCHIVE:?set ARCHIVE}"
: "${BARCODES:?set BARCODES}"
: "${REFERENCE_LAYER:?set REFERENCE_LAYER to the frozen v1 em.mtx}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
gravlax_repo=${GRAVLAX_REPO:-$PROJ_ROOT/work/gravlax}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
threads=${THREADS:-24}
repetitions=${REPETITIONS:-3}

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be positive" >&2; exit 2; }
[[ "$repetitions" =~ ^[1-9][0-9]*$ ]] || { echo "REPETITIONS must be positive" >&2; exit 2; }
[[ -x "$AIE_BIN" ]] || { echo "AIE_BIN is not executable: $AIE_BIN" >&2; exit 2; }
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
for input in "$ARCHIVE" "$gtf" "$BARCODES" "$REFERENCE_LAYER"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done

mkdir -p "$OUT_ROOT"
export RAYON_NUM_THREADS=$threads

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-31\n'
  printf 'dataset\t%s\n' "$DATASET"
  printf 'threads\t%s\n' "$threads"
  printf 'repetitions\t%s\n' "$repetitions"
  printf 'layout\t64 cell shards; 8-byte support records; u32 labels; u64 CSR offsets\n'
  printf 'mask_selection\tSplitMix64(seed,cell,class)\n'
  printf 'gravlax_commit\t%s\n' "$(git -C "$gravlax_repo" rev-parse HEAD)"
  printf 'gravlax_dirty_entries\t%s\n' "$(git -C "$gravlax_repo" status --porcelain | wc -l)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'archive_sha256\t%s\n' "$(sha256sum "$ARCHIVE" | cut -d' ' -f1)"
  printf 'gtf_sha256\t%s\n' "$(sha256sum "$gtf" | cut -d' ' -f1)"
  printf 'piscem_infer_reference_commit\t7182fce639d3a0994fa99bf621fa76d66e92d286\n'
  printf 'salmon_reference_commit\td3be8b31d31242683c01fc5510c48fbd8bc917fb\n'
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

printf 'rep\tarm\n' >"$OUT_ROOT/schedule.tsv"
for rep in $(seq 1 "$repetitions"); do
  masked=$OUT_ROOT/masked-r$rep
  impact=$OUT_ROOT/impact-r$rep
  mkdir -p "$masked" "$impact/emit"
  printf '%s\tmasked\n' "$rep" >>"$OUT_ROOT/schedule.tsv"
  /usr/bin/time -v -o "$masked/time.txt" \
    "$AIE_BIN" em "$ARCHIVE" --gtf "$gtf" --mask 0.2 --seed 7 --alpha 20 \
    >"$masked/stdout.txt" 2>"$masked/stderr.txt"
  if (( rep > 1 )); then
    cmp "$OUT_ROOT/masked-r1/stdout.txt" "$masked/stdout.txt"
  fi

  printf '%s\timpact\n' "$rep" >>"$OUT_ROOT/schedule.tsv"
  /usr/bin/time -v -o "$impact/time.txt" \
    "$AIE_BIN" em "$ARCHIVE" --gtf "$gtf" --mask 0 --emit "$impact/emit" \
      --barcodes "$BARCODES" >"$impact/stdout.txt" 2>"$impact/stderr.txt"
  "$scripts_repo/scripts/133_compare_em_layers.py" \
    "$impact/emit/em.mtx" "$REFERENCE_LAYER" --out "$impact/layer-comparison.json"
  cmp "$impact/emit/features.tsv" "$(dirname "$REFERENCE_LAYER")/features.tsv"
  cmp "$impact/emit/barcodes.tsv" "$(dirname "$REFERENCE_LAYER")/barcodes.tsv"
done

printf 'packed sharded EM gate complete: %s\n' "$OUT_ROOT"

