#!/usr/bin/env bash
# Reproduce the locked compiled-annotation and batched-query gates.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT to the annotation-independent-evidence project root}"
: "${AIE_BIN:?set AIE_BIN to Gravlax commit 4bde878b640fed18460a1defaddc9d84e53153b4}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
archive=${ARCHIVE:-$PROJ_ROOT/runs/archive/d0.aie}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
barcodes=${BARCODES:-$PROJ_ROOT/runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv}
plan=${PLAN:-$script_dir/../experiments/plans/batched-query-d0.tsv}
threads=${THREADS:-24}
commit=${AIE_COMMIT:-4bde878b640fed18460a1defaddc9d84e53153b4}
[[ "$threads" == 24 ]] || {
  echo "the locked gate requires THREADS=24 (received $threads)" >&2
  exit 2
}

for input in "$AIE_BIN" "$archive" "$gtf" "$barcodes" "$plan"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
mkdir -p "$OUT_ROOT"/{annotation,replay-gtf,replay-aic,discovery,annotation-query,batch,independent}
export RAYON_NUM_THREADS=$threads

aic=$OUT_ROOT/annotation/gencode.v49.annotation.aic
aic_repeat=$OUT_ROOT/annotation/gencode.v49.annotation.repeat.aic
/usr/bin/time -f '%e\t%M' -o "$OUT_ROOT/annotation/compile-time.tsv" \
  "$AIE_BIN" compile-annotation "$gtf" --out "$aic" \
  >"$OUT_ROOT/annotation/compile.stdout" 2>"$OUT_ROOT/annotation/compile.stderr"
"$AIE_BIN" compile-annotation "$gtf" --out "$aic_repeat" \
  >"$OUT_ROOT/annotation/compile-repeat.stdout" \
  2>"$OUT_ROOT/annotation/compile-repeat.stderr"
cmp "$aic" "$aic_repeat"

/usr/bin/time -f '%e\t%M' -o "$OUT_ROOT/replay-gtf/time.tsv" \
  "$AIE_BIN" replay-rows "$archive" --gtf "$gtf" --barcodes "$barcodes" \
  --out-dir "$OUT_ROOT/replay-gtf/output" \
  >"$OUT_ROOT/replay-gtf/stdout.txt" 2>"$OUT_ROOT/replay-gtf/stderr.txt"
/usr/bin/time -f '%e\t%M' -o "$OUT_ROOT/replay-aic/time.tsv" \
  "$AIE_BIN" replay-rows "$archive" --gtf "$aic" --barcodes "$barcodes" \
  --out-dir "$OUT_ROOT/replay-aic/output" \
  >"$OUT_ROOT/replay-aic/stdout.txt" 2>"$OUT_ROOT/replay-aic/stderr.txt"
for file in matrix.mtx features.tsv barcodes.tsv; do
  cmp "$OUT_ROOT/replay-gtf/output/$file" "$OUT_ROOT/replay-aic/output/$file"
done

"$AIE_BIN" query "$archive" discover --gtf "$gtf" --tsv \
  >"$OUT_ROOT/discovery/gtf.tsv" 2>"$OUT_ROOT/discovery/gtf.stderr"
"$AIE_BIN" query "$archive" discover --gtf "$aic" --tsv \
  >"$OUT_ROOT/discovery/aic.tsv" 2>"$OUT_ROOT/discovery/aic.stderr"
cmp "$OUT_ROOT/discovery/gtf.tsv" "$OUT_ROOT/discovery/aic.tsv"

: >"$OUT_ROOT/annotation-query/gtf-times.tsv"
: >"$OUT_ROOT/annotation-query/aic-times.tsv"
for run in 1 2 3 4 5; do
  /usr/bin/time -f "$run\t%e\t%M" -o "$OUT_ROOT/annotation-query/gtf-times.tsv" -a \
    "$AIE_BIN" query "$archive" junctions chr11:35138870-35232402 \
    --min-support 20 --gtf "$gtf" --json \
    >"$OUT_ROOT/annotation-query/gtf-$run.json" \
    2>"$OUT_ROOT/annotation-query/gtf-$run.stderr"
  /usr/bin/time -f "$run\t%e\t%M" -o "$OUT_ROOT/annotation-query/aic-times.tsv" -a \
    "$AIE_BIN" query "$archive" junctions chr11:35138870-35232402 \
    --min-support 20 --gtf "$aic" --json \
    >"$OUT_ROOT/annotation-query/aic-$run.json" \
    2>"$OUT_ROOT/annotation-query/aic-$run.stderr"
done
cmp "$OUT_ROOT/annotation-query/gtf-1.json" "$OUT_ROOT/annotation-query/aic-1.json"
for run in 2 3 4 5; do
  cmp "$OUT_ROOT/annotation-query/aic-1.json" "$OUT_ROOT/annotation-query/aic-$run.json"
done

: >"$OUT_ROOT/batch/times.tsv"
: >"$OUT_ROOT/independent/times.tsv"
for run in 1 2 3 4 5 6 7; do
  /usr/bin/time -f "$run\t%e\t%M" -o "$OUT_ROOT/batch/times.tsv" -a \
    "$AIE_BIN" query "$archive" batch --plan "$plan" --top 20 \
    >"$OUT_ROOT/batch/run-$run.json" 2>"$OUT_ROOT/batch/run-$run.stderr"
  AIE_BIN="$AIE_BIN" AIE_ARCHIVE="$archive" AIE_PLAN="$plan" \
    /usr/bin/time -f "$run\t%e\t%M" -o "$OUT_ROOT/independent/times.tsv" -a \
    /usr/bin/bash -c '
      set -euo pipefail
      tail -n +2 "$AIE_PLAN" | while IFS=$'"'"'\t'"'"' read -r id kind locus; do
        if [[ "$kind" == region ]]; then
          "$AIE_BIN" query "$AIE_ARCHIVE" region "$locus" --top 20 >/dev/null 2>/dev/null
        else
          "$AIE_BIN" query "$AIE_ARCHIVE" junction "$locus" --top 20 --json >/dev/null 2>/dev/null
        fi
      done
    '
done
for run in 2 3 4 5 6 7; do
  cmp "$OUT_ROOT/batch/run-1.json" "$OUT_ROOT/batch/run-$run.json"
done

python3 "$script_dir/173_validate_batched_query.py" \
  --aie "$AIE_BIN" --archive "$archive" --plan "$plan" \
  --batch "$OUT_ROOT/batch/run-1.json" --top 20 \
  --out "$OUT_ROOT/batch/correctness.json"

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-31\n'
  printf 'gravlax_commit\t%s\n' "$commit"
  printf 'threads\t%s\n' "$threads"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'archive_sha256\t%s\n' "$(sha256sum "$archive" | cut -d' ' -f1)"
  printf 'gtf_sha256\t%s\n' "$(sha256sum "$gtf" | cut -d' ' -f1)"
  printf 'aic_sha256\t%s\n' "$(sha256sum "$aic" | cut -d' ' -f1)"
  printf 'barcodes_sha256\t%s\n' "$(sha256sum "$barcodes" | cut -d' ' -f1)"
  printf 'plan_sha256\t%s\n' "$(sha256sum "$plan" | cut -d' ' -f1)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

python3 "$script_dir/174_summarize_compiled_annotations_batched_query.py" \
  "$OUT_ROOT" --gtf "$gtf" --commit "$commit" --out "$OUT_ROOT/summary.json"
printf 'compiled-annotation and batched-query gates passed: %s\n' "$OUT_ROOT/summary.json"
