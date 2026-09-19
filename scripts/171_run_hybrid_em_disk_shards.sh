#!/usr/bin/env bash
# Reproduce the locked exact disk-backed support-shard memory gate.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN to Gravlax commit 93def86f}"
: "${PREPARED:?set PREPARED to the frozen D3 preparation from script 160}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

archive=${ARCHIVE:-$PROJ_ROOT/runs/archive/d3.aie}
gtf=${GTF:-$PROJ_ROOT/annotations/gencode.v49.annotation.gtf}
groups=${GROUPS:-$PREPARED/groups/D3.shuffled.tsv}
reference=${REFERENCE:-$PROJ_ROOT/runs/post-v1/hybrid-em-promotion-d3-r1/shuffled/seed-7}
threads=${THREADS:-24}

for input in "$AIE_BIN" "$archive" "$gtf" "$groups" \
  "$reference/metrics.json" "$reference/stdout.txt"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done
[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
mkdir -p "$OUT_ROOT/tmp" "$OUT_ROOT/seed-7"
export RAYON_NUM_THREADS=$threads
export TMPDIR=$OUT_ROOT/tmp

/usr/bin/time -v -o "$OUT_ROOT/seed-7/time.txt" \
  "$AIE_BIN" em "$archive" --gtf "$gtf" --mask 0.2 --seed 7 --alpha 20 \
    --groups "$groups" \
    --convex-cell-weight 0.20 --convex-group-weight 0.60 --convex-group-prior 80 \
    --dirichlet-cell-prior 64 --dirichlet-group-prior 80 \
    --hybrid-depth-scale 8 --hybrid-depth-power 8 --hybrid-only \
    --metrics-json "$OUT_ROOT/seed-7/metrics.json" \
    >"$OUT_ROOT/seed-7/stdout.txt" 2>"$OUT_ROOT/seed-7/stderr.txt"

cmp "$reference/metrics.json" "$OUT_ROOT/seed-7/metrics.json"
cmp "$reference/stdout.txt" "$OUT_ROOT/seed-7/stdout.txt"
if find "$TMPDIR" -mindepth 1 -name 'gravlax-em-*' -print -quit | grep -q .; then
  echo "EM temporary support artifacts remain after exit" >&2
  exit 1
fi

rss=$(awk -F: '/Maximum resident set size/{gsub(/[[:space:]]/, "", $2); print $2}' \
  "$OUT_ROOT/seed-7/time.txt")
wall=$(awk -F': ' '/Elapsed \(wall clock\) time/{print $NF}' "$OUT_ROOT/seed-7/time.txt")
wall_seconds=$(python3 -c '
import sys
p = [float(x) for x in sys.argv[1].split(":")]
print(sum(x * 60**i for i, x in enumerate(reversed(p))))
' "$wall")
(( rss <= 4194304 )) || { echo "peak RSS exceeds 4 GiB: $rss KiB" >&2; exit 1; }
python3 -c 'import sys; raise SystemExit(float(sys.argv[1]) > 60)' "$wall_seconds" || {
  echo "wall time exceeds 60 seconds: $wall_seconds" >&2
  exit 1
}
{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-31\n'
  printf 'scope\tlocked_D3_disk_support_shard_gate\n'
  printf 'threads\t%s\n' "$threads"
  printf 'peak_rss_kib\t%s\n' "$rss"
  printf 'wall_seconds\t%s\n' "$wall_seconds"
  printf 'metrics_sha256\t%s\n' "$(sha256sum "$OUT_ROOT/seed-7/metrics.json" | cut -d' ' -f1)"
  printf 'stdout_sha256\t%s\n' "$(sha256sum "$OUT_ROOT/seed-7/stdout.txt" | cut -d' ' -f1)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'archive_sha256\t%s\n' "$(sha256sum "$archive" | cut -d' ' -f1)"
  printf 'gtf_sha256\t%s\n' "$(sha256sum "$gtf" | cut -d' ' -f1)"
  printf 'groups_sha256\t%s\n' "$(sha256sum "$groups" | cut -d' ' -f1)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"
printf 'disk-backed EM support-shard gate passed: %s KiB RSS\n' "$rss"
