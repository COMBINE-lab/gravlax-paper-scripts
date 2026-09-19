#!/usr/bin/env bash
# Controlled D0 Gene-only timing comparison.
#
# One best-effort specific-file cold observation per arm is followed by five warm replicates in a
# balanced two-arm crossover schedule. The first arm was selected from preregistered seed 20260830;
# the order then alternates so each arm runs first in at least two blocks. Host-wide cache drops are
# deliberately forbidden on the shared server.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT to the project root}"
: "${AIE_BIN:?set AIE_BIN to the Gravlax binary under test}"
star_bin=${STAR_BIN:-STAR}
threads=${THREADS:-24}
bench_root=${BENCH_ROOT:-$PROJ_ROOT/runs/post-v1/fair-baselines/d0/randomized-counts-v1}
gravlax_repo=${GRAVLAX_REPO:-$PROJ_ROOT/work/gravlax}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 2; }
[[ -x "$AIE_BIN" ]] || { echo "AIE_BIN is not executable: $AIE_BIN" >&2; exit 2; }
command -v "$star_bin" >/dev/null || { echo "STAR not found: $star_bin" >&2; exit 2; }
git -C "$gravlax_repo" rev-parse --verify HEAD >/dev/null 2>&1 || {
  echo "GRAVLAX_REPO is not a Git checkout: $gravlax_repo" >&2
  exit 2
}
[[ ! -e "$bench_root" ]] || { echo "refusing to overwrite BENCH_ROOT: $bench_root" >&2; exit 2; }

archive=$PROJ_ROOT/runs/archive/d0.aie
gtf=$PROJ_ROOT/annotations/gencode.v49.annotation.gtf
barcodes=$PROJ_ROOT/runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv
reference=$PROJ_ROOT/runs/replay/full-frombam
r1=$PROJ_ROOT/data/fastq/pbmc_1k_v3_S1_R1_001.fastq.gz
r2=$PROJ_ROOT/data/fastq/pbmc_1k_v3_S1_R2_001.fastq.gz
genome_dir=$PROJ_ROOT/ref/star-49
whitelist=$PROJ_ROOT/data/3M-february-2018.txt

for input in "$archive" "$gtf" "$barcodes" "$reference/matrix.mtx" \
  "$reference/features.tsv" "$reference/barcodes.tsv" "$r1" "$r2" \
  "$genome_dir/Genome" "$genome_dir/SA" "$genome_dir/SAindex" "$whitelist"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done

mkdir -p "$bench_root"
export RAYON_NUM_THREADS=$threads

run_aie() {
  local label=$1
  local out=$bench_root/$label
  mkdir -p "$out"
  /usr/bin/time -v -o "$out/time.txt" "$AIE_BIN" replay-rows "$archive" \
    --gtf "$gtf" --barcodes "$barcodes" --out-dir "$out" \
    >"$out/command.stdout.txt" 2>"$out/command.stderr.txt"
  for artifact in matrix.mtx barcodes.tsv features.tsv; do
    cmp "$out/$artifact" "$reference/$artifact"
  done
}

run_star() {
  local label=$1
  PROJ_ROOT="$PROJ_ROOT" STAR_BIN="$star_bin" THREADS="$threads" \
    STAR_COUNTS_OUT="$bench_root/$label" \
    "$script_dir/112_d0_starsolo_counts_only.sh" gene-only
}

{
  printf 'key\tvalue\n'
  printf 'date\t2026-08-30\n'
  printf 'dataset\tD0 / 10x PBMC 1k v3\n'
  printf 'threads\t%s\n' "$threads"
  printf 'schedule_seed\t20260830\n'
  printf 'schedule_design\tbalanced two-arm alternating crossover; first arm selected by seed\n'
  printf 'gravlax_commit\t%s\n' "$(git -C "$gravlax_repo" rev-parse HEAD)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'star_version\t%s\n' "$("$star_bin" --version)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
  printf 'cold_definition\tPOSIX_FADV_DONTNEED on named arm inputs; advisory and non-destructive\n'
} >"$bench_root/protocol.tsv"

printf 'phase\tblock\tposition\tarm\tlabel\n' >"$bench_root/schedule.tsv"

# Cold observations use disjoint primary inputs; evict each arm's named inputs immediately before
# it runs. They are diagnostics because advisory eviction cannot stop another user recaching data.
python3 "$script_dir/115_fadvise_dontneed.py" --log "$bench_root/cold-aie-cache.tsv" \
  "$archive" "$gtf" "$barcodes"
printf 'cold\t0\t1\taie\tcold-aie\n' >>"$bench_root/schedule.tsv"
run_aie cold-aie

python3 "$script_dir/115_fadvise_dontneed.py" --log "$bench_root/cold-star-cache.tsv" \
  "$r1" "$r2" "$genome_dir" "$whitelist"
printf 'cold\t0\t2\tstar\tcold-star\n' >>"$bench_root/schedule.tsv"
run_star cold-star

# Seed 20260830 selected STAR first; alternating the order is the balanced Latin-square design for
# two treatments. With five blocks, STAR is first three times and AIE twice.
for block in 1 2 3 4 5; do
  if (( block % 2 == 1 )); then
    arms=(star aie)
  else
    arms=(aie star)
  fi
  position=0
  for arm in "${arms[@]}"; do
    position=$((position + 1))
    label=warm-b${block}-p${position}-${arm}
    printf 'warm\t%s\t%s\t%s\t%s\n' "$block" "$position" "$arm" "$label" \
      >>"$bench_root/schedule.tsv"
    if [[ "$arm" == aie ]]; then
      run_aie "$label"
    else
      run_star "$label"
    fi
  done
done

printf 'D0 randomized warm/cold comparison passed: %s\n' "$bench_root"
