#!/usr/bin/env bash
# Genome-wide, annotation-independent event federation for biological signal screening.
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 ROOT AIE OUT" >&2
  exit 2
fi

ROOT=$(realpath "$1")
AIE=$(realpath "$2")
OUT=$(realpath -m "$3")
mkdir -p "$OUT"

PYTHON=${PYTHON:-python3}
THREADS=${THREADS:-24}
WORKERS=${WORKERS:-4}

if [[ ! "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
  echo "THREADS must be a positive integer (got '$THREADS')" >&2
  exit 2
fi
if [[ ! "$WORKERS" =~ ^[1-9][0-9]*$ ]]; then
  echo "WORKERS must be a positive integer (got '$WORKERS')" >&2
  exit 2
fi

# Keep concurrent Rayon pools within the requested aggregate thread budget.
EFFECTIVE_WORKERS=$WORKERS
if (( EFFECTIVE_WORKERS > THREADS )); then
  EFFECTIVE_WORKERS=$THREADS
fi
PER_PROCESS_THREADS=$((THREADS / EFFECTIVE_WORKERS))
DRIVER_STARTED_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)
DRIVER_START_NS=$(date +%s%N)

"$PYTHON" "$ROOT/gravlax-paper-scripts/scripts/179_build_federated_biology_groups.py" \
  --root "$ROOT" --out "$OUT"

AIC="$OUT/gencode.v49.aic"
if [[ ! -s "$AIC" ]]; then
  "$AIE" compile-annotation "$ROOT/annotations/gencode.v49.annotation.gtf" --out "$AIC"
fi
AIE_SHA=$(sha256sum "$AIE" | awk '{print $1}')
AIC_SHA=$(sha256sum "$AIC" | awk '{print $1}')
PBMC_SCOPE_SHA=$(sha256sum "$OUT/D0.TM.tsv" "$OUT/D1.TM.tsv" \
  "$OUT/D3.TM.tsv" "$OUT/D4.TM.tsv" | awk '{print $1}' | sha256sum | awk '{print $1}')
SIX_SCOPE_SHA=$(sha256sum "$OUT/D0.called.tsv" "$OUT/D1.called.tsv" \
  "$OUT/D2.called.tsv" "$OUT/D2p.called.tsv" "$OUT/D3.called.tsv" \
  "$OUT/D4.called.tsv" | awk '{print $1}' | sha256sum | awk '{print $1}')

run_timed_atomic() {
  local output=$1
  local timing=$2
  local stderr=$3
  local signature=$4
  shift 4

  # Use same-directory temporary files so the final rename is atomic. The JSON
  # moves last and is the completion marker used by resumability checks.
  local token="${BASHPID:-$$}"
  local output_tmp="${output}.tmp.${token}"
  local timing_tmp="${timing}.tmp.${token}"
  local stderr_tmp="${stderr}.tmp.${token}"
  local complete="${output}.complete"
  local complete_tmp="${complete}.tmp.${token}"
  trap 'rm -f "$output_tmp" "$timing_tmp" "$stderr_tmp" "$complete_tmp"' RETURN

  if ! RAYON_NUM_THREADS="$PER_PROCESS_THREADS" \
    /usr/bin/time -f 'wall_seconds=%e\tpeak_rss_kib=%M' -o "$timing_tmp" \
      "$@" > "$output_tmp" 2> "$stderr_tmp"; then
    echo "federated query failed for $output" >&2
    return 1
  fi
  if [[ ! -s "$output_tmp" || ! -s "$timing_tmp" ]]; then
    echo "federated query produced an incomplete result for $output" >&2
    return 1
  fi

  mv -f "$stderr_tmp" "$stderr"
  mv -f "$timing_tmp" "$timing"
  mv -f "$output_tmp" "$output"
  printf '%s\n' "$signature" > "$complete_tmp"
  mv -f "$complete_tmp" "$complete"
  trap - RETURN
}

result_complete() {
  local output=$1
  local timing=$2
  local stderr=$3
  local signature=$4
  local complete="${output}.complete"
  local actual_signature=""
  [[ -s "$output" && -s "$timing" && -e "$stderr" && -s "$complete" ]] || return 1
  IFS= read -r actual_signature < "$complete"
  [[ "$actual_signature" == "$signature" ]]
}

run_chromosome() {
  local CHROM=$1
  local LENGTH=$2
  local PBMC="$OUT/pbmc4-tm.$CHROM.json"
  local SIX="$OUT/six-called.$CHROM.json"
  local PBMC_TIME="$OUT/pbmc4-tm.$CHROM.time.txt"
  local PBMC_STDERR="$OUT/pbmc4-tm.$CHROM.stderr.txt"
  local SIX_TIME="$OUT/six-called.$CHROM.time.txt"
  local SIX_STDERR="$OUT/six-called.$CHROM.stderr.txt"
  local PBMC_SIGNATURE="driver=v3|arm=pbmc4-tm|locus=$CHROM:0-$LENGTH|event-types=all|min-support=2|min-samples=3|min-informative=30|min-row=10|max-events=2000000|threads=$PER_PROCESS_THREADS|workers=$EFFECTIVE_WORKERS|aie=$AIE_SHA|aic=$AIC_SHA|scope=$PBMC_SCOPE_SHA"
  local SIX_SIGNATURE="driver=v3|arm=six-called|locus=$CHROM:0-$LENGTH|event-types=all|min-support=2|min-samples=3|min-informative=50|min-row=20|max-events=2000000|threads=$PER_PROCESS_THREADS|workers=$EFFECTIVE_WORKERS|aie=$AIE_SHA|aic=$AIC_SHA|scope=$SIX_SCOPE_SHA"
  if ! result_complete "$PBMC" "$PBMC_TIME" "$PBMC_STDERR" "$PBMC_SIGNATURE"; then
    run_timed_atomic "$PBMC" "$PBMC_TIME" "$PBMC_STDERR" "$PBMC_SIGNATURE" \
      "$AIE" cohort events "$CHROM:0-$LENGTH" \
      --sample "D0=$ROOT/runs/archive/d0.aie" \
      --sample "D1=$ROOT/runs/archive/d1.aie" \
      --sample "D3=$ROOT/runs/archive/d3.aie" \
      --sample "D4=$ROOT/runs/archive/d4.aie" \
      --groups "D0=$OUT/D0.TM.tsv" --groups "D1=$OUT/D1.TM.tsv" \
      --groups "D3=$OUT/D3.TM.tsv" --groups "D4=$OUT/D4.TM.tsv" \
      --min-support 2 --min-samples 3 --min-informative 30 --max-events 2000000 \
      --min-row-informative 10 --gtf "$AIC" --json
  fi
  if ! result_complete "$SIX" "$SIX_TIME" "$SIX_STDERR" "$SIX_SIGNATURE"; then
    run_timed_atomic "$SIX" "$SIX_TIME" "$SIX_STDERR" "$SIX_SIGNATURE" \
      "$AIE" cohort events "$CHROM:0-$LENGTH" \
      --sample "D0=$ROOT/runs/archive/d0.aie" \
      --sample "D1=$ROOT/runs/archive/d1.aie" \
      --sample "D2=$ROOT/runs/archive/d2.aie" \
      --sample "D2p=$ROOT/runs/archive/d2p.aie" \
      --sample "D3=$ROOT/runs/archive/d3.aie" \
      --sample "D4=$ROOT/runs/archive/d4.aie" \
      --groups "D0=$OUT/D0.called.tsv" --groups "D1=$OUT/D1.called.tsv" \
      --groups "D2=$OUT/D2.called.tsv" --groups "D2p=$OUT/D2p.called.tsv" \
      --groups "D3=$OUT/D3.called.tsv" --groups "D4=$OUT/D4.called.tsv" \
      --min-support 2 --min-samples 3 --min-informative 50 --max-events 2000000 \
      --min-row-informative 20 --gtf "$AIC" --json
  fi
}

mapfile -t CHROMOSOMES < <(
  awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print $1, $2}' \
    "$ROOT/ref/GRCh38.primary_assembly.genome.fa.fai"
)
if [[ ${#CHROMOSOMES[@]} -ne 24 ]]; then
  echo "expected 24 primary chromosomes, found ${#CHROMOSOMES[@]}" >&2
  exit 1
fi

RUNNING=0
FAILED=0
for CHROMOSOME in "${CHROMOSOMES[@]}"; do
  read -r CHROM LENGTH <<< "$CHROMOSOME"
  run_chromosome "$CHROM" "$LENGTH" &
  RUNNING=$((RUNNING + 1))
  if (( RUNNING >= EFFECTIVE_WORKERS )); then
    if ! wait -n; then
      FAILED=1
    fi
    RUNNING=$((RUNNING - 1))
  fi
done
while (( RUNNING > 0 )); do
  if ! wait -n; then
    FAILED=1
  fi
  RUNNING=$((RUNNING - 1))
done
if (( FAILED != 0 )); then
  echo "one or more federated chromosome jobs failed; completed results remain resumable" >&2
  exit 1
fi

"$PYTHON" "$ROOT/gravlax-paper-scripts/scripts/181_summarize_federated_biology.py" \
  --root "$ROOT" --screen "$OUT" \
  --out-json "$OUT/summary.json" \
  --out-pbmc "$OUT/pbmc4-consistent.tsv" \
  --out-tissue "$OUT/brain-associated.tsv" \
  --driver-started-utc "$DRIVER_STARTED_UTC" \
  --driver-start-ns "$DRIVER_START_NS" \
  --driver-workers "$EFFECTIVE_WORKERS" \
  --driver-total-threads "$THREADS" \
  --driver-per-process-threads "$PER_PROCESS_THREADS"
