#!/usr/bin/env bash
# Run the remotely registered GSE234790 fixed-event validation with donor-isolated workers.
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 ROOT OUT" >&2
  exit 2
fi

GV_ROOT=$(realpath "$1")
GV_OUT=$(realpath -m "$2")
GV_AIE=${AIE:-$GV_ROOT/env/cargo-target-features/release/aie}
GV_STAR=${STAR:-$GV_ROOT/env/sc/bin/STAR}
GV_SAMTOOLS=${SAMTOOLS:-$GV_ROOT/env/sc/bin/samtools}
GV_PYTHON=${PYTHON:-python3}
GV_THREADS=${THREADS:-20}
GV_SAM_THREADS=${SAM_THREADS:-4}
GV_DONOR_JOBS=${DONOR_JOBS:-8}
[[ "$GV_DONOR_JOBS" =~ ^[1-8]$ ]] || {
  echo "DONOR_JOBS must be an integer from 1 through 8" >&2
  exit 2
}
GV_MANIFEST=$GV_ROOT/gravlax-paper-scripts/experiments/manifests/fnbp1-multidonor-sez-validation.json
GV_FEASIBILITY=$GV_ROOT/gravlax-paper-scripts/results/post-v1-fnbp1-multidonor-sez-feasibility/feasibility.json
GV_GROUP_ROOT=$GV_ROOT/gravlax-paper-scripts/results/post-v1-fnbp1-multidonor-sez-feasibility/groups
GV_FASTQ_ROOT=$GV_ROOT/data/multidonor-fnbp1/gse234790/fastq
GV_ACQUISITION=$GV_FASTQ_ROOT/acquisition.json
GV_ENA=$GV_ROOT/data/multidonor-fnbp1/gse234790/PRJNA983239-ena-runs.tsv
GV_BED=$GV_ROOT/gravlax-paper-scripts/experiments/manifests/fnbp1-window.grch38.bed
GV_FULL_WHITELIST=$GV_ROOT/data/3M-february-2018.txt
GV_GENOME=$GV_ROOT/data/GRCh38.primary_assembly.genome.fa.gz
GV_STAR_INDEX=$GV_ROOT/ref/star-base
GV_AIC=$GV_OUT/gencode.v49.aic
mkdir -p "$GV_OUT"

"$GV_PYTHON" -c '
import hashlib,json,sys
manifest=json.load(open(sys.argv[1])); acquisition=json.load(open(sys.argv[2]))
assert manifest["status"] == "PROSPECTIVE_GATE_FROZEN_BEFORE_ANY_FNBP1_SEQUENCE_OR_EVENT_QUERY"
assert acquisition["status"] == "PASS" and acquisition["event_evidence_inspected"] is False
assert acquisition["manifest_sha256"] == hashlib.sha256(open(sys.argv[1],"rb").read()).hexdigest()
assert acquisition["runs"] == 40 and acquisition["total_bytes"] == 44692024247
' "$GV_MANIFEST" "$GV_ACQUISITION"

if [[ ! -s "$GV_AIC" ]]; then
  "$GV_AIE" compile-annotation "$GV_ROOT/annotations/gencode.v49.annotation.gtf" --out "$GV_AIC"
fi

run_donor() {
  GV_DONOR=$1
  echo "donor $GV_DONOR: start" >&2
  GV_DONOR_OUT=$GV_OUT/donor-$GV_DONOR
  GV_STAR_OUT=$GV_DONOR_OUT/star
  GV_FASTQ=$GV_FASTQ_ROOT/donor-$GV_DONOR
  GV_GROUP_FULL=$GV_GROUP_ROOT/donor-$GV_DONOR.tsv
  mkdir -p "$GV_STAR_OUT"

  GV_EXPECTED_RUNS=$("$GV_PYTHON" -c '
import json,sys
d=json.load(open(sys.argv[1])); print(d["donors"][sys.argv[2]])
' "$GV_ACQUISITION" "$GV_DONOR")
  GV_R1_CSV=$(find "$GV_FASTQ" -maxdepth 1 -type f -name '*_1.fastq.gz' | sort | paste -sd,)
  GV_R2_CSV=$(find "$GV_FASTQ" -maxdepth 1 -type f -name '*_2.fastq.gz' | sort | paste -sd,)
  GV_R1_COUNT=$(tr ',' '\n' <<< "$GV_R1_CSV" | wc -l)
  GV_R2_COUNT=$(tr ',' '\n' <<< "$GV_R2_CSV" | wc -l)
  [[ "$GV_R1_COUNT" -eq "$GV_EXPECTED_RUNS" && "$GV_R2_COUNT" -eq "$GV_EXPECTED_RUNS" ]] || {
    echo "donor $GV_DONOR FASTQ count mismatch: expected $GV_EXPECTED_RUNS runs, got $GV_R1_COUNT/$GV_R2_COUNT" >&2
    exit 1
  }
  if find "$GV_FASTQ" -maxdepth 1 -type f -name '*.aria2' | grep -q .; then
    echo "donor $GV_DONOR retains an incomplete aria2 transfer" >&2
    exit 1
  fi

  GV_TARGET_BAM=$GV_DONOR_OUT/fnbp1.target.bam
  if [[ ! -s "$GV_DONOR_OUT/alignment.complete.sha256" ]]; then
    export GV_STAR GV_SAMTOOLS GV_THREADS GV_SAM_THREADS GV_STAR_INDEX GV_R1_CSV GV_R2_CSV
    export GV_FULL_WHITELIST GV_STAR_OUT GV_TARGET_BAM GV_BED
    /usr/bin/time -v -o "$GV_DONOR_OUT/alignment.time.txt" /usr/bin/bash -c '
set -euo pipefail
"$GV_STAR" --runThreadN "$GV_THREADS" --genomeDir "$GV_STAR_INDEX" \
  --genomeLoad NoSharedMemory --twopassMode Basic \
  --readFilesIn "$GV_R2_CSV" "$GV_R1_CSV" --readFilesCommand zcat \
  --soloType CB_UMI_Simple --soloCBwhitelist "$GV_FULL_WHITELIST" \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 --soloFeatures SJ \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --outSAMattributes NH HI AS nM CR UR CY UY \
  --outSAMtype BAM Unsorted --outStd BAM_Unsorted --outSAMunmapped None \
  --outFileNamePrefix "$GV_STAR_OUT/" 2> "$GV_STAR_OUT.stderr.txt" \
  | "$GV_SAMTOOLS" view -@ "$GV_SAM_THREADS" -b -L "$GV_BED" -o "$GV_TARGET_BAM" -
'
    "$GV_SAMTOOLS" quickcheck -v "$GV_TARGET_BAM"
    sha256sum "$GV_TARGET_BAM" "$GV_STAR_OUT/Log.final.out" \
      "$GV_STAR_OUT/Solo.out/SJ/raw/features.tsv" \
      "$GV_STAR_OUT/Solo.out/SJ/raw/barcodes.tsv" \
      "$GV_STAR_OUT/Solo.out/SJ/raw/matrix.mtx" \
      > "$GV_DONOR_OUT/alignment.complete.sha256"
  fi

  cut -f1 "$GV_GROUP_FULL" > "$GV_DONOR_OUT/called-whitelist.txt"
  if [[ ! -s "$GV_DONOR_OUT/fnbp1.called.aie" ]]; then
    /usr/bin/time -v -o "$GV_DONOR_OUT/ingest.time.txt" \
      "$GV_AIE" ingest-archive "$GV_TARGET_BAM" \
        --whitelist "$GV_DONOR_OUT/called-whitelist.txt" --genome "$GV_GENOME" \
        --out "$GV_DONOR_OUT/fnbp1.called.aie"
  fi
  "$GV_AIE" query "$GV_DONOR_OUT/fnbp1.called.aie" region chr9:0-138394717 \
    --top 100000 --json > "$GV_DONOR_OUT/archive-cells.json"
  "$GV_PYTHON" "$GV_ROOT/gravlax-paper-scripts/scripts/188_filter_fnbp1_archive_groups.py" \
    --full-map "$GV_GROUP_FULL" --archive-cells "$GV_DONOR_OUT/archive-cells.json" \
    --out "$GV_DONOR_OUT/groups.archive.tsv" > "$GV_DONOR_OUT/scope.stdout.json"
  /usr/bin/time -v -o "$GV_DONOR_OUT/query.time.txt" \
    "$GV_AIE" query "$GV_DONOR_OUT/fnbp1.called.aie" jset \
      --include chr9:129908999-129923843 --include chr9:129924026-129924959 \
      --exclude chr9:129908999-129924959 \
      --groups "$GV_DONOR_OUT/groups.archive.tsv" --agg group --json \
      > "$GV_DONOR_OUT/jset-by-group.json"
  sha256sum "$GV_AIE" "$GV_AIC" "$GV_GROUP_FULL" "$GV_DONOR_OUT/groups.archive.tsv" \
    "$GV_DONOR_OUT/fnbp1.called.aie" "$GV_DONOR_OUT/jset-by-group.json" \
    > "$GV_DONOR_OUT/execution-inputs.sha256"
  echo "donor $GV_DONOR: complete" >&2
}

GV_PIDS=()
GV_PID_DONORS=()
GV_FAILED=0
for GV_DONOR in A B C D E F G H; do
  run_donor "$GV_DONOR" &
  GV_PIDS+=("$!")
  GV_PID_DONORS+=("$GV_DONOR")
  if [[ "${#GV_PIDS[@]}" -ge "$GV_DONOR_JOBS" ]]; then
    if ! wait "${GV_PIDS[0]}"; then
      echo "donor ${GV_PID_DONORS[0]} worker failed" >&2
      GV_FAILED=1
    fi
    GV_PIDS=("${GV_PIDS[@]:1}")
    GV_PID_DONORS=("${GV_PID_DONORS[@]:1}")
  fi
done
for GV_INDEX in "${!GV_PIDS[@]}"; do
  if ! wait "${GV_PIDS[$GV_INDEX]}"; then
    echo "donor ${GV_PID_DONORS[$GV_INDEX]} worker failed" >&2
    GV_FAILED=1
  fi
done
[[ "$GV_FAILED" -eq 0 ]] || exit 1

"$GV_PYTHON" "$GV_ROOT/gravlax-paper-scripts/scripts/193_summarize_fnbp1_multidonor_validation.py" \
  --manifest "$GV_MANIFEST" --feasibility "$GV_FEASIBILITY" --acquisition "$GV_ACQUISITION" \
  --run "$GV_OUT" --executable "$GV_AIE" \
  --out-json "$GV_OUT/summary.json" --out-tsv "$GV_OUT/donor-group-counts.tsv"
