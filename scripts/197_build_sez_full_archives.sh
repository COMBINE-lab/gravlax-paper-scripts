#!/usr/bin/env bash
# Build one annotation-free, all-published-QC-cell archive per GSE234790 donor.
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 PROJECT_ROOT OUT_ROOT" >&2
  exit 2
fi

GV_ROOT=$(realpath "$1")
GV_OUT=$(realpath -m "$2")
GV_AIE=${AIE:-$GV_ROOT/env/cargo-target-features/release/aie}
GV_STAR=${STAR:-$GV_ROOT/env/sc/bin/STAR}
GV_SAMTOOLS=${SAMTOOLS:-$GV_ROOT/env/sc/bin/samtools}
GV_THREADS=${THREADS:-20}
GV_DONOR_JOBS=${DONOR_JOBS:-4}
GV_FASTQ_ROOT=$GV_ROOT/data/multidonor-fnbp1/gse234790/fastq
GV_GROUP_ROOT=$GV_ROOT/gravlax-paper-scripts/results/post-v1-sez-transcript-end-feasibility/groups
GV_ACQUISITION=$GV_FASTQ_ROOT/acquisition.json
GV_GENOME=$GV_ROOT/data/GRCh38.primary_assembly.genome.fa.gz
GV_STAR_INDEX=$GV_ROOT/ref/star-base
GV_FULL_WHITELIST=$GV_ROOT/data/3M-february-2018.txt
GV_GATE=$GV_ROOT/gravlax-paper-scripts/experiments/manifests/sez-transcript-end-atlas-v1.json

[[ "$GV_DONOR_JOBS" =~ ^[1-8]$ ]] || {
  echo "DONOR_JOBS must be an integer from 1 through 8" >&2
  exit 2
}
mkdir -p "$GV_OUT"

python3 -c '
import hashlib,json,sys
gate=json.load(open(sys.argv[1])); acq=json.load(open(sys.argv[2]))
assert gate["status"] == "FROZEN_BEFORE_WHOLE_GENOME_ENDPOINT_INSPECTION"
assert gate["registration"]["event_blind_feasibility_sha256"] == "ad9b0c1e1d0acb63eff4a4c671abee38b47edd0ad5fbceb5774e9535ffbe351b"
assert acq["status"] == "PASS" and acq["event_evidence_inspected"] is False
assert acq["runs"] == 40 and acq["total_bytes"] == 44692024247
' "$GV_GATE" "$GV_ACQUISITION"

run_donor() {
  local donor=$1
  local donor_out=$GV_OUT/donor-$donor
  local star_out=$donor_out/star
  local fq_root=$GV_FASTQ_ROOT/donor-$donor
  local groups=$GV_GROUP_ROOT/donor-$donor.tsv
  local bam=$donor_out/Aligned.out.bam
  local archive=$donor_out/sez.called.aie
  mkdir -p "$star_out" "$donor_out/tmp"
  cut -f1 "$groups" > "$donor_out/called-whitelist.txt"

  local expected r1_csv r2_csv r1_count r2_count
  expected=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["donors"][sys.argv[2]])' "$GV_ACQUISITION" "$donor")
  r1_csv=$(find "$fq_root" -maxdepth 1 -type f -name '*_1.fastq.gz' | sort | paste -sd,)
  r2_csv=$(find "$fq_root" -maxdepth 1 -type f -name '*_2.fastq.gz' | sort | paste -sd,)
  r1_count=$(tr ',' '\n' <<< "$r1_csv" | wc -l)
  r2_count=$(tr ',' '\n' <<< "$r2_csv" | wc -l)
  [[ "$r1_count" -eq "$expected" && "$r2_count" -eq "$expected" ]] || {
    echo "donor $donor FASTQ count mismatch" >&2
    return 1
  }

  if [[ ! -s "$donor_out/alignment.complete.sha256" ]]; then
    echo "donor $donor: annotation-free alignment" >&2
    TMPDIR=$donor_out/tmp /usr/bin/time -v -o "$donor_out/alignment.time.txt" \
      "$GV_STAR" --runThreadN "$GV_THREADS" --genomeDir "$GV_STAR_INDEX" \
        --genomeLoad NoSharedMemory --twopassMode Basic \
        --readFilesIn "$r2_csv" "$r1_csv" --readFilesCommand zcat \
        --soloType CB_UMI_Simple --soloCBwhitelist "$GV_FULL_WHITELIST" \
        --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
        --soloBarcodeReadLength 0 --soloFeatures SJ \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --outSAMattributes NH HI AS nM CR UR CY UY \
        --outSAMtype BAM Unsorted --outSAMunmapped None \
        --outFileNamePrefix "$star_out/" \
        > "$donor_out/star.stdout.txt" 2> "$donor_out/star.stderr.txt"
    mv "$star_out/Aligned.out.bam" "$bam"
    "$GV_SAMTOOLS" quickcheck -v "$bam"
    sha256sum "$bam" "$star_out/Log.final.out" > "$donor_out/alignment.complete.sha256"
  fi

  if [[ ! -s "$archive" ]]; then
    echo "donor $donor: archive ingest" >&2
    TMPDIR=$donor_out/tmp /usr/bin/time -v -o "$donor_out/ingest.time.txt" \
      "$GV_AIE" ingest-archive "$bam" --whitelist "$donor_out/called-whitelist.txt" \
        --genome "$GV_GENOME" --out "$archive" \
        > "$donor_out/ingest.stdout.txt" 2> "$donor_out/ingest.stderr.txt"
  fi
  sha256sum "$GV_AIE" "$groups" "$donor_out/called-whitelist.txt" "$archive" \
    > "$donor_out/archive.complete.sha256"
  echo "donor $donor: complete" >&2
}

pids=()
donors=()
failed=0
for donor in A B C D E F G H; do
  run_donor "$donor" &
  pids+=("$!")
  donors+=("$donor")
  if [[ ${#pids[@]} -ge $GV_DONOR_JOBS ]]; then
    if ! wait "${pids[0]}"; then
      echo "donor ${donors[0]} failed" >&2
      failed=1
    fi
    pids=("${pids[@]:1}")
    donors=("${donors[@]:1}")
  fi
done
for index in "${!pids[@]}"; do
  if ! wait "${pids[$index]}"; then
    echo "donor ${donors[$index]} failed" >&2
    failed=1
  fi
done
[[ $failed -eq 0 ]] || exit 1

sha256sum "$GV_GATE" "$GV_AIE" > "$GV_OUT/workflow-inputs.sha256"
echo "all eight full-genome archives complete" >&2
