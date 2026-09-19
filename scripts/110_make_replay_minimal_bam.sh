#!/usr/bin/env bash
# Construct the ingest-equivalent BAM used by the post-v1 fair-baseline pilot.
#
# Gravlax's current BAM extraction consumes reference/start, flags, CIGAR,
# QNAME (to group NH>1 alignments), and the CR/CY/UR/NH tags.  This transform
# deliberately retains those fields and removes information not consumed by
# that path: MAPQ, read sequence, read quality, and all other auxiliary tags.
# Acceptance still requires byte-identical Gravlax archives and matrices; this
# field audit alone is not treated as proof of equivalence.
#
# Usage:
#   scripts/110_make_replay_minimal_bam.sh INPUT.bam OUTPUT.bam [THREADS]
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "usage: $0 INPUT.bam OUTPUT.bam [THREADS]" >&2
  exit 2
fi

input=$1
output=$2
threads=${3:-${AIE_THREADS:-24}}
samtools_bin=${SAMTOOLS:-samtools}

[[ -r "$input" ]] || { echo "input is not readable: $input" >&2; exit 2; }
[[ "$input" != "$output" ]] || { echo "input and output must differ" >&2; exit 2; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 2; }
command -v "$samtools_bin" >/dev/null || { echo "samtools not found: $samtools_bin" >&2; exit 2; }
command -v awk >/dev/null || { echo "awk not found" >&2; exit 2; }

mkdir -p "$(dirname "$output")"
partial="${output}.partial.$$"
trap 'rm -f "$partial"' EXIT

"$samtools_bin" quickcheck -v "$input"
"$samtools_bin" view -@ "$threads" -h --keep-tag CR,CY,UR,NH "$input" \
  | awk 'BEGIN { OFS="\t" }
         /^@/ { print; next }
         NF < 11 { print "malformed SAM record with fewer than 11 fields" > "/dev/stderr"; exit 3 }
         { $5=0; $10="*"; $11="*"; print }' \
  | "$samtools_bin" view -@ "$threads" -b -o "$partial" -

"$samtools_bin" quickcheck -v "$partial"
mv "$partial" "$output"
trap - EXIT

printf 'input\t%s\t%s\n' "$input" "$(stat -c%s "$input")"
printf 'replay_minimal_bam\t%s\t%s\n' "$output" "$(stat -c%s "$output")"
