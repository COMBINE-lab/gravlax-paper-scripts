#!/usr/bin/env bash
# Track B: the cross-sample junction catalogue — the annotation-free fix for the G-C alignment
# ceiling, in the spirit of intropolis/Snaptron (docs/prior-art.md §2).
#
# Pools splice junctions DISCOVERED FROM DATA by annotation-free STAR runs, and feeds the pooled
# list back into the ingest alignment via --sjdbFileChrStartEnd, which STAR accepts without any
# GTF — so the ingest remains strictly annotation-free. Junctions from sjdb-aware runs are never
# pooled: their SJ.out.tab contains annotated junctions and would leak the GTF into the catalogue.
#
# Usage: scripts/60_junction_catalogue.sh <out_prefix> <SJ.out.tab> [more SJ.out.tab ...]
set -euo pipefail
OUT="${1:?usage: $0 <out_prefix> <SJ.out.tab> [...]}"
shift
: "${PROJ_ROOT:?source env/setup.sh first}"
cd "$PROJ_ROOT"

# SJ.out.tab columns: chrom, intron start (1-based), intron end, strand (0/1/2), motif,
# annotated-in-sjdb flag, unique reads, multi reads, max overhang.
# Pool by junction, summing unique-read support across runs; keep junctions with pooled unique
# support >= 3. The threshold is deliberately applied AFTER pooling — a junction weakly supported
# in one run but seen across runs/samples is exactly the rescue case the catalogue exists for.
MIN_UNIQ=3

awk -v OFS='\t' '
  { key = $1 OFS $2 OFS $3 OFS $4; uniq[key] += $7 }
  END {
    for (k in uniq) if (uniq[k] >= '"$MIN_UNIQ"') print k, uniq[k]
  }
' "$@" | sort -k1,1 -k2,2n > "${OUT}.pooled.tsv"

# STAR --sjdbFileChrStartEnd wants: chrom, intron start, intron end, strand as +/-/.
awk -v OFS='\t' '{ s = ($4==1 ? "+" : ($4==2 ? "-" : ".")); print $1, $2, $3, s }' \
  "${OUT}.pooled.tsv" > "${OUT}.sjdb"

echo "catalogue: $(wc -l < "${OUT}.sjdb") junctions (pooled unique support >= ${MIN_UNIQ}) -> ${OUT}.sjdb"
