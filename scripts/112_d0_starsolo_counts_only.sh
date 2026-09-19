#!/usr/bin/env bash
# D0 STARsolo no-BAM comparator for the post-v1 fair-baseline campaign.
#
# Usage: scripts/112_d0_starsolo_counts_only.sh same-features|gene-only
#
# `same-features` changes only BAM emission relative to the frozen oracle:
# Gene, GeneFull, SJ, and Velocyto remain enabled. `gene-only` additionally
# matches the current Gravlax replay product. Results from these arms must not
# be pooled or described with the same capability label.
set -euo pipefail

arm=${1:-}
case "$arm" in
  same-features) features=(Gene GeneFull SJ Velocyto) ;;
  gene-only) features=(Gene) ;;
  *) echo "usage: $0 same-features|gene-only" >&2; exit 2 ;;
esac

: "${PROJ_ROOT:?source environment/setup.sh or set PROJ_ROOT explicitly}"
star_bin=${STAR_BIN:-STAR}
threads=${THREADS:-24}
out_dir=${STAR_COUNTS_OUT:-$PROJ_ROOT/runs/post-v1/fair-baselines/d0/starsolo-counts-only-$arm}

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 2; }
command -v "$star_bin" >/dev/null || { echo "STAR not found: $star_bin" >&2; exit 2; }
verify_only=${VERIFY_ONLY:-0}
if [[ "$verify_only" == 1 ]]; then
  [[ -d "$out_dir/Solo.out" ]] || { echo "VERIFY_ONLY output is incomplete: $out_dir" >&2; exit 2; }
elif [[ -e "$out_dir" ]]; then
  echo "refusing to overwrite STAR_COUNTS_OUT: $out_dir" >&2
  exit 2
fi

r1=$PROJ_ROOT/data/fastq/pbmc_1k_v3_S1_R1_001.fastq.gz
r2=$PROJ_ROOT/data/fastq/pbmc_1k_v3_S1_R2_001.fastq.gz
genome_dir=$PROJ_ROOT/ref/star-49
whitelist=$PROJ_ROOT/data/3M-february-2018.txt
oracle=$PROJ_ROOT/runs/oracle/full/v49/Solo.out

for input in "$r1" "$r2" "$genome_dir/Genome" "$genome_dir/SA" "$genome_dir/SAindex" \
  "$whitelist"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done

if [[ "$verify_only" != 1 ]]; then
  mkdir -p "$out_dir"
  {
    printf 'arm\t%s\nthreads\t%s\n' "$arm" "$threads"
    "$star_bin" --version
  } >"$out_dir/tool-and-arm.txt"

  /usr/bin/time -v -o "$out_dir/time.txt" "$star_bin" \
    --runThreadN "$threads" \
    --genomeDir "$genome_dir" \
    --readFilesIn "$r2" "$r1" \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$whitelist" \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 \
    --soloFeatures "${features[@]}" \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIdedup 1MM_CR \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloCellFilter EmptyDrops_CR \
    --outSAMtype None \
    --outSAMunmapped Within \
    --outFileNamePrefix "$out_dir/" \
    >"$out_dir/command.stdout.txt" 2>"$out_dir/command.stderr.txt"
fi

while IFS= read -r artifact; do
  cmp "$out_dir/Solo.out/$artifact" "$oracle/$artifact"
done < <(cd "$out_dir/Solo.out" && find . -type f -printf '%P\n' | sort)

if find "$out_dir" -maxdepth 1 -type f \( -name '*.bam' -o -name '*.sam' \) -print -quit \
  | grep -q .; then
  echo "counts-only arm unexpectedly emitted an alignment file" >&2
  exit 1
fi

printf 'D0 STARsolo counts-only %s passed oracle equivalence: %s\n' "$arm" "$out_dir"
