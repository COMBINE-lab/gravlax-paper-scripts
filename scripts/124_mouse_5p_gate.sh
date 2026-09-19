#!/usr/bin/env bash
# End-to-end mouse 5' v2 gate: fresh STARsolo, annotation-free ingest, AIE, exact source replay,
# fresh-pipeline deviation, velocity, and wrong-strand negative control.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"
: "${AIE_BIN:?set AIE_BIN built from the recorded Gravlax commit}"
: "${OUT_ROOT:?set OUT_ROOT to a new directory}"

scripts_repo=${SCRIPTS_REPO:-$PROJ_ROOT/gravlax-paper-scripts}
gravlax_repo=${GRAVLAX_REPO:-$PROJ_ROOT/work/gravlax}
star=${STAR_BIN:-$PROJ_ROOT/env/sc/bin/STAR-avx2}
samtools=${SAMTOOLS:-$PROJ_ROOT/env/sc/bin/samtools}
python=${PYTHON_BIN:-python3}
threads=${THREADS:-24}
cpuset=${CPUSET:-0-23}
nofile=${NOFILE_SOFT:-65535}
run_date=${RUN_DATE:-$(date +%F)}
summary_out=${SUMMARY_OUT:-$OUT_ROOT/result.json}

data_dir=$PROJ_ROOT/data/mouse5p-1k
fastq_dir=$data_dir/sc5p_v2_mm_c57bl6_splenocyte_1k_5gex_fastqs
r1=$fastq_dir/sc5p_v2_mm_c57bl6_splenocyte_1k_5gex_S1_L003_R1_001.fastq.gz
r2=$fastq_dir/sc5p_v2_mm_c57bl6_splenocyte_1k_5gex_S1_L003_R2_001.fastq.gz
whitelist=$data_dir/737K-august-2016.txt
gtf=$PROJ_ROOT/ref/mouse-m39/gencode.vM39.primary_assembly.annotation.gtf
fasta=$PROJ_ROOT/ref/mouse-m39/GRCm39.primary_assembly.genome.fa
index=${STAR_INDEX:-$PROJ_ROOT/ref/star-m39-base}

[[ ! -e "$OUT_ROOT" ]] || { echo "refusing to overwrite OUT_ROOT: $OUT_ROOT" >&2; exit 2; }
for input in "$AIE_BIN" "$star" "$samtools" "$r1" "$r2" "$whitelist" "$gtf" "$fasta"; do
  [[ -r "$input" ]] || { echo "required input is not readable: $input" >&2; exit 2; }
done
[[ -x "$AIE_BIN" && -x "$star" && -x "$samtools" ]] || {
  echo "AIE_BIN, STAR_BIN, and SAMTOOLS must be executable" >&2
  exit 2
}
command -v "$python" >/dev/null || { echo "Python not found: $python" >&2; exit 2; }

mkdir -p "$OUT_ROOT"/{index,oracle,ingest,archive,replay-archive,replay-bam,replay-forward,replay-velocity,replay-velocity-bam,comparisons}
export RAYON_NUM_THREADS=$threads

if [[ ! -s "$index/SA" || ! -s "$index/Genome" || ! -s "$index/SAindex" ]]; then
  mkdir -p "$index"
  /usr/bin/time -v -o "$OUT_ROOT/index/time.txt" \
    taskset -c "$cpuset" "$star" --runMode genomeGenerate --runThreadN "$threads" \
      --genomeDir "$index" --genomeFastaFiles "$fasta" --limitGenomeGenerateRAM 80000000000 \
      >"$OUT_ROOT/index/command.stdout.txt" 2>"$OUT_ROOT/index/command.stderr.txt"
fi

solo_common=(
  --runThreadN "$threads" --genomeDir "$index"
  --readFilesIn "$r2" "$r1" --readFilesCommand zcat
  --soloType CB_UMI_Simple --soloCBwhitelist "$whitelist"
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10
  --soloBarcodeReadLength 0 --soloStrand Reverse
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
  --soloUMIdedup 1MM_CR --soloUMIfiltering MultiGeneUMI_CR
  --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30
  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 32000000000
  --outSAMunmapped Within
)

/usr/bin/time -v -o "$OUT_ROOT/oracle/time.txt" \
  prlimit --nofile="$nofile":262144 taskset -c "$cpuset" "$star" \
    "${solo_common[@]}" --sjdbGTFfile "$gtf" --sjdbOverhang 89 \
    --soloFeatures Gene GeneFull SJ Velocyto \
    --outSAMattributes NH HI AS nM CR UR CB UB GX GN \
    --outFileNamePrefix "$OUT_ROOT/oracle/" \
    >"$OUT_ROOT/oracle/command.stdout.txt" 2>"$OUT_ROOT/oracle/command.stderr.txt"

/usr/bin/time -v -o "$OUT_ROOT/ingest/time.txt" \
  prlimit --nofile="$nofile":262144 taskset -c "$cpuset" "$star" \
    "${solo_common[@]}" --twopassMode Basic --soloFeatures SJ \
    --outSAMattributes NH HI AS nM CR UR CY UY \
    --outFileNamePrefix "$OUT_ROOT/ingest/" \
    >"$OUT_ROOT/ingest/command.stdout.txt" 2>"$OUT_ROOT/ingest/command.stderr.txt"

/usr/bin/time -v -o "$OUT_ROOT/archive/time.txt" \
  taskset -c "$cpuset" "$AIE_BIN" ingest-archive \
    "$OUT_ROOT/ingest/Aligned.sortedByCoord.out.bam" --whitelist "$whitelist" \
    --out "$OUT_ROOT/archive/mouse5p-1k.aie" --zstd-level 19 --genome "$fasta" \
    >"$OUT_ROOT/archive/command.stdout.txt" 2>"$OUT_ROOT/archive/command.stderr.txt"
taskset -c "$cpuset" "$AIE_BIN" debug "$OUT_ROOT/archive/mouse5p-1k.aie" \
  >"$OUT_ROOT/archive/debug.txt" 2>"$OUT_ROOT/archive/debug.stderr.txt"

barcodes=$OUT_ROOT/oracle/Solo.out/Gene/raw/barcodes.tsv
for source in archive bam; do
  extra=()
  input=$OUT_ROOT/archive/mouse5p-1k.aie
  if [[ "$source" == bam ]]; then
    input=$OUT_ROOT/ingest/Aligned.sortedByCoord.out.bam
    extra=(--from-bam --whitelist "$whitelist")
  fi
  /usr/bin/time -v -o "$OUT_ROOT/replay-$source/time.txt" \
    taskset -c "$cpuset" "$AIE_BIN" replay-rows "$input" "${extra[@]}" \
      --gtf "$gtf" --barcodes "$barcodes" --out-dir "$OUT_ROOT/replay-$source" \
      --solo-strand reverse \
      >"$OUT_ROOT/replay-$source/command.stdout.txt" \
      2>"$OUT_ROOT/replay-$source/command.stderr.txt"
done

for artifact in matrix.mtx features.tsv barcodes.tsv; do
  cmp "$OUT_ROOT/replay-archive/$artifact" "$OUT_ROOT/replay-bam/$artifact"
done

/usr/bin/time -v -o "$OUT_ROOT/replay-forward/time.txt" \
  taskset -c "$cpuset" "$AIE_BIN" replay-rows "$OUT_ROOT/archive/mouse5p-1k.aie" \
    --gtf "$gtf" --barcodes "$barcodes" --out-dir "$OUT_ROOT/replay-forward" \
    >"$OUT_ROOT/replay-forward/command.stdout.txt" \
    2>"$OUT_ROOT/replay-forward/command.stderr.txt"

/usr/bin/time -v -o "$OUT_ROOT/replay-velocity/time.txt" \
  taskset -c "$cpuset" "$AIE_BIN" replay-rows "$OUT_ROOT/archive/mouse5p-1k.aie" \
    --gtf "$gtf" --barcodes "$barcodes" --out-dir "$OUT_ROOT/replay-velocity" \
    --solo-strand reverse --velocity \
    >"$OUT_ROOT/replay-velocity/command.stdout.txt" \
    2>"$OUT_ROOT/replay-velocity/command.stderr.txt"

/usr/bin/time -v -o "$OUT_ROOT/replay-velocity-bam/time.txt" \
  taskset -c "$cpuset" "$AIE_BIN" replay-rows \
    "$OUT_ROOT/ingest/Aligned.sortedByCoord.out.bam" --from-bam --whitelist "$whitelist" \
    --gtf "$gtf" --barcodes "$barcodes" --out-dir "$OUT_ROOT/replay-velocity-bam" \
    --solo-strand reverse --velocity \
    >"$OUT_ROOT/replay-velocity-bam/command.stdout.txt" \
    2>"$OUT_ROOT/replay-velocity-bam/command.stderr.txt"
for artifact in spliced.mtx unspliced.mtx ambiguous.mtx features.tsv barcodes.tsv entries.complete.tsv; do
  cmp "$OUT_ROOT/replay-velocity/$artifact" "$OUT_ROOT/replay-velocity-bam/$artifact"
done

"$python" "$scripts_repo/scripts/30_gateb_compare.py" \
  "$OUT_ROOT/oracle/Solo.out/Gene/raw" "$OUT_ROOT/replay-archive" fresh_STARsolo Gravlax_replay \
  --cells "$OUT_ROOT/oracle/Solo.out/Gene/filtered/barcodes.tsv" \
  --out "$OUT_ROOT/comparisons/gene-fidelity.json" \
  >"$OUT_ROOT/comparisons/gene-fidelity.stdout.txt"
"$python" "$scripts_repo/scripts/91_velo_gate.py" \
  "$OUT_ROOT/oracle/Solo.out/Velocyto/raw" "$OUT_ROOT/replay-velocity" \
  "$OUT_ROOT/oracle/Solo.out/Gene/filtered/barcodes.tsv" \
  >"$OUT_ROOT/comparisons/velocity-fidelity.txt"

if [[ ! -s "$fasta.fai" ]]; then
  "$samtools" faidx "$fasta"
fi
PROJ_ROOT="$PROJ_ROOT" AIE_BIN="$AIE_BIN" \
  ARCHIVE="$OUT_ROOT/archive/mouse5p-1k.aie" DATASET=mouse5p-1k \
  BARCODES="$barcodes" REFERENCE_MATRIX="$OUT_ROOT/replay-archive" \
  OUT_DIR="$OUT_ROOT/matched" GENOME="$fasta" GTF="$gtf" SOLO_STRAND=reverse \
  SAMTOOLS="$samtools" THREADS="$threads" GRAVLAX_REPO="$gravlax_repo" \
  "$scripts_repo/scripts/118_post_correction_storage.sh"
PROJ_ROOT="$PROJ_ROOT" AIE_BIN="$AIE_BIN" \
  MOLECULE_CRAM="$OUT_ROOT/matched/post-correction.molecules.cram" DATASET=mouse5p-1k \
  BARCODES="$barcodes" REFERENCE_MATRIX="$OUT_ROOT/replay-archive" \
  OUT_ROOT="$OUT_ROOT/matched-cram-replay" GENOME="$fasta" GTF="$gtf" \
  SOLO_STRAND=reverse SAMTOOLS="$samtools" THREADS="$threads" CPUSET="$cpuset" \
  WARM_REPS=5 GRAVLAX_REPO="$gravlax_repo" \
  "$scripts_repo/scripts/117_stream_molecule_cram_replay.sh"

{
  printf 'key\tvalue\n'
  printf 'date\t%s\n' "$run_date"
  printf 'dataset\tsc5p_v2_mm_c57bl6_splenocyte_1k\n'
  printf 'threads\t%s\n' "$threads"
  printf 'cpuset\t%s\n' "$cpuset"
  printf 'nofile_soft\t%s\n' "$nofile"
  printf 'gravlax_commit\t%s\n' "$(git -C "$gravlax_repo" rev-parse HEAD)"
  printf 'dirty_worktree\t%s\n' "$(git -C "$gravlax_repo" status --porcelain | wc -l)"
  printf 'aie_binary_sha256\t%s\n' "$(sha256sum "$AIE_BIN" | cut -d' ' -f1)"
  printf 'r1_sha256\t%s\n' "$(sha256sum "$r1" | cut -d' ' -f1)"
  printf 'r2_sha256\t%s\n' "$(sha256sum "$r2" | cut -d' ' -f1)"
  printf 'whitelist_sha256\t%s\n' "$(sha256sum "$whitelist" | cut -d' ' -f1)"
  printf 'gtf_sha256\t%s\n' "$(sha256sum "$gtf" | cut -d' ' -f1)"
  printf 'fasta_sha256\t%s\n' "$(sha256sum "$fasta" | cut -d' ' -f1)"
  printf 'archive_sha256\t%s\n' "$(sha256sum "$OUT_ROOT/archive/mouse5p-1k.aie" | cut -d' ' -f1)"
  printf 'archive_bam_replay_identical\ttrue\n'
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'kernel\t%s\n' "$(uname -srmo)"
} >"$OUT_ROOT/protocol.tsv"

normalized_stream_sha=$(awk -F '\t' '$1 == "normalized_sam_stream" { print $3 }' \
  "$OUT_ROOT/matched/artifact-metrics.tsv")
[[ "$normalized_stream_sha" =~ ^[0-9a-f]{64}$ ]] || {
  echo "invalid normalized SAM-stream digest: $normalized_stream_sha" >&2
  exit 1
}
summary_args=(
  --index-time "$OUT_ROOT/index/time.txt"
  --oracle "$OUT_ROOT/oracle"
  --ingest "$OUT_ROOT/ingest"
  --archive "$OUT_ROOT/archive/mouse5p-1k.aie"
  --archive-time "$OUT_ROOT/archive/time.txt"
  --archive-debug "$OUT_ROOT/archive/debug.txt"
  --replay "$OUT_ROOT/replay-archive"
  --replay-bam "$OUT_ROOT/replay-bam"
  --forward "$OUT_ROOT/replay-forward"
  --velocity "$OUT_ROOT/replay-velocity"
  --velocity-bam "$OUT_ROOT/replay-velocity-bam"
  --matched "$OUT_ROOT/matched"
  --cram-replay "$OUT_ROOT/matched-cram-replay"
  --normalized-stream-sha256 "$normalized_stream_sha"
  --gene-fidelity "$OUT_ROOT/comparisons/gene-fidelity.json"
  --gravlax-commit "$(git -C "$gravlax_repo" rev-parse HEAD)"
  --aie-binary "$AIE_BIN"
  --out "$summary_out"
)
if [[ -n "${ONEPASS_FIDELITY:-}" ]]; then
  [[ -r "$ONEPASS_FIDELITY" ]] || {
    echo "ONEPASS_FIDELITY is not readable: $ONEPASS_FIDELITY" >&2
    exit 2
  }
  summary_args+=(--onepass-fidelity "$ONEPASS_FIDELITY")
fi
"$python" "$scripts_repo/scripts/125_summarize_mouse_5p.py" "${summary_args[@]}"

printf 'mouse 5-prime gate completed: %s\n' "$OUT_ROOT"
printf 'compact result: %s\n' "$summary_out"
