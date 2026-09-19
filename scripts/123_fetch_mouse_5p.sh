#!/usr/bin/env bash
# Fetch the public 10x mouse 5' v2 pilot, its whitelist, and current GENCODE M39.
set -euo pipefail

: "${PROJ_ROOT:?set PROJ_ROOT}"

data_dir=${DATA_DIR:-$PROJ_ROOT/data/mouse5p-1k}
ref_dir=${REF_DIR:-$PROJ_ROOT/ref/mouse-m39}
bundle=$data_dir/sc5p_v2_mm_c57bl6_splenocyte_1k_fastqs.tar
fastq_dir=$data_dir/sc5p_v2_mm_c57bl6_splenocyte_1k_5gex_fastqs
whitelist=$data_dir/737K-august-2016.txt
gtf_gz=$ref_dir/gencode.vM39.primary_assembly.annotation.gtf.gz
fasta_gz=$ref_dir/GRCm39.primary_assembly.genome.fa.gz

mkdir -p "$data_dir" "$ref_dir"

fetch() {
  local url=$1
  local out=$2
  local expected=$3
  if [[ ! -e "$out" ]]; then
    curl -fL --retry 5 --retry-delay 5 -o "$out" "$url"
  fi
  printf '%s  %s\n' "$expected" "$out" | sha256sum --check --status || {
    echo "checksum mismatch: $out" >&2
    exit 1
  }
}

fetch \
  https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_mm_c57bl6_splenocyte_1k/sc5p_v2_mm_c57bl6_splenocyte_1k_fastqs.tar \
  "$bundle" c9da2d7af946853f30dc122d16d23de93d6ae1aa958c4dba3c2bcaa095212eac
fetch \
  https://raw.githubusercontent.com/10XGenomics/supernova/master/tenkit/lib/python/tenkit/barcodes/737K-august-2016.txt \
  "$whitelist" b0dda4b114d8fc7ea8def94fb472e6205c1831a72e3ac245ed13cc332d9f20d1
fetch \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39/gencode.vM39.primary_assembly.annotation.gtf.gz \
  "$gtf_gz" d6da97913ce30f99883fc1216b111569f9947cf203886a5afb607b59228574d4
fetch \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39/GRCm39.primary_assembly.genome.fa.gz \
  "$fasta_gz" 82724ccc3726f7bbc3d1cf46a86bbda753f05c1c1b27772e3ff87641fc71669f

if [[ ! -d "$fastq_dir" ]]; then
  if tar -tf "$bundle" | awk '$0 ~ /^\// || $0 ~ /(^|\/)\.\.($|\/)/ { bad=1 } END { exit bad ? 0 : 1 }'; then
    echo "unsafe path in FASTQ tar" >&2
    exit 1
  fi
  tar -xf "$bundle" -C "$data_dir"
fi

gzip -t "$gtf_gz" "$fasta_gz" "$fastq_dir"/*_R1_001.fastq.gz "$fastq_dir"/*_R2_001.fastq.gz

if [[ ! -e "${gtf_gz%.gz}" ]]; then
  gzip -dc "$gtf_gz" >"${gtf_gz%.gz}.partial"
  mv "${gtf_gz%.gz}.partial" "${gtf_gz%.gz}"
fi
if [[ ! -e "${fasta_gz%.gz}" ]]; then
  gzip -dc "$fasta_gz" >"${fasta_gz%.gz}.partial"
  mv "${fasta_gz%.gz}.partial" "${fasta_gz%.gz}"
fi

printf 'mouse 5-prime inputs verified under %s and %s\n' "$data_dir" "$ref_dir"
