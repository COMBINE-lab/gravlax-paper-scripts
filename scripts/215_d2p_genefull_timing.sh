#!/usr/bin/env bash
# Archived from the analysis host as scripts/150_d2p_genefull_timing.sh (2026-09-11); renumbered here to follow the repository sequence.
# D2' GeneFull timing: STARsolo (GeneFull only, no alignment output) versus native GeneFull replay.
# Same protocol as gravlax-paper-scripts/scripts/120: one advisory-cold observation per arm, then
# five warm blocks in a balanced alternating crossover (STAR first in odd blocks), pinned to 24 CPUs.
set -euo pipefail
P=${GRAVLAX_PROJECT_ROOT:?set GRAVLAX_PROJECT_ROOT to the analysis project root (see README)}
AIE=$P/runs/gravlax-release-0.2.2-20260910/aie-0.2.2
STAR=$P/env/sc/bin/STAR
PY=${GRAVLAX_PYTHON:-python3}
OUT=$P/runs/paper-followup-20260911/d2p-genefull-timing
ARCHIVE=$P/runs/archive/d2p.aie
GTF=$P/annotations/gencode.v49.annotation.gtf
GENOME=$P/ref/star-49
WL=$P/data/3M-february-2018.txt
BARC=$P/runs/oracle/d2p/v49/Solo.out/Gene/raw/barcodes.tsv
ORACLE=$P/runs/oracle/d2p/v49/Solo.out
REF_REPLAY=$P/runs/genefull-20260910/replay-v49-GeneFull
FQ=$P/data/d2p/Brain_3p_fastqs
R1=$FQ/mm_ST_7500_S6_L001_R1_001.fastq.gz,$FQ/mm_ST_7500_S6_L002_R1_001.fastq.gz,$FQ/mm_ST_7500_S6_L003_R1_001.fastq.gz,$FQ/mm_ST_7500_S6_L004_R1_001.fastq.gz
R2=$FQ/mm_ST_7500_S6_L001_R2_001.fastq.gz,$FQ/mm_ST_7500_S6_L002_R2_001.fastq.gz,$FQ/mm_ST_7500_S6_L003_R2_001.fastq.gz,$FQ/mm_ST_7500_S6_L004_R2_001.fastq.gz
THREADS=24; CPUSET=0-23; BLOCKS=5
[[ ! -e $OUT ]] || { echo "refusing to overwrite $OUT" >&2; exit 2; }
mkdir -p $OUT; export RAYON_NUM_THREADS=$THREADS
{
  printf 'date\t%s\ndataset\tD2-prime brain nuclei\ncounting_model\tGeneFull\nthreads\t%s\ncpuset\t%s\nwarm_blocks\t%s\n' "$(date -I)" $THREADS $CPUSET $BLOCKS
  printf 'aie_version\t%s\naie_sha256\t%s\narchive_sha256\t%s\nstar_version\t%s\nstar_features\tGeneFull\nstar_alignment_output\tNone\nannotation\t%s\n' \
    "$($AIE --version)" "$(sha256sum $AIE | cut -d' ' -f1)" "$(sha256sum $ARCHIVE | cut -d' ' -f1)" "$($STAR --version)" "$GTF"
} > $OUT/protocol.tsv
printf 'phase\tblock\tposition\tarm\tlabel\tstatus\n' > $OUT/schedule.tsv

run_aie() {
  local label=$1 d=$OUT/$1; mkdir -p $d/replay
  /usr/bin/time -v -o $d/time.txt taskset -c $CPUSET $AIE replay-rows $ARCHIVE --gene-full \
    --gtf $GTF --barcodes $BARC --out-dir $d/replay --report-format json --report-output $d/report.json \
    > $d/stdout.txt 2> $d/stderr.txt
  local st=exact
  for a in matrix.mtx barcodes.tsv features.tsv; do cmp -s $d/replay/$a $REF_REPLAY/$a || st=DIFFERS_FROM_REFERENCE; done
  echo $st > $d/status.txt; echo "$label $st"
}
run_star() {
  local label=$1 d=$OUT/$1; mkdir -p $d
  /usr/bin/time -v -o $d/time.txt taskset -c $CPUSET $STAR --runThreadN $THREADS --genomeDir $GENOME \
    --readFilesIn $R2 $R1 --readFilesCommand zcat --soloType CB_UMI_Simple --soloCBwhitelist $WL \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
    --soloFeatures GeneFull --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR \
    --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter EmptyDrops_CR --outSAMtype None --outSAMunmapped Within \
    --outFileNamePrefix $d/ > $d/stdout.txt 2> $d/stderr.txt
  local st=exact
  for a in raw/matrix.mtx raw/features.tsv raw/barcodes.tsv filtered/matrix.mtx filtered/barcodes.tsv; do
    cmp -s $d/Solo.out/GeneFull/$a $ORACLE/GeneFull/$a || st=DIFFERS_FROM_ORACLE; done
  rm -rf $d/_STARtmp; echo $st > $d/status.txt; echo "$label $st"
}
rec() { printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$(cat $OUT/$5/status.txt)" >> $OUT/schedule.tsv; }

$PY $P/gravlax-paper-scripts/scripts/115_fadvise_dontneed.py --log $OUT/cold-aie-cache.tsv $ARCHIVE $GTF $BARC
run_aie cold-aie; rec cold 0 1 aie cold-aie
$PY $P/gravlax-paper-scripts/scripts/115_fadvise_dontneed.py --log $OUT/cold-star-cache.tsv ${R1//,/ } ${R2//,/ } $GENOME/Genome $GENOME/SA $GENOME/SAindex
run_star cold-star; rec cold 0 1 star cold-star
for b in $(seq 1 $BLOCKS); do
  if (( b % 2 == 1 )); then first=star; second=aie; else first=aie; second=star; fi
  for pos in 1 2; do
    arm=$first; (( pos == 2 )) && arm=$second
    label=warm-b$b-p$pos-$arm
    if [[ $arm == star ]]; then run_star $label; else run_aie $label; fi
    rec warm $b $pos $arm $label
  done
done
$PY - "$OUT" <<'PYEOF'
import sys,glob,os,statistics,re
out=sys.argv[1]
def wall(f):
    t=open(f).read(); m=re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.*)",t); h=m.group(1).split(":")
    s=float(h[-1])+60*float(h[-2])+(3600*float(h[-3]) if len(h)==3 else 0)
    r=re.search(r"Maximum resident set size \(kbytes\): (\d+)",t); return s,int(r.group(1))/1048576
rows={}
for d in sorted(glob.glob(f"{out}/warm-*")):
    arm=d.split("-")[-1]; s,g=wall(f"{d}/time.txt"); rows.setdefault(arm,[]).append((s,g))
with open(f"{out}/summary.tsv","w") as f:
    f.write("arm\tn\tmedian_wall_s\tmedian_peak_gib\twalls\n")
    for arm,v in rows.items():
        f.write(f"{arm}\t{len(v)}\t{statistics.median(x for x,_ in v):.2f}\t{statistics.median(g for _,g in v):.2f}\t{','.join(f'{x:.2f}' for x,_ in v)}\n")
    if "star" in rows and "aie" in rows:
        f.write(f"speedup_median\t{statistics.median(x for x,_ in rows['star'])/statistics.median(x for x,_ in rows['aie']):.2f}\n")
print(open(f"{out}/summary.tsv").read())
PYEOF
echo "D2P GENEFULL TIMING DONE"
