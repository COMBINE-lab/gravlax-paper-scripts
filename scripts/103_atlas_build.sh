set -x
P=${GRAVLAX_PROJECT_ROOT:?set GRAVLAX_PROJECT_ROOT to the analysis project root (see README)}
cd $P
AIE=src/target/release/aie
for D in d4 d3; do
  if [ $D = d3 ]; then FQ=data/d3/pbmc_10k_v3_fastqs; else FQ=data/d4/500_PBMC_3p_LT_Chromium_X_fastqs; fi
  R1=$(ls $FQ/*_R1_*.fastq.gz | paste -sd,)
  R2=$(ls $FQ/*_R2_*.fastq.gz | paste -sd,)
  mkdir -p runs/ingest/$D
  STAR --runThreadN 24 --genomeDir ref/star-base --twopassMode Basic \
    --readFilesIn "$R2" "$R1" --readFilesCommand zcat \
    --soloType CB_UMI_Simple --soloCBwhitelist data/3M-february-2018.txt \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 --soloFeatures SJ \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --outSAMattributes NH HI AS nM CR UR CY UY \
    --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 32000000000 \
    --outSAMunmapped Within --outFileNamePrefix runs/ingest/$D/ \
    > logs/$D-ingest.log 2>&1 || { echo "$D STAR FAILED"; exit 1; }
  $AIE ingest-archive runs/ingest/$D/Aligned.sortedByCoord.out.bam \
    --whitelist data/3M-february-2018.txt --out runs/archive/$D.aie --zstd-level 19 \
    2>&1 | grep extracted
  stat -c "${D}_SIZE %s" runs/archive/$D.aie
  echo "${D}_ARCHIVE_DONE"
done
echo "---- atlas federate: 6 archives ----"
/usr/bin/time -f "federate6_wall=%e" $AIE federate \
  runs/archive/d0.aie runs/archive/d1.aie runs/archive/d3.aie runs/archive/d4.aie \
  runs/archive/d2.aie runs/archive/d2p.aie chr16:89562391-89562883 --top 2
echo ATLAS_BUILD_DONE
