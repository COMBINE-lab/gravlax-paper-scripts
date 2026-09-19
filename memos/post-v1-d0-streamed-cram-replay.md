# Quantification directly from the reduced molecule CRAM

Date: 2026-08-30  
Result: `results/post-v1-d0-streamed-cram-replay.json`  
Protocol: `scripts/117_stream_molecule_cram_replay.sh`

## What is required

The 258.2 MB post-correction molecule CRAM is sufficient for the tested Gene quantification, but
not by generic CRAM semantics alone. Replay requires a decoder that preserves record order and the
Gravlax local tags, the exact reference FASTA used for CRAM, the tag-aware molecule importer, a
compatible GTF, and the output barcode ordering. It does not require FASTQ, a STAR index, raw
barcode correction, a whitelist, nucleotide UMI strings, or a materialized BAM.

The benchmark streams `samtools view -b` into `aie replay-rows --from-molecule-bam` over a named
pipe. Decoder and quantifier share CPUs 0--23, enforcing one combined 24-hardware-thread ceiling.
Every output is compared byte-for-byte with the frozen D0 matrix, features, and barcodes.

## Runtime result

Five warm end-to-end runs take 20.35--21.25 s, median **20.65 s** (CV 1.63%). The corresponding
native `.aie` median is 2.44 s, so quantification from molecule CRAM is **8.46× slower**. It is
still **4.05× faster** than the 83.71 s counts-only Gene STARsolo median.

The operational context for the storage result is therefore:

| representation | bytes | warm quantification | relative to `.aie` |
|---|---:|---:|---:|
| `.aie` | 111,087,381 | 2.44 s | 1× size, 1× time |
| post-correction molecule CRAM | 258,151,959 | 20.65 s | 2.32× size, 8.46× time |

Streaming is materially better than the naive two-step workflow: the separately measured CRAM
decode plus materialized-BAM replay totals about 29.5 s, whereas overlap through the pipe finishes
in 20.65 s and writes no 685 MB temporary BAM. It remains 21% slower than replay from an already
materialized molecule BAM (17.03 s).

## Where the time goes

The replay logs report roughly 18.0--18.4 s in load/reconstruction and only about 1.7 s for GTF
loading plus assignment/counting. `samtools` itself runs for a median 18.02 s and overlaps the
consumer. Thus the main cost is expanding and validating 48.2 million tag/CIGAR records into
22.2 million molecules, not the biological quantification step. A native CRAM reader could remove
the pipe and external process, but it cannot remove the semantic verbosity of the explicit record
representation; that distinction should guide whether such engineering is worth pursuing.

One input-specific cold run takes 20.74 s after successful advisory eviction of 6.85 GB across
five named inputs. Both processes report nonzero input-I/O counters. The cold/warm similarity is
plausible because CRAM decompression and record reconstruction dominate, but this remains one
advisory shared-host observation.

## Claim decision

Pair the D0 2.32× storage result with the 8.46× replay-runtime result. Do not imply that reduced
CRAM is merely a somewhat larger archive with comparable interactivity. Conversely, do not call it
non-replayable: the streaming experiment proves exact quantification without raw reads or an
intermediate BAM. Generalize both axes together on D1 and D2′.
