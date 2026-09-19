# Matched molecule-CRAM storage and quantification across three datasets

Date: 2026-08-30  
Result: `results/post-v1-matched-cram-generalization.json`  
Protocols: `scripts/117_stream_molecule_cram_replay.sh`,
`scripts/118_post_correction_storage.sh`, and `scripts/119_native_aie_replay.sh`

## Verdict

The D0 finding generalizes. A post-correction molecule CRAM can drive exact Gene quantification
without FASTQ, realignment, raw barcode correction, or a materialized BAM. It is not a
near-equivalent substitute for the native index: across D0, D1, and D2′, CRAM is **2.32–2.55×
larger** and its streamed quantification is **8.84–13.69× slower**.

| dataset | `.aie` | molecule CRAM | size ratio | native `.aie` replay | streamed CRAM replay | time ratio |
|---|---:|---:|---:|---:|---:|---:|
| D0 PBMC 1k | 111.1 MB | 258.2 MB | 2.32× | 2.34 s | 20.65 s | 8.84× |
| D1 PBMC 5k | 544.0 MB | 1,293.8 MB | 2.38× | 7.13 s | 97.02 s | 13.61× |
| D2′ brain nuclei | 595.2 MB | 1,519.1 MB | 2.55× | 6.29 s | 86.19 s | 13.69× |

Each timing is the median of five warm runs under the same 24-CPU affinity set. One advisory-cold
run per arm is similar to warm (CRAM: 20.74, 97.68, and 86.35 s; `.aie`: 2.63, 7.54, and 6.75 s).
All 36 timed outputs reproduce the frozen matrix, features, and barcodes byte-for-byte. D1 and D2′
BAM/CRAM normalized record streams also hash identically.

## What CRAM replay actually requires

The operational path is `samtools view -b` into `aie replay-rows --from-molecule-bam` over a named
pipe. It requires the exact reference FASTA used by CRAM, a decoder that preserves record order and
the documented Gravlax local tags, the tag-aware molecule importer, a compatible GTF, and the
barcode output ordering. It does not require a STAR index, FASTQ, raw whitelist correction,
nucleotide UMI strings, or a BAM temporary file.

This is a standards-compliant CRAM container but a Gravlax-specific evidence contract. Generic
tools expose records and tags; they do not interpret molecule groups, opaque UMI classes, weighted
placements, or UMI-graph edges as a quantification-ready abstraction.

## Where the latency comes from

CRAM replay expands 48.2M, 253.1M, and 224.3M records to reconstruct 22.2M, 108.7M, and 130.3M
molecules. Its warm cost is strikingly stable at 0.38–0.43 seconds per million records. Median
load/reconstruction alone takes 18.3, 91.3, and 80.6 s; native `.aie` load takes 0.4, 2.5, and
1.8 s. The GTF and assignment tail is only a few seconds. This identifies explicit-record decode
and reconstruction—not the biological counting policy—as the bottleneck.

The separate decoder and replay-process RSS maxima must not be mistaken for a synchronized peak.
Their sums are conservative upper bounds of 14.6, 28.5, and 31.6 GB, versus native replay medians
of 6.94, 18.45, and 19.01 GB. A future benchmark should sample the process group if a pipeline peak
RSS claim is needed.

## Claim and engineering decisions

- Use **2.32–2.55× smaller and 8.84–13.69× faster** as the three-dataset function-matched result.
- Preserve the larger 9.0–12.7× comparison only for capability-qualified, tag-preserving ingest
  CRAM; it retains substantially more information.
- Say that reduced CRAM is sufficient and replayable, not that CRAM is incapable of the task.
- Do not prioritize a native CRAM reader merely to remove the pipe: it might trim orchestration
  overhead, but the measured critical path is explicit record expansion. A native reader is worth
  pursuing only if profiling demonstrates a structural way to avoid that reconstruction.
- The next binding fair-baseline gate is D1/D2′ counts-only STARsolo versus native `.aie` replay.

