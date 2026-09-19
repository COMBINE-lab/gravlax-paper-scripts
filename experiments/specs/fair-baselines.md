# Experiment 1: function-matched storage and runtime baselines

Status: D0 correctness, randomized five-warm-run Gene-only timing, a non-destructive specific-file
cold observation, post-correction molecule BAM/CRAM storage plus streamed quantification, and
counts-only STARsolo timing on D0/D1/D2' passed on 2026-08-30. See
`results/post-v1-d0-fair-baseline-pilot.json`,
`results/post-v1-d0-counts-only-pilot.json`, and
`results/post-v1-d0-post-correction-baseline.json`; the accepted timing result is
`results/post-v1-d0-randomized-counts.json`, and the matched-CRAM generalization is
`results/post-v1-matched-cram-generalization.json`. Scale counts-only timing is in
`results/post-v1-scale-randomized-counts.json`.

## Hypothesis

At equal retained replay capability, Gravlax remains materially smaller than a minimized BAM/CRAM
and materially faster than replay from those files. Against fresh annotation-aware counting with
unrequested BAM output disabled, Gravlax still has a large end-to-end runtime advantage, but the
published 150–370× ratio will shrink.

## Datasets

- Primary: D0 PBMC 1k (fast iteration and the existing Malva match).
- Scale: D1 PBMC 5k (deep/saturated).
- Generalization: D2' brain nuclei (high archive bits/read and sparse cells).
- Optional confirmation: D2 glioblastoma after the protocol is frozen.

Use the exact FASTQ, genome, whitelist, GTF, and STAR identities recorded in `PROVENANCE.md` and
`experiments/manifests/gates.yaml`.

## Artifact tiers

1. **FASTQ:** original R1+R2 gzip files.
2. **All-tag BAM/CRAM:** the existing annotation-free ingest BAM and CRAM 3.1 archive-mode file.
   This is an information-richer archival comparison, not function matched.
3. **Ingest-equivalent minimized BAM/CRAM:** retain reference id/start, flags, CIGAR, QNAME where
   multimapper grouping needs it, `CR`, `CY`, `UR`, and `NH`; remove bases, cDNA qualities, and
   unused tags. Validate that Gravlax extraction and every current replay/query result are
   unchanged before accepting the tier.
4. **Post-correction replay baseline:** ACCEPTED ON D0, D1, AND D2'. A sequence-free BAM with documented local
   tags represents opaque corrected classes, molecule groups, placements, weights, and explicit
   UMI edges without inventing nucleotide UMIs. CRAM 3.1 is 258.2 MB versus 111.1 MB for `.aie`
   (2.324×); the corresponding D1 and D2' ratios are 2.378× and 2.552×. Call this a
   Gravlax-tagged post-correction molecule CRAM: its replay evidence is
   matched, but generic tools do not understand its relations, it has no `.aie` query indexes, and
   the prototype records rather than embeds the archive's genome signature.
5. **`.aie`:** current v1 archive and any later streaming-reader-compatible archive, with identical
   molecule evidence.

For each tier, publish an explicit capability vector: realignment, barcode re-correction,
sequence/allele query, arbitrary sequence search, compatible-GTF replay, velocity, junction/region
query, and EM. Size ratios are only called function matched when these vectors agree for the
claimed function.

The accepted D0 transform is `scripts/110_make_replay_minimal_bam.sh`: it keeps QNAME for all
records (not only multimappers), sets MAPQ to zero, removes sequence/quality, and retains only
`CR`, `CY`, `UR`, and `NH`. The full reproduction and equivalence checks are in
`scripts/111_d0_fair_baseline_pilot.sh`. This tier is called ingest-equivalent rather than fully
function matched because it still retains raw correction evidence and every alignment instance.

## Runtime arms

- Fresh STARsolo, annotation aware, current published command including sorted BAM output.
- Fresh STARsolo counts-only (`--outSAMtype None`) with otherwise identical alignment/counting
  parameters.
- Gravlax replay from `.aie`.
- Gravlax `replay-rows --from-bam` from all-tag BAM.
- Gravlax replay from ingest-equivalent minimized BAM and CRAM, if CRAM decoding can feed the same
  row path without conversion; otherwise measure and report the conversion explicitly.
- Optional: precompiled annotation plus `.aie` to separate GTF compilation from molecule replay.

Use 24 threads, pinned tool versions, local scratch, and a randomized Latin-square order. Run one
cold-cache replicate after an explicit, recorded cache-control procedure and five warm-cache
replicates. Do not use destructive host-wide cache drops on the shared server. Record wall time,
CPU time, maximum RSS, bytes read/written, output bytes, and exit status with `/usr/bin/time -v`
plus the existing command logs.

The accepted D0 implementation uses a balanced two-arm alternating crossover (seed 20260830
selects the first arm) and `POSIX_FADV_DONTNEED` on only the named inputs before one cold
observation per arm. Every advisory call must succeed, and `/usr/bin/time` must report nonzero
file-system inputs for each cold observation; advisory cache control remains labeled best-effort
because another shared-host process can recache a file.

The accepted D1/D2' implementation adds CPU affinity (CPUs 0--23) to the same five-block schedule.
D1 STAR/native medians are 472.90/5.80 s (81.53×); D2' medians are 326.43/4.76 s (68.58×).
All 24 scale runs are exact and no STAR arm emits alignments. Warm means no explicit eviction:
D1's 57 GB working set was not always fully resident, so retain per-run input counters and the
advisory-cold results (58.36× D1, 61.63× D2') rather than claiming every run was cache-resident.

For quantification from post-correction CRAM, the accepted path is a streaming decoder feeding
`replay-rows --from-molecule-bam` over a named pipe; a materialized BAM is forbidden. Decoder and
quantifier share a 24-CPU affinity set so composition does not double the hardware budget. Warm
medians for CRAM versus matched native `.aie` are 20.65/2.34 s on D0, 97.02/7.13 s on D1, and
86.19/6.29 s on D2': **8.84–13.69×** slower from CRAM. Every matrix artifact is byte-identical.
Report these runtimes beside the **2.32–2.55×** storage ratios; see
`results/post-v1-matched-cram-generalization.json`. D0 CRAM remains 4.05× faster than counts-only
STARsolo under the earlier randomized comparator.

## Correctness controls

- SHA-256 the input identities and every output artifact.
- Require byte-identical archive-versus-BAM Gravlax matrices for the same reduced rows.
- Compare minimized-versus-all-tag extraction summaries and matrices byte-for-byte.
- Compare every arm with fresh STARsolo using moved UMI mass, cell/gene totals, nonzero-set
  disagreement, and a per-feature tail report; never relabel this comparison as exact.
- Fail closed if a minimization step removes information used by any current command.

## Metrics and plots

- Bytes/read, bytes/archived molecule, bytes/called-cell molecule, and compression/decompression
  throughput.
- Wall, CPU, peak RSS, and physical I/O for ingest and replay separately.
- Storage-versus-capability table rather than a single scalar ladder.
- Warm/cold runtime distributions with all replicates shown.
- Scaling line through D0, D1, and D2', with dataset identity visible rather than treating depth as
  the only varying factor. The accepted matched-CRAM measurements are 0.38–0.43 streamed seconds
  per million expanded records.

## Decision gates

- Replace the manuscript storage headline with the function-matched ratio if tier 3 passes all
  correctness controls; otherwise retain only absolute sizes and capability-qualified ratios.
- Replace 150–370× with the accepted 34.31×/81.53×/68.58× counts-only ratios as the main runtime
  comparison. The former ratio may remain only as an end-to-end avoided-work result.
- Advance a custom v2 codec only if it beats the smallest matched baseline by at least 2× or enables
  a capability/RSS/query improvement unavailable to that baseline.
- If minimized CRAM is within 25% of `.aie`, shift the paper's storage story from compression to
  indexed molecular queryability and recomputation.

## Failure interpretation

- A much smaller counts-only speedup invalidates the current headline but not the avoided-work or
  interactive-query story.
- A near-sized minimized CRAM means general-purpose compression is already strong once content is
  matched; the scientific contribution must rest on the evidence abstraction and indexes.
- A minimized BAM/CRAM that cannot preserve the capability vector is not a failed baseline; it is
  evidence that the formats encode a different abstraction and must be described as such.

## Compute estimate and paper value

Protocol pilot on D0: 0.5–1 day. Full D0/D1/D2' run: 2–4 host-days plus roughly 2× the existing
alignment storage temporarily. Paper value is critical: this experiment converts the most obvious
reviewer objection into either a stronger matched claim or an honest reframing before submission.
