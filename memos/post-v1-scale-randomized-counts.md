# Counts-only STARsolo versus replay at D1/D2′ scale

Date: 2026-08-30  
Result: `results/post-v1-scale-randomized-counts.json`  
Protocol: `scripts/120_scale_randomized_counts_comparison.sh`

## Verdict

The fair runtime result strengthens at scale. In five alternating two-arm blocks at an equal
24-CPU affinity budget, with STARsolo restricted to Gene counts and `--outSAMtype None`:

| dataset | counts-only STARsolo | native `.aie` replay | wall speedup | peak-RSS ratio |
|---|---:|---:|---:|---:|
| D0 PBMC 1k | 83.71 s | 2.44 s | 34.31× | 5.03× |
| D1 PBMC 5k | 472.90 s | 5.80 s | **81.53×** | 1.92× |
| D2′ brain nuclei | 326.43 s | 4.76 s | **68.58×** | 1.81× |

D0 is the previously accepted randomized result. D1 and D2′ each add five warm/no-explicit-
eviction paired blocks plus one advisory-cold observation per arm. All 36 runs across the three
datasets reproduce their arm-specific frozen references; STAR emits no BAM or SAM.

## Cold and cache interpretation

After accepted `POSIX_FADV_DONTNEED` calls on named inputs, D1 takes 476.24/8.16 s (58.36×) and
D2′ takes 327.28/5.31 s (61.63×) for STAR/replay. Both arms report physical input I/O.

“Warm” must be qualified as *no explicit eviction*. D2′ has negligible input counters after the
cold run, but D1's 57 GB STAR working set is not reliably resident on the shared host: some warm
runs perform physical reads. This does not weaken the measured distribution—the STAR CV is only
1.25%—but it forbids a blanket fully-warm-cache label. The cold observation and randomized series
should be reported together.

## Context for the reduced-CRAM result

Using the separately accepted streamed molecule-CRAM medians, CRAM quantification remains faster
than fresh counts-only STARsolo: 4.05× on D0, 4.87× on D1, and 3.79× on D2′. The native index is
another 8.84–13.69× faster than the content-matched CRAM path while using 2.32–2.55× fewer bytes.
This supplies a clean three-level story: fresh alignment is slowest; reduced CRAM avoids alignment
but pays explicit-record decode/reconstruction; `.aie` directly loads the molecular abstraction.

## Claim decision

- Replace the D0-only main runtime sentence with **34–82× faster across three datasets** and retain
  the exact medians in the Results section.
- Keep memory claims dataset-specific: the relative benefit is 5.03× on D0 but 1.81–1.92× at
  scale because native replay memory grows with the molecule set while STAR's genome index is
  nearly fixed.
- The fair-baseline experiment is no longer binding. The next paper-changing gate is controlled
  sequencing-depth scaling, which can disentangle read depth, molecule count, saturation, archive
  size, replay time, and fresh-count time within one biological sample.

