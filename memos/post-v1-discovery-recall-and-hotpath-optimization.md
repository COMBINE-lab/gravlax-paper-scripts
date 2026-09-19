# Post-v1 discovery recall and hot-path optimization

**Date:** 2026-08-31  
**Protocol locks:** `c6edf01`, `86193f5`, `4ca25df`, `023eef4`  
**Implementation:** `c1ef9c4d8a4f65c762bec118c98a79f0d746d512`  
**Verdict:** residual-site recall **PASS**; default discovery and junction-query optimization
**PASS**; replay optimization **NO CHANGE**

## Complete-denominator discovery

The original span-claiming rule remains the compatibility default and reproduces 164/367
(44.7%) genes in the complete v32-to-v49 denominator. Its failure is informative: any overlap
with a v32 transcript suppresses a molecule, even when its stored exon/junction evidence is
incompatible with that transcript.

Three claiming ablations did not solve the problem cleanly. Strand-aware spans recover 200/367
but emit 39,969 D0 candidates, above the locked threefold volume guard. Direct compatibility
claiming produces broad transitive components and recovers only 57/367 under deterministic
one-to-one matching. Both results remain in the record.

The successful mode keeps every original 10-UMI span candidate, then adds a second channel from
span-suppressed but transcript-incompatible molecules. It clusters their transcript-oriented
terminal bases within a non-transitive 1-kb window. Candidate geometry is the observed terminal
extent, not the full alignment span. The diagnostic 10-UMI residual channel was too large, so the
predeclared support grid was evaluated without interpolation:

| residual support | D0 candidates | ratio | genes | recall | UMI recall | D1 recurrence | guard |
|---:|---:|---:|---:|---:|---:|---:|---|
| 25 | 79,356 | 6.41x | 331/367 | 90.2% | 93.0% | 99.5% | fail volume |
| 50 | 39,578 | 3.19x | 307/367 | 83.7% | 88.4% | 99.2% | fail volume |
| **75** | **26,787** | **2.16x** | **266/367** | **72.5%** | **81.5%** | **99.0%** | **pass** |
| 100 | 21,054 | 1.70x | 237/367 | 64.6% | 73.3% | 98.9% | pass, lower recall |

At the selected threshold, recall is 51/52 clean, 52/60 antisense-only, and 163/255 same-strand.
This is a +27.8 percentage-point gain over the conservative rule while satisfying the locked
volume and recurrence safeguards. All 14,398 added D0 residual candidates are at most 1,001 bp
(median 926 bp; p95 1,001 bp), preventing long intervals from gaming interval-overlap recall.

This is a recall result, not a precision result. GENCODE v49 is used as later-annotation truth for
the 367-gene denominator. Lack of v49 overlap does not make another emitted locus a biological
false positive; it remains unresolved. The selected residual threshold is targeted to the human
3' data and should not silently become a universal default for other assays.

## Query and discovery performance

Sharing immutable archive shape/pattern tables instead of deep-cloning them improves the D1
whole-chromosome junction query with exact cells from 0.51 to 0.44 s median and reduces median
peak RSS from 1,189,836 to 1,048,624 KiB. The 5,089-row output is byte-identical. Small D0 region
and junction probes improve in the same direction but remain below 0.2 s and are not paper timing
claims.

For default discovery, parallel decode/classification now pairs sparse access units but processes
dense units singly. A pair contains at most 500,000 molecules, so the extra decode working set is
bounded rather than growing with archive length. D0 improves from 3.36 to 2.30 s (+63 MiB median
RSS); D1 improves from 15.07 to 10.67 s (+81 MiB). The speedups are 31.5% and 29.2%, and both TSVs
are byte-identical to the baseline. The D0 percentage RSS change exceeds the original 10% guard,
but the absolute, explicitly capped 63 MiB trade was accepted after review.

## Streaming re-quantification

No new replay change was retained. Fusing decode with classification left D0 unchanged and
regressed D1 from 7.12 to 8.42 s. In-place reducer compaction changed wall time by about 1% with no
stable memory reduction. Both changes were reverted. The existing bounded replay remains the
best measured implementation: 2.81 s / 1,358,860 KiB on D0 and 7.12 s / 5,375,848 KiB on D1 for
the fixed replay gate. Further replay work should begin with a different representation or
parallel decomposition, not more rearrangement of the current inner loop.

## Reproduction

Build release binaries at baseline `0efabed0` and candidate `c1ef9c4d`, then run:

```sh
AIE_BIN=/path/to/candidate/aie \
  bash scripts/140_discovery_recall_optimization.sh

BASELINE_AIE_BIN=/path/to/baseline/aie \
AIE_BIN=/path/to/candidate/aie \
  bash scripts/142_hotpath_optimization_gate.sh
```

The compact records are `results/post-v1-discovery-recall-optimization.json` and
`results/post-v1-hotpath-optimization.json`. The selected full evaluator output, truth audit, and
candidate table are committed beside them; raw arm TSVs and GNU time logs remain regenerated
run artifacts.
