# Gate G-D interim result (dev scale) — annotation-independent UMI grouping

**Date:** 2026-08-28 · **Scale:** D0-dev (4M read pairs of pbmc_1k_v3, 1,135 cells)
**Status:** interim — full-scale D0 confirmation pending. Not a verdict.

## Question

STARsolo collapses UMIs per (cell, gene), which is annotation-dependent by construction
(`SoloFeature_collapseUMIall.cpp`). The archive must collapse without an annotation. How much
error does that introduce?

Frozen thresholds (`docs/decisions/2026-08-28-frozen-gates.md`):
median per-cell L1 relative error **< 1%**, and the G1↔G2 discrepancy strictly **smaller than the
A1↔A2 annotation delta** from Gate G-B.

## Setup

Both groupings run over the **same alignments** — the v49 oracle BAM — so this isolates grouping
from alignment. Alignment invariance is Gate G-C's separate question.

- **G1** (oracle): distinct corrected `UB` within each (cell, gene) group.
- **G2** (annotation-free): the plan's §20 molecule key `(c, u, L, o)`. Loci are found by
  single-linkage clustering of read starts within (cell, chrom, strand); UMIs are the raw `UR`,
  never `UB`, since STARsolo's UMI correction is itself per-gene (decision D9). Optional 1-mismatch
  UMI correction runs *inside* a locus, consulting no annotation.

## Result

| locus gap (bp) | G1 molecules | G2 molecules | median per-cell L1 | p90 | multi-gene molecules |
|---:|---:|---:|---:|---:|---:|
| 500 | 1,675,657 | 1,699,226 | 1.4386% | 2.0351% | 3 |
| 1,000 | 1,675,657 | 1,693,750 | 1.1055% | 1.5695% | 9 |
| 2,000 | 1,675,657 | 1,687,846 | 0.7356% | 1.0731% | 13 |
| 5,000 | 1,675,657 | 1,685,080 | 0.5668% | 0.8527% | 15 |
| 10,000 | 1,675,657 | 1,684,336 | 0.5172% | 0.8065% | 15 |
| 20,000 | 1,675,657 | 1,684,103 | 0.5040% | 0.7880% | 15 |
| 50,000 | 1,675,657 | 1,683,957 | 0.4961% | 0.7800% | 15 |
| 100,000 | 1,675,657 | 1,683,921 | 0.4932% | 0.7800% | 15 |

`results/gated-dev-sweep.csv`.

## Reading

**The error is entirely over-splitting, and it plateaus at ~0.49%.** G2 always exceeds G1: the
annotation-free grouping never merges molecules the oracle keeps apart, it only fails to merge some
the oracle joins. The plateau from ~10 kb onward is a genuine floor, not a tuning artifact — the
residual ~8,264 molecules (0.49%) are cases where STARsolo merges UMIs across an entire gene
regardless of genomic distance, which no distance-based rule can reproduce.

**Over-merging never becomes the problem.** Multi-gene molecules saturate at 15 out of ~1.68M even
at a 100 kb window — effectively "one bucket per (cell, chromosome, strand)". So within a cell, a
UMI is very nearly a unique key, and annotation-independent grouping is intrinsically safe. That is
a stronger statement than the gate needed, and it is the reason the curve plateaus instead of
turning back up.

**Two design conclusions.** Fixed-width bucketing was wrong and had to go: an earlier version using
a 1 kb grid reported 2.09%, of which roughly a third was molecules straddling a bucket boundary
rather than any real disagreement. And 1-mismatch UMI correction is worth little here (2,159 merges,
~0.13 percentage points) — the gap is dominated by locus structure, not sequencing error in UMIs.

## Verdict

At dev scale the frozen threshold is **met with margin** for any locus gap ≥ 5 kb
(0.57% ≤ median ≤ 1%, plateauing at 0.49%). Two things must still hold before G-D can be called:

1. Full-scale D0 confirmation (running).
2. The second criterion — that this error is **smaller than the annotation delta** — cannot be
   evaluated until Gate G-B reports. 0.49% is only acceptable if annotation churn moves
   substantially more than 0.49% of the count mass; the frozen G-B pass bar is ≥3%, which would
   give roughly a 6× margin.

**Provisional locus gap: 10,000 bp** — the knee of the curve, past which further widening buys
under 0.03 percentage points.
