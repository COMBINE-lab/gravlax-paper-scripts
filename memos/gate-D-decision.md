# Gate G-D decision — annotation-independent UMI grouping

**Date:** 2026-08-28 · **Host:** analysis-host · **Scale:** full D0 (66,601,887 reads; 1,225 cells; 70.0% saturation)
**Provenance:** `experiments/manifests/gates.yaml` · **Frozen thresholds:** `docs/decisions/2026-08-28-frozen-gates.md`
**Supersedes:** `memos/gate-D-interim-dev.md`, whose dev-scale conclusions were wrong in two ways — recorded below rather than quietly corrected.

## Recommendation: **PASS**, at locus gap 50 kb with UMI error correction

## Question

STARsolo collapses UMIs per (cell, gene) — annotation-dependent by construction
(`SoloFeature_collapseUMIall.cpp`). The archive must collapse without an annotation. How much error
does that introduce, and does it stay below the annotation signal it claims to capture?

Both groupings run on the **same alignments** (the v49 oracle BAM), isolating grouping from
alignment. G2 uses raw `UR`, never STARsolo's `UB`, whose correction is itself per-gene (decision D9).

## Full-scale results

| locus gap (bp) | UMI 1MM corr. | G1 molecules | G2 molecules | median per-cell L1 | p90 | UMIs merged | multi-gene |
|---:|:--|---:|---:|---:|---:|---:|---:|
| 2,000 | yes | 10,191,963 | 10,565,453 | 3.945% | 4.654% | 558,691 | 1,678 |
| 10,000 | **no** | 10,191,963 | 10,875,482 | **7.141%** | 7.956% | 0 | 24 |
| 10,000 | yes | 10,191,963 | 10,308,775 | 1.215% | 1.438% | 566,707 | 1,697 |
| 50,000 | yes | 10,191,963 | 10,282,543 | **0.939%** | 1.139% | 567,635 | 1,819 |
| 500,000 | yes | 10,191,963 | 10,278,423 | 0.900% | 1.102% | 568,096 | 2,178 |

`results/gated-full-*.json`.

## Verdict against the frozen criteria

| Criterion | Threshold | Observed (50 kb + correction) | Verdict |
|---|---|---:|---|
| Median per-cell L1 relative error | PASS < 1% | 0.939% | **PASS** |
| Error vs the A1↔A2 annotation delta (G-B far pair, 2.432%) | strictly smaller | 0.939% (2.6× margin) | **PASS** |
| STOP trigger (error ≥ annotation delta) | — | does not hold | **no STOP** |

## Two dev-scale conclusions that did not survive

The 4M-read subsample was at **17.2% sequencing saturation**; full D0 is at **70.0%**. That single
difference reversed both of the interim memo's secondary findings, and it is the reason a dev-scale
pass is not evidence of a full-scale pass:

1. **UMI error correction went from negligible to essential.** At dev scale it bought 0.13
   percentage points and the memo called it unimportant. At full depth it is the difference between
   **7.141% and 1.215%** at the same locus gap — a 6× reduction. At 3.4 reads per UMI there are
   simply far more UMI sequencing errors to absorb; at 1.2 reads per UMI there were almost none.
2. **The recommended operating point was wrong.** Dev scale plateaued by 10 kb and the memo proposed
   10 kb; at full depth 10 kb gives 1.215% (a fail) and the curve does not flatten until ~50 kb.

The headline error itself roughly doubled, 0.49% → 0.94%. Anyone reading the dev number as a
prediction of the full-scale number would have been wrong by 2×.

## What the numbers say

**The error remains entirely over-splitting.** G2 exceeds G1 at every setting; annotation-free
grouping never merges molecules the oracle separates, it only fails to merge some the oracle joins.
The residual ~90,580 molecules at 50 kb (0.89%) are cases where STARsolo merges UMIs across an
entire gene irrespective of genomic distance — no distance-based rule reproduces that.

**Over-merging never becomes binding.** Multi-gene molecules stay at 1,819 of 10.3M (0.018%) at
50 kb, and only reach 2,178 (0.021%) at 500 kb — effectively one bucket per (cell, chromosome,
strand). Within a cell a UMI is very nearly a unique key, which is why the curve flattens rather
than turning back up, and why widening the window is safe rather than reckless.

**50 kb is the honest operating point, not 500 kb.** 500 kb is marginally better (0.900% vs 0.939%)
but is no longer meaningfully a "locus" — it spans whole chromosome arms. 50 kb comfortably contains
a typical human gene (median ~24 kb), so it retains a defensible biological reading while sitting
under the threshold. **Adopted: locus gap 50 kb, 1-mismatch UMI correction on.**

## Caveats

- One dataset at one depth. Since the error is saturation-driven, a deeper library should be worse.
  **This must be re-measured on D1/D2 rather than assumed**, and saturation reported alongside it.
- The 1-mismatch correction is a simple most-abundant-first greedy. A directional/adjacency method
  (UMI-tools) would likely close part of the residual; that is the obvious next improvement if
  more margin is wanted.
- Measured against STARsolo semantics only.

## Consequence

The margin against the annotation signal is **2.6×**, not the ~5× the interim memo projected. That
is still the right side of the frozen line, but it is thin enough that the near-pair regime
(v48→v49, 0.369% annotation delta) is now definitively **out of scope**: our grouping error alone
would exceed the entire signal there. Gate G-B's conclusion that the value proposition lives in
distant annotation pairs is reinforced, and is now load-bearing rather than a caveat.
