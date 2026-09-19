# Gate G-B decision — does annotation version actually move the count matrix?

**Date:** 2026-08-28 · **Host:** analysis-host · **Scale:** full D0 (pbmc_1k_v3, 66,601,887 reads, 1,225 cells)
**Provenance:** `experiments/manifests/gates.yaml` · **Frozen thresholds:** `docs/decisions/2026-08-28-frozen-gates.md`

## Recommendation: **GO**, with the value proposition narrowed to *distant* annotation pairs

No STOP trigger fires. One frozen PASS criterion is met with ~4× margin; the other lands at 2.43%
against a 3% bar. The headline finding is not the marginal number but the **shape** of the effect,
which is decision-relevant and, as far as the prior-art scan found, unpublished for scRNA-seq.

## Results

All figures restricted to the 1,225 real cells (the unfiltered 238k-barcode set is mostly empty
droplets of ambient RNA, on which the two annotations agree trivially).

| Metric | v32 → v49 (far, 2019→2026) | v48 → v49 (near, adjacent) | Frozen bar |
|---|---:|---:|---|
| UMI mass moved, gene union | **2.432%** | **0.369%** | PASS ≥3%, STOP <1% |
| — gained by newer | 3.171% | 0.722% | |
| — lost by older | 1.692% | 0.017% | |
| — genuinely reassigned | 1.692% | 0.017% | |
| — net gain | 1.479% | 0.705% | |
| Expressed genes changing >10% | **21.798%** (5,208/23,892) | **3.866%** (1,073/27,754) | PASS ≥5%, STOP <2% |
| Cells changing >1% of total count | 275/1,225 (22.4%) | 124/1,225 (10.1%) | reported |
| Genes only in the newer release | 19,443 | 14 | reported |
| UMI mass in newer-only genes | 1.030% | 0.000% | reported |
| Genes detected | 23,356 → 27,713 | 27,677 → 27,713 | reported |

`results/gateb-v32-v49-cells.json`, `results/gateb-v48-v49-cells.json`.

## Verdict against the frozen criteria

| Criterion | Threshold | Observed (far pair) | Verdict |
|---|---|---:|---|
| UMI mass moved | PASS ≥3% | 2.432% | **MARGINAL** |
| Expressed genes changing >10% | PASS ≥5% | 21.798% | **PASS** (4.4×) |
| STOP trigger (both <1% moved *and* <2% genes) | — | neither holds | **no STOP** |

The frozen rule makes STOP conjunctive, and neither conjunct holds. Proceeding.

## What the numbers actually say

**Annotation churn is real but concentrated, not diffuse.** Only 2.4% of UMI mass moves between a
2019 and a 2026 annotation, yet **21.8% of expressed genes change by more than 10%**. The effect is
not spread thinly across the transcriptome; it lands hard on a fifth of genes. A 2.4% aggregate
figure would be easy to dismiss, and would be the wrong reading — the per-gene view is what a
biologist would notice, and it is the number that justifies the capability.

**Annotation updates are predominantly additive, not corrective.** In both pairs, mass gained by the
newer release exceeds mass lost by a wide margin (far: 3.17% vs 1.69%; near: 0.72% vs 0.02%). GENCODE
grows into previously unassigned territory far more than it relabels existing assignments. For the
project this is *favourable*: newly annotated genes claim reads that were previously intergenic or
intronic, and those reads are exactly what a count matrix threw away but an annotation-independent
archive still holds. The 19,443 genes existing only in v49, carrying 1.03% of mass, are unreachable
from a v32 count matrix by any means whatsoever.

**Adjacent releases are nearly free.** v48 → v49 moves 0.369% of mass despite adding 31% more
transcripts (385,669 → 507,365) — below the 1% STOP floor on its own. Transcript-level growth barely
touches a 3'-biased gene-level matrix. This bounds the pitch honestly: the capability is worth
having for archives that outlive several years of annotation drift, not for a one-version bump.

## Consequence for Gate G-D

G-D's frozen second criterion is that annotation-independent grouping error must sit **below** the
annotation signal it claims to capture. With the dev-scale G-D result of **0.49%** against the far
pair's **2.43%**, the margin is **~5×**. That criterion is satisfied. Against the near pair (0.369%)
it would *not* be, which is a further reason to state the far pair as the operating regime.

## Caveats

- One dataset, one tissue (PBMC). Annotation churn is tissue-dependent, and PBMC is among the
  best-annotated tissues — so 2.43% is plausibly a **lower** bound. Brain or testis (high lncRNA
  churn) should move considerably more; that is dataset D2's job.
- Gene-level `Gene` (exonic) semantics only. `GeneFull` and Velocyto spliced/unspliced features were
  produced by the same runs and are not yet analysed; intronic-inclusive semantics should be more
  annotation-sensitive, not less.
- STARsolo semantics only. The alevin-fry oracle family is deferred past Gate -1.

## Next

G-C (alignment invariance) is the remaining cheap kill-shot and is now unblocked.
