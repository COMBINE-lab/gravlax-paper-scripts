# Hierarchical EM gate result

**Date:** 2026-08-31  
**Verdict:** STOP for production/default use; retain the optional evaluator and the negative result  
**Implementation:** Gravlax `7813c0b1c755eea07abb0d21e2296f8cc17091ee`  
**Registered run scripts:** `ec0c09f8d2350174d70079ec990a8d609abf442a`  
**Compact result:** `results/post-v1-hierarchical-em.json`

> **Post-hoc qualification (2026-08-31):** this memo records the registered bounded-grid result
> and its STOP verdict unchanged. A subsequent scale audit showed that the selected 80/20 masses
> were diluted by whole-transcriptome normalization and did not test the intended strength of
> partial pooling. See `memos/post-v1-hierarchical-em-scale-audit.md`; hierarchy itself did not
> fail, and the recommendation below for a structurally different model was too strong.

## What was tested

The packed evaluator gained optional query-time barcode/group metadata, group-only and additive
cell/group/global modes, called-cell scoring masks, exact one-group collapse controls, candidate
gene export, and JSON output with top-1, truth probability, negative log loss, multiclass Brier,
and calibration counts. Production `--mask 0 --emit` remains pooled.

Groups were constructed from the base Gravlax Gene matrix with no EM layer. Across mask seeds 7,
17, and 29, every gene occurring in an evaluation candidate set was excluded before replay-only
HVG selection and PCA. This removed 5,663 D0, 8,500 D1, and 8,800 D2-prime genes while leaving
12,508, 17,130, and 20,734 eligible detected genes. D0 and D1 group ensembles were stable (mean
medoid ARI 0.978 and 0.971); D2-prime was less stable (0.657), which independently limits how
strongly a cell-type hierarchy can be interpreted there.

D0 selected `alpha_group=80`, `alpha_global=20` by the registered mean-log-loss rule. Those values
were carried unchanged into D1 and D2-prime.

## Result

| Dataset | Best non-hierarchical | Baseline top-1 | Hierarchical top-1 | Baseline NLL | Hierarchical NLL | Group-only NLL |
|---|---:|---:|---:|---:|---:|---:|
| D0 development | pooled | 98.182% | 98.233% | 0.05973 | 0.08028 | 0.06540 |
| D1 confirmation | pooled | 98.296% | 98.351% | 0.05381 | 0.07612 | 0.05417 |
| D2-prime confirmation | pooled | 88.554% | 88.008% | 0.27846 | 0.43126 | 0.30011 |

Combined confirmation NLL was **45.4% worse**, rather than at least 5% better. D2-prime lost
0.546 top-1 percentage point and worsened NLL by 54.9%. Real groups were modestly better than
size-preserving shuffled groups, so the labels contain assignment signal, but the registered
additive shrinkage converts that signal into worse tail probabilities. Group-only sharing was
also no better than pooled by NLL on any dataset.

The implementation checks passed: target denominators were identical across real/shuffled/
collapse controls; collapse matched pooled exactly in all aggregate metrics; D1 median wall time
was 13.83 s (1.375x the prior packed evaluator); and D1 peak RSS was 4,548,216 KB (~4.34 GiB),
below the 6 GiB limit.

## Decision

Do not make hierarchy the production default and do not use it as a positive paper claim. The
negative result supports a cleaner conclusion: for this masked recovery task, whole-sample pooled
sharing is difficult to improve and is better calibrated than the tested additive hierarchy.
The optional group evaluator is worth retaining because it provides a compact test bed for future
models and its all-one collapse is a useful invariant.

The registered grid cannot support a claim that hierarchical sharing fails. A follow-up must
parameterize effective candidate-level mixture strength, use cross-fitting, and validate on data
not inspected in either the registered gate or the post-hoc scale audit. Candidate-normalized
convex partial pooling is preferable to another raw-alpha grid because its weights have the same
meaning for common and rare candidate sets.
