# Hierarchical EM: post-hoc scale diagnosis

**Date:** 2026-08-31  
**Status:** exploratory; not independent confirmation  
**Compact result:** `results/post-v1-hierarchical-em-scale-audit.json`

## Verdict

The registered STOP was not an implementation failure and is not evidence that hierarchy fails.
It was principally a parameter-scale failure. At registered `alpha_group=80` and
`alpha_global=20`, the D0 prior contribution across candidate/group pairs had median
`3.60e-5`, 90th percentile `0.00571`, and 99th percentile `0.2516`. A single local count therefore
dominated more than 99% of pairs; the grid never exercised the intended strong partial pooling.

## Exploratory result

D0 selected `12,800/12,800` after expanding the scale. All values below are target-weighted over
seeds 7, 17, and 29.

| Dataset | Pooled NLL | Calibrated cell/global NLL | Real hierarchy NLL | Shuffled hierarchy NLL |
|---|---:|---:|---:|---:|
| D0 | 0.05973 | 0.05768 | 0.05740 | 0.05767 |
| D1 | 0.05381 | 0.05225 | 0.05150 | 0.05194 |
| D2-prime | 0.27845 | 0.27677 | 0.27573 | 0.27709 |

Across D1+D2-prime, pooled NLL is `0.0706295`, calibrated cell/global is `0.0690617`, and the real
hierarchy is `0.0682954`. Thus stronger cell/global calibration supplies a 2.22% improvement over
pooled; adding real groups brings the total to 3.30%, a 1.11% increment beyond the calibrated
no-group model. Real labels also beat size-preserving shuffled labels on every dataset. There is
a small, consistent group signal, but most of the rescue is better calibration rather than hard
group structure.

These values are post hoc because D1 and D2-prime had already been inspected. They change the
failure diagnosis, not the registered verdict, and are not yet a positive paper claim.

## Better next model

Replace whole-transcriptome pseudo-count addition with candidate-normalized convex partial
pooling. For each target candidate set `C`, normalize cell, group, and sample abundances within
`C`, then learn simplex-constrained mixture weights using cross-fitted negative log loss. Shrink
the group distribution toward the sample distribution with leave-one-cell-out counts; allow the
model to return exactly to the calibrated no-group arm. This makes pooling strength comparable
between rare and abundant candidate sets and prevents the five-digit-alpha pathology.

Only if a hard-group candidate-normalized model adds repeatable held-out value should we invest in
a soft k-nearest-neighbor prior or a full hierarchical Dirichlet model. Confirmation must use a
new dataset or independently generated ambiguity truth, because D0, D1, and D2-prime have now all
influenced model choice.
