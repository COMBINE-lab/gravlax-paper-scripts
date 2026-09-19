# Memo — EM multimapper recovery with cross-cell sharing: PASS on both datasets

**Date:** 2026-08-30 · **Gates:** `docs/decisions/2026-08-30-em-gates.md` (frozen first; funded
by D33's +6.25% discarded-classes bound). Instrument: `aie em <archive> --gtf` (masked-evidence
protocol; `--mask 0` = impact report).

## The experiment

Labeled truth from real data, no simulation: 20% of MIXED classes (unique evidence pins the gene)
are masked down to their multi-gene evidence — making them indistinguishable from the 691k
classes discarded today — and each sharing mode must recover the hidden gene. Same candidate
sets, same truth, all modes.

## EM-1 — masked recovery (top-1 / expected accuracy)

| mode | D0 (n=173,109) | D1 (n=903,904) |
|---|---|---|
| uniform baseline | 39.1 / 41.6% | 35.4 / 40.7% |
| cell-only EM | 94.4 / 94.5% | 93.4 / 94.1% |
| **pooled (cross-cell)** | **98.2 / 97.1%** | **98.3 / 97.4%** |
| blend (α=20) | 98.3 / 97.3% | 98.4 / 97.7% |

**PASS as frozen on both datasets** (bar: ≥75% top-1, ≥15 pts over uniform, beats cell-only).
Cross-cell sharing cuts per-cell EM's residual error by ~70–75% (5.6→1.7% at D0, 6.6→1.6% at
D1). α-insensitive (5/20/80 identical) and seed-stable; pooled is the adopted mode by parsimony.
Truth-coverage caveat: 2.5% (D0) / 1.7% (D1) of masked classes lose their true gene from the
candidate set (the unique row was the only evidence for it) — recovery is scored on the covered
set, coverage reported.

## EM-2 — matrix impact (pooled mode, D0)

**+688,558 recovered classes = +6.28% of the counted matrix, 85.0% assigned with r > 0.8.** The
recipients are exactly where multimapping lives — ribosomal protein paralogs (RPL41, RPL17,
RPL22, RPS10/25/27…), mitochondrial genes (MT-CO1, MT-ATP6), EEF1G, H3-3A — no anomalous
recipients. Recovered counts are a separate additive layer, never mixed silently into the
byte-exact STARsolo-replay matrices (the fidelity contract stays intact; EM counts are an opt-in
product of the same archive).

## Why the archive is the right substrate (the claim)

A per-cell quantifier sees one cell's evidence; the archive holds every cell's
equivalence-class evidence in one pass — candidate sets are the stored paralog patterns, unique
support is the stored matrix, and sample-wide π needs no second pass over reads. The
discover→replay loop (D32) and this EM are the two working instances of archive-enabled
information sharing; STARsolo `--soloMultiMappers EM` replication (EM-0) stays deferred as a
separate fidelity gate.

## Addendum — EM-0: STARsolo EM replicated (PASS)

Port of `SoloFeature_collapseUMIall.cpp`'s multimapper EM (intersection candidate sets per UMI,
skip-if-seen-among-unique — which is precisely our mixed/multi-only class split — raw UMIs,
init = unique+uniform, 0.01 zeroing, max-abs-change < 0.01 convergence), unique side = the
byte-exact Gene replay. Against a fresh `--soloMultiMappers EM` oracle (filtered cells,
UniqueAndMult-EM vs UniqueAndMult-EM): **rel. L1 = 0.581% — PASS** (bar ≤1%). The argument is
now complete in both directions: the archive reproduces STARsolo's own per-cell EM (EM-0), and
the pooled cross-cell EM the archive uniquely enables beats that per-cell design by ~4 points
top-1 on labeled masked evidence (EM-1). Product artifacts: `aie em --star` (replication) and
`aie em --emit` (the additive recovered layer, 594,290 fractional entries at D0).
