# Frozen gates — capability campaign: velocity, discover→replay, APA-diff, federation

**Date:** 2026-08-30 · **Status: FROZEN before implementation.** Order per direction: velocity,
then discover→replay, then APA differential, with federation investigation alongside.

## VELO — RNA-velocity semantics replay

Port of STARsolo Velocyto counting, read from source (`SoloFeature_countVelocyto.cpp`,
`Transcriptome_classifyAlign.cpp::alignToTranscriptMinOverlap`, minOverlapMinusOne=6): unique
mappers only; per (cell, UMI) — cell-scoped, so same-cell same-UMI molecules at different loci
intersect their transcript sets (usually to empty, dropping the UMI); transcript-wise
intersection across reads with type-bit OR; span bits imply exonic+intronic; per-UMI model rules
(exonModel/intronModel/mixedModel/spanModel) → spliced / unspliced / ambiguous; multigene UMIs
dropped. Our evidence substitutes the two span-extreme representatives for "all reads" — the
2-rep sufficiency frontier was measured for Gene (D19: 0.248%); velocity re-measures it under a
harder statistic, and that frontier number is itself a contribution.

| Gate | Test | PASS | MARGINAL | STOP |
|---|---|---|---|---|
| **VELO-1** | Replay S/U/A matrices from `d0.aie` (v49) vs fresh STARsolo Velocyto oracle, filtered cells | per-component mass-weighted rel. L1 ≤1.5% and per-component total mass within 5% | ≤5% | >10% |
| **VELO-2** | If not PASS: decompose disagreement (extreme-vs-all-reads, intersection semantics, UMI collapse) | reported | — | — |

## DISCR — the discover→refine→replay loop

`aie query discover --emit-gtf` writes candidate loci as single-exon transcripts; replay under
that GTF quantifies them like any annotation. Labeled test: w1 candidates; on the 45
cleanly-withheld genes (no residual-w1 overlap, oracle ≥50 UMIs):

| Gate | Test | PASS | STOP |
|---|---|---|---|
| **DISCR-1** | recall through the full loop (candidate → GTF → replay matrix) | 45/45 loci quantified ≥1 UMI | any loss |
| **DISCR-2** | replay counts of matched loci vs GeneFull oracle | median rel err ≤35% (single-exon models lose junction-spanning molecules — the gap is reported, not hidden) | >60% |

## APA-D — differential 3′ usage between cell groups

`apa --groups <tsv>` (cell → group): per-site per-group UMI counts + a per-gene usage-shift
statistic. Demonstration gate: on D0 with marker-based groups (e.g., monocyte vs T-cell barcode
sets from the oracle matrix), report top shifted sites; **negative control frozen in**: random
within-group halves must show a shift distribution centered on zero, and the observed
between-group shifts must exceed the 95th percentile of the within-group null for reported sites.

## FED — cross-archive federation (the atlas story)

`aie federate <archive...> junction|region <locus>`: one query over N archives, per-sample
per-cell results. D0+D1 is the miniature atlas.

| Gate | Test | PASS |
|---|---|---|
| **FED-1 consistency** | federated per-archive results ≡ single-archive query outputs | identical counts |
| **FED-2 latency** | federated 2-archive junction/region query wall | ≤2× single-archive |
| **FED-3 amortization (measurement, not gate)** | overlap of junction catalogues, shape/chain dictionaries across D0/D1 | reported — the shared-dictionary atlas argument with numbers |

## Deferred with a written design (not gated yet)

Cross-cell information sharing for quantification (EM-style): the archive's locus-level molecule
evidence across all cells enables (a) archive-driven equivalence classes for an alevin-fry-style
EM, and (b) resolving the multi-gene UMI mass STARsolo's MultiGeneUMI filter currently discards,
using cross-cell locus evidence as the prior. First measurement before any design freeze: how
much matrix mass the MultiGeneUMI filter discards (the upside bound).
