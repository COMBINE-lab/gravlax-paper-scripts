# Track B — the cross-sample junction catalogue does NOT raise the alignment ceiling

**Date:** 2026-08-28 · **Scale:** full D0 vs sjdb_v49 reference, gene-assigned reads (34,628,655)
**Verdict: negative result. Deprioritize the catalogue; the ceiling's matrix-level footprint is
already small and this lever moves it the wrong way.**

## Setup

Hypothesis (from `memos/gate-C-decision.md`): seeding the annotation-free ingest's second pass with
junctions pooled across samples — intropolis/Snaptron's move, via `--sjdbFileChrStartEnd`, no GTF —
rescues junctions weakly supported within one sample and closes the G-C gap toward 99%.

Catalogue: `scripts/60_junction_catalogue.sh`, pooled unique-read support ≥ 3, annotation-free
sources only. D0's own runs give 94,074 junctions; adding D1 (5k PBMC v3, 385M reads, aligned
annotation-free) gives **185,411, of which 91,337 are new** — the cross-sample increment existed.

## Result

| Ingest configuration | Identical evidence vs sjdb_v49 | junction divergence |
|---|---:|---:|
| plain 2-pass (baseline) | **98.6965%** | 0.8083% |
| + D0-only catalogue (control) | 98.6770% | — |
| + D0+D1 catalogue | **98.2067%** | 1.2105% |

`results/gatec-catd0d1-v49-geneassigned.json`. The control was flat (own junctions were already in
the 2-pass); the cross-sample catalogue is **half a point worse**, with the damage concentrated in
the junction category (0.81% → 1.21%) plus multiplicity (0.38% → 0.42%).

## Why, and what it teaches

The metric is agreement with an **annotation-aware** aligner. Handing the annotation-free aligner
91k additional real, data-supported junctions gives it splice placements the v49-aware reference
does not consider — including junctions that are genuinely present in the data but absent from
v49. Every read pulled onto such a junction scores as divergence. **More discovered junctions
means more disagreement with an annotation-bounded oracle, not less.** The intropolis analogy
fails in one direction: recount-style catalogues serve *quantification against a chosen
annotation later*, not *agreement with an annotation-aware alignment now*.

Two things keep this from being a loss:

1. **The fidelity budget no longer needs this lever.** Gate -1A's end-to-end residual is 0.217%
   with a symmetric profile; the alignment ceiling's matrix-level footprint (~0.15pp at dev scale)
   is a fraction of its evidence-tuple headline (1.13%), because most divergent alignments still
   land in the same gene. The catalogue was insurance we turned out not to need.
2. **The negative result sharpens a paper claim.** The archive's annotation-free alignment sits at
   a measured, stable ~98.7% agreement with annotation-aware alignment, and pushing junction
   discovery harder moves *away* from any specific annotation. That is an honest characterization
   of what "annotation-free" costs and buys — and it implies the ceiling is annotation-*richness*
   driven (see the v32-vs-v49 finding in `memos/gate-C-decision.md`), not fixable by more data.

## If it is ever revisited

Only with a support-weighted design: per-junction sample counts and expression-weighted insertion
(so a junction seen in 3 reads of one sample cannot outrank v49's canonical form of the same
intron), evaluated at matrix level rather than evidence-tuple level. Not before D2/D3 exist.
