# Gate -1A — end-to-end annotation replay vs fresh STARsolo

**Date:** 2026-08-28 · **Host:** analysis-host · **Scale:** full D0 (66,601,887 reads; 1,225 cells)
**Code:** `anno` crate (GTF compiler + `alignToTranscript` port), `aie replay`, `aie assign-diff`
**Results:** `results/gate1a-full-v49-pseudo.json` (headline), `results/gate1a-full-v49.json`
(simple-barcode variant), `results/assigndiff-dev.json`

## Headline

**The capability exists.** One annotation-free ingest, plus a GTF supplied at query time, reproduces
fresh STARsolo processing to within **0.217% of UMI mass** — with **zero of 1,225 cells changing by
more than 1%** and the total within 0.018% — in **83 seconds** against roughly 16 minutes of fresh
processing. This is the first end-to-end demonstration; every earlier result was component-level.

| | oracle | replay | delta |
|---|---:|---:|---:|
| UMI total (1,225 cells) | 9,440,606 | 9,438,863 | **−0.018%** |
| UMI mass moved (union) | — | **0.217%** | 0.208% reassigned; net −0.018% |
| Expressed genes changing >10% | — | 3.94% | |
| Cells changing >1% of total | — | **0 / 1,225** | |
| Replay wall time | ~16 min (align+count) | **82.6 s** (23.7 s of it the barcode-frequency pass) | ~12× |

Against the operating-regime annotation signal (v32→v49, 2.43% of mass): **an 11.2× margin,
demonstrated end-to-end** rather than budgeted from components. An earlier iteration with the
simple exact-or-unique-1MM barcode rule scored 0.552% with a −0.72% systematic under-count
(`results/gate1a-full-v49.json`); implementing STARsolo's `1MM_multi_Nbase_pseudocounts`
(frequency-posterior with pseudocounts, quality-aware, N-tolerant — still purely annotation-free)
removed the under-count entirely and is now the default path.

## The discovery that made it work: multimappers carry countable evidence

The first replay was 3.5% under the oracle with no explanation left — per-read gene verdicts agreed
*exactly* with the oracle's `GX` tags (`aie assign-diff`: zero disagreement classes other than
reads both sides discard). The resolution came from entry-level matrix-vs-tag diffing: the oracle's
matrix **strictly exceeds** every tag-derived count, by 37,590 UMIs (2.2%) at dev scale, with zero
entries below — concentrated on multimapping genes.

**STARsolo's `Gene` feature counts unique-GENE reads, not unique-alignment reads.** A read mapping
to several loci is counted when the concordant genes of all its alignments reduce to one gene; the
primary record's `GX` stays `-` and the evidence lives on secondary records. Verified on EEF1A1
(decision D16): one cell, matrix count 102, primary-tag count 4, 118 secondary alignments holding
the difference. At full scale, 1,332,588 reads are counted through this rule.

Two consequences beyond the fix itself:

1. **Alternative placements are quantification-relevant evidence, not optional completeness.** The
   plan's §14.4/E4 "all alternative placements" field carries 2.2% of matrix mass. The archive
   format must keep alternative placements for multimapped molecules; the storage increment is a
   new, measurable line item for the format design.
2. **Earlier baselines were computed against a slightly wrong target.** Gate G-D's "G1" and the
   graph experiment's reference were distinct-`UB`-per-(cell,gene) from tags — 2.2% under the true
   matrix. Their relative conclusions (policy composition, graph sufficiency) stand; their absolute
   numbers must be re-based before appearing in any paper table.

## The error decomposition (dev-scale ladder, identical pipeline)

| Rung | Isolates | UMI mass moved |
|---|---|---:|
| L1: oracle alignments + oracle barcodes | assignment + collapse port | **0.032%** (0 cells off >1%) |
| L2: + our barcode correction | barcode step | 0.715% |
| L3: + annotation-free alignment | alignment ceiling | 0.869% |

The port itself is essentially exact. Barcode correction was the largest component (+0.68pp) until
STARsolo's `1MM_multi_Nbase_pseudocounts` was implemented (`build::BcCorrector`) — dev-scale
end-to-end dropped from 0.869% to **0.191%** (one cell of 1,135 above 1%), and full-scale from
0.552% to 0.217%. What remains is symmetric reassignment consistent with the G-C alignment
ceiling, which is exactly what the Track B junction catalogue exists to raise.

## What replay actually runs

`aie replay <ingest.bam> --gtf <any GTF>`: GTF compiled to transcript exon models + a per-chromosome
stabbing index (`anno`, 1.8 s for v49's 507k transcripts); each read classified by the
`alignToTranscript` port (blocks/junctions/strand only — decision D4 holds); multimapper union rule
(D16); UMI collapse with tie-merging 1MM + `MultiGeneUMI_CR` (validated in
`memos/umi-adjacency-graph.md`); matrix emitted in STARsolo raw shape so
`scripts/30_gateb_compare.py` diffs replay vs oracle unchanged.

## Honest caveats

- Molecules are currently rebuilt from the ingest **BAM** in memory; the compact archive file
  format does not exist yet. Fidelity is unaffected (same evidence), but the 55 s replay time reads
  a 4.2 GB BAM — the real archive should be faster still, and Gate F's claim should be made only
  from the real format.
- Read-level replay is what ran here. Molecule-level replay (the archive's economics) must
  approximate the per-read union rule from per-molecule evidence — the representative-placement
  question returns as a format-design constraint, now with the multimapper rule in scope.
- Matrix-level metrics only; the source doc's molecule-by-molecule agreement and disagreement
  taxonomy are the next measurement, and A2/A3 (cross-annotation, withheld-annotation) replay —
  the actual campaign — has not run yet. This memo is the -1A baseline that unblocks it.
- One dataset, STARsolo semantics, gene-level `Gene` feature only.
