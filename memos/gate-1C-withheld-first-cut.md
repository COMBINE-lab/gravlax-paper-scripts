# Gate -1C, first cut — recovery of withheld annotation features

**Date:** 2026-08-28 · **Scale:** full D0, 1,225 cells
**Panel:** `annotations/withheld/w1.{gtf,manifest.tsv}` — 900 real, expressed features held out of
GENCODE v49 by `scripts/70_withhold_annotation.py`: 300 whole genes, 300 isoforms of multi-isoform
genes, 300 genes with every transcript's 3'-terminal exon trimmed 500 bp. Expression-stratified
(terciles of the oracle's gene sums, minimum 20 UMIs), deterministic under seed 20260828.

This panel supersedes the source plan's synthetic A3 edits, per the held-out-truth design the user
proposed (oarfish genome-projection / TranSigner style): withheld features are real biology with
known oracle counts, and the withholding fraction is a dial.

## Result — measured from already-existing matrices, no new compute

Replay(A1) vs fresh oracle(A1), restricted to the withheld features. The archive is annotation-free,
so these are features it had never seen in any form at ingest.

| Class | n | Oracle UMIs | Replay UMIs | Mass recall | Genes within 10% | Median rel. error |
|---|---:|---:|---:|---:|---:|---:|
| whole gene | 300 | 170,938 | 170,737 | **99.88%** | **300/300** | 0.0000 |
| isoform | 300 | 90,130 | 90,082 | 99.95% | 298/300 | 0.0000 |
| 3' trim | 300 | 213,887 | 213,881 | 100.00% | 298/300 | 0.0000 |

A median relative error of zero means the majority of withheld genes are recovered **exactly**,
UMI for UMI. The whole-gene class is the headline: these are genes that a lab whose annotation was
`w1` would have **no trace of in their count matrix** — 170,938 UMIs across 300 genes, unreachable
by any matrix-side method, recovered at 99.88% by replaying the archive with the newer annotation.

## What this does and does not show

- It **does** show unseen-feature recovery for a genuinely annotation-free archive: nothing about
  these features was available at ingest, by construction of the whole pipeline.
- It does **not** yet include the matrix-only negative control run concretely (fresh STARsolo under
  `w1`, showing where the withheld genes' reads end up — unassigned or absorbed by neighbours —
  and that no matrix operation recovers them). That run also enables the sharper question for the
  isoform and 3'-trim classes, whose gene-level counts barely move by design: their signal lives in
  transcript-level / spliced-unspliced semantics, which the replay engine does not yet emit.
- Fidelity on withheld features matches overall fidelity (Gate -1A, 0.217%) — recovery is not
  degraded for unseen features, which is the point the doc's Gate -1C exists to test.

## Next for the campaign

1. Fresh STARsolo oracle under `w1` (one run) → the matrix-only baseline and the misassignment map.
2. Replay under `w1` as well, diffed against that oracle — cross-annotation replay (Gate -1B shape)
   with labeled ground truth on exactly the features that differ.
3. Sweep the withholding fraction and repeat under a second seed.
4. Extend replay to `GeneFull`/Velocyto semantics so the isoform and 3'-trim classes get a
   measurement that can actually move.
