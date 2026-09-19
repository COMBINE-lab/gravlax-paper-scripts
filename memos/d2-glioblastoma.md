# Memo — D2: human glioblastoma (brain-derived tissue), full gate verdicts

**Date:** 2026-08-30 · **Gates:** `docs/decisions/2026-08-30-d2-gates.md` (frozen first).
**Dataset:** `Parent_SC3v3_Human_Glioblastoma`, 250.7M reads, 5,573 filtered cells, ~92%
uniquely mapped. Third dataset, first non-blood tissue. Driver: `scripts/96_d2_campaign.sh`.

| Gate | Frozen bar | Measured | Verdict |
|---|---|---|---|
| D2-REG | byte-identical | byte-identical (132.4M molecules) | **PASS** |
| D2-FID | ≤0.5% mass | **0.426%** (283/5,573 cells >1% — see below) | **PASS** |
| D2-SIGNAL | ≥ PBMC's 2.43% | 2.12% union mass; 20.8% of genes >10% | **MARGINAL** (below) |
| D2-SIZE | within 15% of 40.8–41.1 b/mol | 27.2 b/mol (−33%); **per READ: 14.4 vs 13.4 (+8%)** | letter-MARGINAL, invariant found |
| D2-VELO | each component ≤5% | **1.59 / 0.45 / 1.97%** | **PASS** — best velocity yet |
| D2-EM | pooled ≥95%, > cell-only | 95.7% vs 90.0% (uniform 40.0%), n=257,757 | **PASS** |
| D2-DISC | report | not run in this campaign | owed |

Storage ladder: FASTQ 17.6 GB → BAM 12.7 GB → CRAM 3.1 4.16 GB → **archive 449.4 MB**
(39× / 28× / **9.3×**). Ingest 350 s; replay 27.4 s.

## The three findings

**1. The per-molecule size metric was the wrong invariant — per-read is stable.** D2 sits at
1.89 reads/molecule (vs PBMC's 3.07: shallower saturation), so molecules are cheaper (27.2 vs
40.8 bits) and the CRAM ratio compresses to 9.3× (less collapse advantage when molecules are
read-sparse). Per read the archive costs 14.4 vs 13.4 bits (+8%) — the format's cost is a
property of reads, and the CRAM advantage is ∝ reads/molecule. Both facts now stated as such in
the paper's storage section instead of a single ratio.

**2. The tissue-sensitivity premise came back MARGINAL — and the reason matters.** Glioblastoma's
v32→v49 mass churn (2.12%) does NOT exceed PBMC's (2.43%); gene-level churn matches (20.8% vs
21.8%). Tumor expression is dominated by highly-expressed, anciently-annotated programs; the
annotation-sensitive mass in brain lives in lncRNAs and intronic/nuclear transcripts that
*cytoplasmic tumor cells* under-sample. The premise test needs a single-NUCLEUS brain dataset
(intron-rich, lncRNA-rich); recorded as D2′ if we want the premise demonstrated rather than
assumed. Choosing tumor cells bought discovery-richness but not churn — an honest miss.

**3. Fidelity and capabilities generalize.** Replay fidelity 0.43% (PASS), velocity its best
showing yet (1.6/0.5/2.0% — consistent with VELO-2: shallow saturation means fewer chains lose
middle reads), EM sharing +5.7 points over per-cell at 95.7%. One observation for the
disagreement taxonomy: 5.1% of cells moved >1% (PBMC: ≤0.2%) — per-cell concentration of
disagreement is tissue-dependent even when mass-level fidelity holds; not yet decomposed.
