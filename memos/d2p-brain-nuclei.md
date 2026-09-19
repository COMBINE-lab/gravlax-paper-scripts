# Memo — D2′: human brain nuclei (Brain_3p), full gate verdicts

**Date:** 2026-08-31 · **Gates:** D2 table + D2′ addendum (frozen before measurement).
**Dataset:** `Brain_3p`, 263.4M reads, 6,460 nuclei, 90.9% uniquely mapped. Fourth dataset;
first single-nucleus; the annotation-sensitivity premise test.

| Gate | Frozen bar | Measured | Verdict |
|---|---|---|---|
| D2′-SIGNAL | churn > PBMC's 2.43% | **4.90% union (3.33% shared); 23.7% of genes >10%** | **PASS — premise demonstrated, 2× blood** |
| D2′-REG | byte-identical | byte-identical (130.3M molecules) | **PASS** |
| D2′-FID | ≤0.5% | 0.454% (165/6,460 cells >1%) | **PASS** |
| D2′-VELO | each component ≤5% | 2.79 / 0.57 / 3.95% | **PASS** (the ≤1.5% PBMC-era label prints MARGINAL; the frozen D2/D2′ bar is ≤5%) |
| D2′-EM | pooled ≥95 / ≥85 | pooled **90.9%** vs cell-only **59.4%** vs uniform 33.5 (n=99,852) | **MARGINAL as frozen — with the campaign's biggest finding (below)** |
| D2′-SIZE | within 15% of PBMC b/mol | 38.6 b/mol; **19.1 bits/read** | outside band — complexity-dependent, audited below |

Storage ladder: FASTQ 19.0 GB → BAM 16.4 GB → CRAM 3.1 5.36 GB → **archive 629.0 MB**
(30.2× / 26.1× / **8.5×**). Ingest 396 s; replay 28.5 s. Velocity inversion confirmed: nuclei
count 48.1M unspliced vs 13.7M spliced UMIs — and unspliced is the replay's best component.

## The three findings

**1. The premise the archive is built on is now demonstrated, not argued.** Brain-nuclei
annotation churn (v32→v49) is 4.90% of UMI mass — double PBMC's 2.43% — and Gene counting
assigns only 44.3M of 130.3M molecules (34%): two-thirds of nuclei evidence is invisible to the
standard matrix and fully retained by the archive. The two numbers together are the paper's
motivation section.

**2. Cross-cell EM sharing matters most exactly where cells are sparsest.** Per-cell EM
*collapses* on nuclei (59.4% — shallow per-cell priors), while pooled sharing holds 90.9%:
a +31.5-point gap, versus +4–6 on PBMC/tumor. MARGINAL against the frozen 95% bar, and the
strongest evidence yet for the sharing design; the bar was calibrated on deep cells and the gap
is the story.

**3. Per-read cost tracks transcriptomic complexity.** 19.1 bits/read (vs 11.6–14.3 elsewhere):
94% of molecules are fresh singleton classes, doubling shapes (3.3M) and patterns (4.1M), and
cell-of-class becomes 37% of the file. Audit findings (separate probe): only 7.7% of classes are
ambient (aggregation tier ≈ −21 MB only — the cost is *real* nuclei bookkeeping ≈ its 13.1
bit/class entropy); but the **delta codec runs 15% above the iid bound here** (230.9 vs 200.8 MB
— locality vanishes at 6,460 diluted cells), inverting the D0 delta-vs-frequency verdict.
Identified levers, deferred until after this campaign: per-block coc codec choice (delta vs
frequency-rANS, ≈ −25–30 MB), transcript-relative second representatives (rep.pos at 8.5 b/v,
≈ −10–15 MB). Pattern restructuring records its fourth negative (2.05× sharing).
