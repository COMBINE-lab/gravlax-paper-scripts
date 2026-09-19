# Gates -1B and -1C — cross-annotation and withheld-annotation replay

**Date:** 2026-08-28 · **Host:** analysis-host · **Scale:** full D0 (66.6M reads, 1,225 cells)
**Results:** `results/gate1b-full-{v32,v48,w1}.json`, `results/delta-capture-v32-v49.json`,
`results/gate1a-full-v49-pseudo.json`, `memos/gate-1C-withheld-first-cut.md`

## Recommendation: **PASS on both gates.** The source doc's hard continuation gates A–D are now all demonstrated.

## Gate -1B — one archive, four annotations, no rebuild

The same annotation-free ingest, replayed against each annotation's own fresh STARsolo oracle:

| Annotation | Year / nature | UMI mass moved | Cells changing >1% | Total UMI delta |
|---|---|---:|---:|---:|
| GENCODE v32 | 2019 (CellRanger 2020-A era) | **0.207%** | 0 / 1,225 | −0.005% |
| GENCODE v48 | 2025 | 0.212% | 0 / 1,225 | −0.005% |
| GENCODE v49 | 2026 | 0.217% | 0 / 1,225 | −0.018% |
| w1 | v49 minus 900 withheld features | 0.246% | 0 / 1,225 | −0.067% |

**Fidelity is flat across seven years of annotation drift and across a deliberately perturbed
annotation.** Each replay costs ~80–90 s against ~16 min of fresh processing per annotation — and
the fresh path also pays a per-annotation index build, which the archive never does.

### The decisive metric: capturing the change, not the consensus

The doc (§11) is explicit that overall agreement can hide failure on exactly the molecules the
capability exists for. On the **3,283 annotation-sensitive genes** (oracle count changes >10%
between v32 and v49, ≥20 UMIs):

| Metric | Value |
|---|---:|
| Sign agreement of the per-gene change | **3,281 / 3,283 (99.94%)** |
| Replayed delta within 10% of the oracle delta | 97.01% |
| Replayed delta within 25% | 98.78% |
| Total change mass, replay error | **1.17%** of 360,556 UMIs |
| (supplementary Pearson on deltas) | 0.99993 |

The archive reproduces the annotation-driven *change itself* to ~1% of the change mass. This is
the continuation-gate C requirement, met with a large margin.

## Gate -1C — withheld features, with the negative control now concrete

Recovery (from `memos/gate-1C-withheld-first-cut.md`): the 900 held-out expressed features are
recovered from the annotation-free archive at **99.88–100% mass recall**, 300/300 whole genes
within 10%, median per-gene error zero.

The matrix-only negative control, now measured with a fresh STARsolo run under `w1`:

- The w1 matrix contains **0 of the 300 withheld genes** — 170,938 UMIs of signal with no trace.
- Only **8.8%** of that mass is absorbed by neighbouring genes; **91% simply vanishes** from the
  matrix. No matrix-side computation can bring it back; replay recovers 99.88% of it.
- The 3′-trim class is not inert at gene level after all: trimming 500 bp off terminal exons costs
  those 300 genes **9.2%** of their mass under w1 (reads in the trimmed window are no longer
  exonic-concordant). Replay under full v49 recovers that mass at 100.00% — the "3′ end moved"
  case, the highest-signal class for 3′ chemistry, is demonstrated and not just constructed.
- A stale-annotation cost number for the motivation section: this modest 900-feature withholding
  costs the w1 matrix **1.88%** of its total counts (deficit 177,277 UMIs), invisibly.

## Standing of the source doc's hard continuation gates

| Gate | Requirement | Status |
|---|---|---|
| A — capability | one archive supports A1/A2/A3 without rebuilding | **demonstrated** (4 annotations) |
| B — fidelity | replay ≈ fresh on ordinary quantification | **0.21–0.25% mass, 0 cells >1%** |
| C — annotation-sensitive fidelity | changes captured accurately | **99.94% sign, 1.17% delta error** |
| D — genuinely new features | withheld features recovered with high precision/recall | **99.88% mass recall, 300/300** |
| E — storage headroom | archive ≪ raw representations | 24.06× vs CRAM 3.1 (`memos/gate-E-decision.md`) |
| F — replay speed | avoids mapping, ≫ faster | ~12× unoptimized from a BAM; real format pending |

## Caveats that stay attached to these numbers

- Matrix-level metrics; the doc's molecule-by-molecule disagreement taxonomy remains to be built,
  and will matter for the paper's fidelity-characterization spine.
- One dataset (PBMC), one chemistry (3′ v3), STARsolo `Gene` semantics. D2/D3 and the alevin-fry
  oracle family are still the generalization tests.
- Replay currently reads the ingest BAM as an archive stand-in; Gate F's headline must come from
  the real block-structured format, which is now the principal unbuilt artifact.
- The isoform class's gene-level counts move little by design; its real test needs
  transcript-level / velocity semantics in the replay engine.
