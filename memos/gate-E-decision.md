# Gate G-E decision — storage

**Date:** 2026-08-28 · **Host:** analysis-host · **Scale:** full D0 (pbmc_1k_v3, 66,601,887 reads)
**Provenance:** `experiments/manifests/gates.yaml` · **Frozen thresholds:** `docs/decisions/2026-08-28-frozen-gates.md`

## Recommendation: **PASS**

This is the gate that killed the predecessor project, and it was front-loaded into week 1 for that
reason (decision D2). It clears the frozen bar in the **conservative** variant — the one that keeps
UMIs and therefore keeps the option to redo the per-gene collapse under a future annotation.

## Archive

Built by `aie build` from the annotation-free two-pass ingest alignment: whitelist barcode
correction, single-linkage locus clustering at Gate G-D's adopted 50 kb, 1-mismatch UMI collapse.
No gene model is consulted at any point.

- mapped primary reads: 63,683,421 · with usable barcode: 61,055,878 · dropped: 2,579,055 (4.2%)
- distinct barcodes: 317,891 · **molecules: 19,918,212** · **reads per molecule: 3.065**

| Variant | Size | bits/molecule | vs CRAM 3.1 (ingest) | vs CRAM 3.1 (oracle) | vs FASTQ.gz | vs Malva |
|---|---:|---:|---:|---:|---:|---:|
| E0 with UMI | 117.4 MB | 47.17 | 11.96× | 13.47× | 42.9× | 3.34× |
| E0 no UMI | 50.2 MB | 20.17 | 27.98× | 31.49× | 100.4× | 7.82× |
| **E1 with UMI** | **119.5 MB** | **48.00** | **11.76×** | **13.23×** | **42.2×** | **3.29×** |
| E1 no UMI | 52.3 MB | 21.00 | 26.87× | 30.25× | 96.4× | 7.51× |

Baselines, all measured on this dataset: CRAM 3.1 of the annotation-free ingest BAM with
`--capture-all-tags` = 1.405 GB; CRAM 3.1 of the oracle BAM = 1.581 GB; compressed FASTQ = 5.039 GB;
Malva's index = 393 MB. `results/archive-full.json`, `results/gatee-baselines-full.tsv`.

## Verdict against the frozen criteria

| Criterion | Threshold | Observed (E1, UMIs retained) | Verdict |
|---|---|---:|---|
| Archive vs CRAM 3.1 | PASS ≥10×, MARGINAL 5–10×, STOP <5× | **11.76×** | **PASS** |
| Archive vs Malva index | PASS ≤393 MB | 119.5 MB (3.29× smaller) | **PASS** |

## What the numbers say

**Splice structure is nearly free.** E1 costs **0.83 bits/molecule more than E0** (48.00 vs 47.17),
a 1.8% size increase, to carry every aligned block and splice junction. The fidelity work showed E1
is the evidence class that matters; it turns out to cost almost nothing. There is no interesting
accuracy/size trade-off along the E0→E1 axis, and the planned ablation should say so plainly rather
than present a curve implying a decision exists.

**The UMI is the archive.** It is 67.2 MB of 119.5 MB — **56% of the bytes, 27.00 of 48.00 bits per
molecule** — exactly as `docs/decisions/2026-08-28-archive-entropy-budget.md` predicted before any
measurement. Dropping it gives a further 2.29×. Whether that is legitimate is decided by Gate G-D,
not by Gate G-E: the UMI is only needed if a future annotation must redo the per-gene collapse.
Since G-D passed, freezing the collapse at ingest is defensible, and **both variants should be
reported in the paper** rather than the favourable one.

| Stream | Share | bits/molecule |
|---|---:|---:|
| UMI | 56.25% | 27.00 |
| cell barcode | 23.55% | 11.30 |
| position delta | 7.47% | 3.59 |
| read count | 6.04% | 2.90 |
| shape id | 3.74% | 1.80 |
| shape dictionary | 1.28% | 0.62 |
| flags | 0.60% | 0.29 |
| barcode dictionary | 1.06% | — |
| chromosome | ~0% | ~0 |

**The archive holds roughly twice the molecules gene-level quantification keeps.** 19.9M against
STARsolo's 10.2M gene-assigned molecules — because it retains intronic, intergenic and antisense
molecules that no current annotation claims, which is the entire point: a future annotation may
claim them. The plan flagged this as an under-appreciated cost; it is now measured at ~1.95×, and
the archive clears the bar anyway.

**The pre-registered pessimistic projection was beaten.** The frozen-thresholds document estimated
"~30M records × 12 B ≈ 360 MB, roughly Malva's 375 MB — the bar is tight and that is the expected
outcome." The actual result is 19.9M molecules at 6.0 B each = 119.5 MB, comfortably inside. Two
things account for it: barcode interning (storing a dense index rather than the 32-bit packed
barcode cut the cell stream from 10.6 MB to 4.7 MB at dev scale) and shape interning plus positional
delta coding, which together reduce geometry to 2.4 bits/molecule across 385,001 distinct shapes.

## Honest caveats

- **Not a like-for-like comparison with CRAM, and the paper must say so.** CRAM is lossless: it
  keeps sequence, quality, read names and every alignment. The archive is a lossy summary that
  retains only what annotation replay needs. The right claim is "an order of magnitude smaller than
  the barcoded BAM people are told to keep, while still answering annotation-replay queries the
  count matrix cannot" — not "a better compressor."
- **Boiler (NAR 2016) reported 39–56× versus BAM** for bulk, coverage-domain, molecule-destroying
  archives (`docs/prior-art.md` §5). Our 11.76× vs CRAM 3.1 is not a compression record and must
  not be presented as one. The distinguishing claim is molecule resolution, which Boiler discards.
- One dataset at one depth. Reads per molecule (3.065) drives the win and scales with sequencing
  saturation, so a shallower library compresses **worse**, not better. Must be re-measured on D1/D2
  with saturation reported alongside — the same discipline Gate G-D's dev-vs-full reversal forced.
- The barcode stream is 11.30 bits/molecule for 317,891 distinct barcodes (~18.3 bits raw). Sorting
  molecules by cell rather than position would likely cut this further, at some cost to the position
  stream. Not pursued: the gate already passes and the trade is not obviously favourable.

## Consequence

All five week-1 gates are now decided, none fired a STOP, and the storage question that ended the
predecessor project is answered affirmatively. The project's binding constraint is **not** storage —
it is the fidelity budget in `memos/gate-C-decision.md`, where a 1.13% alignment ceiling and a 0.94%
grouping error sit against a 2.43% annotation signal.
