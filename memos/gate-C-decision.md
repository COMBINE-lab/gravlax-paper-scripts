# Gate G-C decision — alignment invariance under annotation change

**Date:** 2026-08-28 · **Host:** analysis-host · **Scale:** full D0 (66,601,887 reads)
**Provenance:** `experiments/manifests/gates.yaml` · **Frozen thresholds:** `docs/decisions/2026-08-28-frozen-gates.md`

## Recommendation: **MARGINAL — continue**, with the alignment ceiling recorded as a first-class limit

No STOP trigger fires, but the frozen PASS bar is missed. This is the project's binding fidelity
constraint and must be carried into every downstream claim rather than treated as noise.

## Question

The archive is built from an annotation-free alignment; the oracle it must reproduce uses an
annotation-aware one (STAR inserts the GTF's junctions into the index). If those disagree on where
reads land, the archive has a fidelity ceiling that no amount of stored evidence can repair.

Compared: the per-read annotation-independent evidence tuple — chromosome, strand, aligned blocks,
splice junctions, `NH` — i.e. exactly the fields STAR's own gene assignment consults
(`Transcriptome_classifyAlign.cpp`), and nothing else.

## Results

Restricted to the **34,628,655 reads the oracle assigns to exactly one gene**, which is how the
frozen threshold is defined: a read no annotation claims cannot change any count.

| Comparison | Identical evidence | Verdict |
|---|---:|---|
| **(i) annotation-free 2-pass vs (ii) sjdb v49** (A1, baseline) | **98.6965%** | PASS ≥99.0%, STOP <97% → **MARGINAL** |
| **(i) annotation-free 2-pass vs (iii) sjdb v32** (A2-far, the operating regime) | **98.8682%** | **MARGINAL**, closer to the bar |

Divergence taxonomy on the same population:

| Category | Reads | Share |
|---|---:|---:|
| junction — splice structure differs | 279,906 | 0.8083% |
| multiplicity — `NH` changed | 130,118 | 0.3758% |
| locus — different chromosome or strand | 36,648 | 0.1058% |
| presence — aligned in one config only | 3,466 | 0.0100% |
| block_boundary — same splice, shifted ends | 1,250 | 0.0036% |

Against **v32** (A2-far) the same taxonomy is uniformly milder: junction 241,702 (0.6980%),
multiplicity 117,958 (0.3406%), locus 27,585 (0.0797%), presence 3,542 (0.0102%),
block_boundary 1,152 (0.0033%).

Over **all** reads rather than gene-assigned ones, and with the single-pass control:

| Comparison | Identical | junction | locus | multiplicity | presence |
|---|---:|---:|---:|---:|---:|
| annotation-free **2-pass** vs sjdb v49 | 98.3602% | 0.5485% | 0.4157% | 0.3927% | 0.2722% |
| annotation-free **1-pass** vs sjdb v49 | 93.5206% | 1.5936% | 4.1524% | 0.5439% | 0.1869% |

`results/gatec-2pass-v49-geneassigned.json`, `results/gatec-2pass-v32-geneassigned.json`,
`results/gatec-{2pass,1pass}-v49-all.json`.

## What the numbers say

**Two-pass discovery is not optional — it recovers three-quarters of the gap.** Single-pass
annotation-free alignment agrees with the annotation-aware one on 93.52% of reads; two-pass raises
that to 98.36%, closing 6.48% → 1.64%. The mechanism is visible in the taxonomy: without a second
pass, `locus` divergence dominates at 4.15% — reads spanning junctions the index does not know get
soft-clipped and placed somewhere else entirely. Two-pass discovery drops that to 0.42%. This
answers the sub-question the gate was designed around: **annotation-free two-pass alignment does
substantially, but not completely, recover what sjdb insertion provides.**

**The residual is splice junctions, and restricting to gene-assigned reads concentrates it.** Going
from all reads to gene-assigned reads, `locus` falls 0.42% → 0.11% and `presence` all but vanishes
(0.27% → 0.01%), confirming that the gene-assigned population is the clean one. But `junction`
*rises*, 0.55% → 0.81%, because spliced reads are over-represented among reads a gene claims. So
the 1.30% shortfall is not diffuse noise: it is overwhelmingly reads whose splice placement depends
on knowing the annotation in advance.

**The ceiling tracks how much annotation is inserted, not how distant it is.** This is the result I
did not expect. Agreement with the *distant* v32 (98.87%) is **better** than with the *near* v49
(98.70%), and every divergence category is smaller. The reason is that v49 carries 507,365
transcripts against v32's 227,462, so sjdb insertion hands the aligner far more junctions to pull
reads onto — junctions the annotation-free run had no way to discover. The gap is a function of
annotation *richness*, not annotation *age*.

That is convenient for this project, because Gate G-B established the operating regime as distant
annotation pairs, and the ceiling in exactly that regime is the better of the two: **1.13%, not
1.30%.** It also predicts the ceiling will worsen as GENCODE keeps growing, which is a claim the
paper should make explicitly rather than let a reviewer discover.

## Why this matters more than the 0.3 percentage points suggest

Replay can never be more faithful than its alignment. The error budget now reads:

| Source | Magnitude |
|---|---:|
| Alignment ceiling (G-C) | 1.30% vs v49; **1.13% vs v32**, the operating regime |
| Grouping error (G-D) | 0.94% of per-cell count mass |
| **Annotation signal to be captured (G-B, far pair)** | **2.43% of UMI mass** |

The two error sources are on different denominators and are not simply additive, but the signal is
of the same order as the errors, not orders of magnitude above them. **The comfortable
interpretation — that replay fidelity is a solved problem because E1 is a sufficient statistic
(decision D4) — is wrong.** Sufficiency of the evidence *given the alignment* was never the
question; the alignment itself is annotation-sensitive, and that is where the fidelity budget goes.

## Improvement path, and why it is principled rather than a patch

The junction residual is addressable without reintroducing an annotation: seed the ingest's second
pass with a **union of junctions discovered across many samples**, rather than only those found in
the sample being processed. This is exactly what intropolis and Snaptron do for bulk RNA-seq
(`docs/prior-art.md` §2) — a junction catalogue built from data, not from a gene model, so it stays
annotation-free. A rare junction unsupported within one shallow sample is well supported across
thousands. This is worth testing before Gate -1, since it directly raises the fidelity ceiling that
bounds every replay claim, and it strengthens rather than weakens the "single-cell Monorail"
framing.

## Verdict against the frozen criteria

| Criterion | Threshold | Observed | Verdict |
|---|---|---:|---|
| Identical tuples, annotation-free vs A1, gene-assigned reads | PASS ≥99.0% | 98.6965% | **MARGINAL** |
| STOP trigger | <97% | not met | **no STOP** |

Continue. Both comparisons are now complete. The binding number for the project is the v32 one,
**98.87% / 1.13% ceiling**, since that is the regime G-B showed the capability is worth having in.
