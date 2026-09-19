# Memo — capability campaign: velocity, discover→replay, APA-diff, federation

**Date:** 2026-08-30 · **Gates:** `docs/decisions/2026-08-30-capability-gates.md` (frozen first).
Scripts: `91_velo_gate.py`, `92_discr_gate.py`, `93_apa_diff.py`. All full-scale.

## VELO — RNA-velocity replay: MARGINAL, with a four-iteration taxonomy worth publishing

Ported from source: `alignToTranscriptMinOverlap` (velocyto mode, TOL=6, giant-intron rule,
implicit junction concordance; 8 new unit tests) + `countVelocyto`'s cell-scoped UMI
intersection-with-bit-OR and model rules. `replay-rows --velocity` emits
spliced/unspliced/ambiguous in STARsolo's layout. Each iteration's failure was diagnosed against
the oracle and is the paper's taxonomy:

| iter | change | S / U / A rel-L1 |
|---|---|---|
| 1 | no UMI correction | 7.2 / 4.3 / 10.0% — surplus mass: uncorrected 1MM duplicates |
| 2 | 1MM merge over ALL classes | 3.3 / **10.0** / 4.5% — over-merged: STAR never corrected intronic-only UMIs |
| 3 | merge = exactly Gene's canon map (`gene_canon_map`) | 3.4 / 4.3 / 4.4% — uniform deficit |
| 4 | + empty transcript sets are TRANSPARENT (a no-transcript read is absent from STAR's stream, it does not blank the intersection) | **1.81 / 0.52 (PASS) / 3.06%** |

D1 confirms at 4.9× scale: 2.67 / 0.86 (PASS) / 6.13%. **Verdict: MARGINAL as frozen** (PASS bar
1.5%). The residual is spliced↔ambiguous churn — the signature of the two-representative-for-
all-reads substitution (middle reads constrain STAR's intersection; ours is looser). That is the
measured price of 2-rep evidence under velocity semantics — the D19 frontier re-measured under a
harder statistic, stated as such. VELO-2 stratification by chain weight (molecules with ≤2 reads
have no missing middles) is the designed follow-up if we want the attribution airtight.

## DISCR — the discover→refine→replay loop: PASS, and the loop IMPROVES quantification

`discover --emit-gtf` writes candidates as single-exon genes (AIENOVEL ids); `replay-rows` under
that GTF closes the loop. On the 45 cleanly-withheld genes: **DISCR-1 45/45 quantified; DISCR-2
median rel err 7.6% vs GeneFull** — versus 26.3% for raw discover cluster counts. Replaying
candidates applies real assignment + 1MM collapse + best-gene filtering, so the loop is not just
closed, it is the *better quantifier* of novel loci. (The feared single-exon junction-molecule
loss did not materialize at gate scale.)

## APA-D — differential 3′ usage between cell populations: demonstrated with a frozen null

`apa --groups <tsv>` emits per-site per-group counts. Demo: T cells (466, CD3E+LYZ−) vs
monocytes (365, LYZ+CD3E−), labels from the oracle matrix (annotation labels cells; the 3′
evidence itself is annotation-free), 239 top-expressed genes measurable in both runs. Real median
shift 30 bp vs null (random T-cell halves) median 14 bp; **35/239 genes (15%) exceed the null
P95 (5% expected)** — and the top hits are the right biology: PTPRC (CD45), CD44, LCP1, PSAP,
ANXA1. Multi-kb shifts read as differential terminal-exon/3′-isoform usage, not only tandem UTR
APA — presented as "differential 3′-end usage".

## FED — cross-archive federation: the miniature atlas works

`aie federate d0.aie d1.aie chr16:89562391-89562883`: per-sample per-cell counts in **0.66 s**.
FED-1 PASS (per-archive numbers byte-equal to single-archive queries — shared code path,
verified: 60,156/4,511 and 315,157/29,789). FED-2 PASS (0.66 s < 2× single). FED-3 measured:
**44.2% of D0's 170k junctions recur in D1** (which holds 672k at its depth) — the shared-core
dictionary amortization argument for archive collections, with numbers. The atlas story: N
archives, one query surface, molecule resolution — Snaptron-at-molecule-level across samples.

## Deferred by design: cross-cell information sharing for quantification

Two concrete routes, recorded for the design phase: (a) archive-driven equivalence classes
feeding an alevin-fry-style EM — the archive already stores, per multimapper molecule, the exact
alternative-placement pattern, i.e. the equivalence class, across every cell at once; (b)
recovering the multi-gene UMI mass STARsolo's MultiGeneUMI filter discards, using cross-cell
locus evidence as the prior. Pre-measurement still owed (frozen in the gates doc): instrument the
discarded multi-gene mass to bound the upside before any design freeze. The discover→replay loop
(above) is the first working instance of archive-enabled information sharing.

## Addendum — the EM upside bound, and VELO-2 attribution (measured)

**EM upside (`replay-rows --audit-multigene`, D0 v49).** Read-weighted rows: 40.24% match no
gene (the annotation-independence thesis in one number — the archive keeps all of it), 55.94%
exactly one gene (counted), **3.82% two or more genes (discarded today)**. At UMI level:
10.15M classes counted single-gene, 0.91M counted-at-best-gene with multi-gene evidence ignored,
and **691,030 classes are gene-informative yet fully discarded — a +6.25% pool of countable
molecules** an equivalence-class EM could fractionally assign (plus refinement of the 0.91M
mixed). That is a larger delta than most quantifier comparisons in this space; the EM design is
funded by this number.

**VELO-2 (entry-level completeness strata, `entries.complete.tsv`).** The two-representative
hypothesis is confirmed for ONE component and refuted for the rest:

| component | complete entries | incomplete entries |
|---|---:|---:|
| spliced | 1.41% | 1.78% |
| unspliced | 0.84% | 0.27% |
| ambiguous | **1.16%** | **3.25%** |

Missing middle reads drive the *ambiguous* disagreement (2.8×) — middles constrain the
transcript intersection that decides mixed-model status — but the spliced residual (~1.4–1.8%)
is completeness-independent, so it is NOT the 2-rep payload: the remaining suspects are the
class-level (vs STAR's read-level) UMI-correction approximation and min-overlap classifier edge
cases. Caveats recorded: "complete" entries skew small (10% of mass; an entry is complete only if
every counted class is), and the inverted unspliced stratum is consistent with single-read
intronic molecules dominating there. Net: D19's two-representative payload survives velocity
semantics better than feared — its cost is concentrated in the ambiguous component.
