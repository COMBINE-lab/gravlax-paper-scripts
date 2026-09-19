# Memo — coding-stack audit: entropy layer vs structure, and factorization candidates

**Date:** 2026-08-29 · **Instrument:** `aie debug <archive> [--dump-dir]` (new; read-only) plus
per-stream recompression and value-entropy analysis (`scratchpad/{order1,cellclass2}.py`).
**Archive:** `runs/archive/d0.aie` v1.1 (121.41 MB). All numbers full D0.

## Question 1 — is the entropy layer leaving money on the table? **No.**

Per-stream accounting (chunk streams concatenated across all 822 chunks, then zstd-19'd whole;
"order-0 bound" = entropy of the decoded varint *values*, the ideal memoryless value coder):

| stream | raw MB | zstd-19 MB | order-0 value bound MB | verdict |
|---|---:|---:|---:|---|
| cell | 35.83 | 28.80 | 27.91 | at bound — fat is the *model*, not the coder |
| rep.pos | 38.48 | 20.00 | 19.54 | at bound |
| rep.shape | 45.90 | 13.77 | 19.96 | zstd beats memoryless (locus repetition) |
| class | 28.26 | 9.09 | 8.30 | near bound |
| anchor | 23.49 | 9.02 | 10.33 | beats memoryless |
| weight | 21.60 | 7.88 | 7.37 | near bound (~0.5 MB rANS crumbs) |
| mm.pattern | 9.70 | 3.98 | 7.70 | beats memoryless |
| layout | 66.47 | 3.04 | 6.27 | beats memoryless 2× (runs of 68 distinct triples) |
| mm.pos | 4.87 | 1.87 | 1.89 | at bound |
| mm.shape | 5.88 | 1.32 | 2.15 | beats memoryless |
| mm.weight | 4.51 | 0.64 | 0.37 | ~0.3 MB rANS crumbs |

Sum of independently-compressed streams: 99.4 MB vs 101.3 MB inside the actual 4 Mb chunk frames —
**random access costs +1.9%**, cheap. On five streams zstd already codes *below* the memoryless
value bound (it exploits cross-value repetition); on the rest it sits essentially at the bound.
A rANS/arithmetic swap would harvest ≈0.5–1 MB total (<1% of file). Order-1 conditional entropies
confirm structure exists in class/rep.shape/mm.pattern, but the order-1 numbers with 0.5–1.5M
contexts over ≤33M values are overfit (model cost uncharged) and are upper bounds on hope, not
plans. **Conclusion: the entropy coder is done; every remaining megabyte is a modeling/structure
question.** (This also retroactively strengthens the D25 lz4 rejection.)

## Question 2 — factorable patterns and equivalence classes, measured

**1. Cell is a pure function of the UMI class — verified, 0 violations over 22.16M molecules.**
Classes are global (cell, value) equivalence classes (D20), so every backref molecule's cell is
*derivable from its class token*. Only fresh-class molecules (19.27M of 22.16M; 87.0%) need a cell
entry. Measured bounds: fresh-only delta stream order-0 = 25.21 MB vs current 28.80 MB zstd →
**≈3.5 MB saved, exact, zero fidelity risk** (regression-protected like everything else). The
biggest single win found, and it is pure derivation, not coding.

**2. rep.pos leading zeros: 58.0% of all 33.1M rep offsets are 0** — the molecule's first stored
representative *is* the anchor. Order-0 prices those zeros at ~1.9 MB compressed. Elide them
behind a 1-bit "first-rep-at-anchor" flag in the layout strand byte (7 bits currently unused;
layout compresses 22:1 so the flag is ~free). **≈1.5–2 MB, exact.**

**3. Class token stream**: 87.0% of tokens are fresh (value 0); zstd gets 9.09 MB vs the 8.30
memoryless bound. Fresh-run RLE + separate backref-distance stream is the experiment; the honest
expectation is 1–3 MB, not the overfit order-1 fantasy (4.77 MB).

**4. Shapes → junction-chain equivalence classes** (terminal block lengths factored out): 809,069
shapes collapse to 117,343 distinct chains (**6.89× sharing** — the equivalence classes are real),
but the factored encoding saves only ~0.5 MB raw (residual chain-id+extensions per shape eat the
gain), on a section that compresses to 3.11 MB anyway. Hop (intron-length) vocabulary inside the
chain dict: 8.46× sharing, ceiling ~1 MB. Real but small; worth bundling only if we touch the
shape dict for other reasons. The *semantic* value (chains ≈ splice isoform skeletons, a queryable
object) may exceed the storage value — noted for the paper, not for bytes.

**5. Pattern-alt vocabulary — negative result, recorded.** 2.51M alt entries share only **2.14×**
across 1.17M distinct (chrom, offset, flip) triples, and 1.14M of the offsets are unique.
A vocabulary indirection would break even at best. Not doing it.

**6. Edges (5.52 MB compressed, 6.34 raw)** barely compress — svarint deltas over cell-scoped
class pairs are near-random by construction (class ids are molecule-ordered). No factoring found;
this is what 1MM adjacency information costs. Cells dict (1.30 MB) is packed random barcodes —
incompressible, correctly so.

## Recommendation

Implement (1) cell-from-class elision and (2) rep.pos anchor-zero elision — both exact
derivations, together ≈5–5.5 MB (−4.5%, → ~116 MB, ≈12.1× vs CRAM) with byte-identity regression
as the gate. Run (3) class RLE as a measured experiment afterwards. Skip (4) unless the shape dict
is opened anyway; (5) is closed negative; leave (6). Do not swap entropy coders.

## Outcome addendum (implemented same day — D27, format v1.2)

Both elisions landed, with one instructive correction. (1) as designed: cell → global
`cellofclass` section, 28.8 → 26.0 MB. (2)'s first version (flag bit in the layout byte) was a
**wash**: rep.pos −2.2 MB but the varying flag destroyed layout's run structure (3.04 → 5.62 MB) —
a per-molecule flag costs as much as the zeros it removes; elision only pays when the flag is
*implied*. Fixed by making it structural: reps store the span-minimum first, chains are
position-sorted, and "first stored rep == anchor" became an enforced writer invariant (the bail
caught the wrong rep order on the first attempt) — elided by rule (n_chains > 0), no flag, layout
restored. Class RLE (3) measured 8.89 vs 9.09 MB and is **rejected** (recorded negative). Final:
**113.73 MB, 12.35× vs CRAM** (from 121.41 / 11.57×), regression byte-identical, query outputs
identical; open 0.07 → 0.17 s from eager cellofclass decode (lazy dictionary loading remains the
recorded lever).
