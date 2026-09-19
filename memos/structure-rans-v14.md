# Memo — meta-color-style structural probes, and the rANS stage (v1.4)

**Date:** 2026-08-30 · **Regression:** byte-identical, both datasets, all paths. Probes:
`scratchpad/{metacolor_probes,family_probe}.py` on the current D0 archive.

## Structural probes (the Mac-dBG partial/meta-color analogy, measured)

The meta-colored dBG idea factors recurring color *sets* through a two-level indirection:
partition the id space, extract shared partial sets, spell each set as meta-references. Three
candidate mappings onto our objects, all now measured:

1. **Paralog patterns as color sets — the closest analog. Sub-set sharing is weak:** 20.4% of
   dictionary alts are shareable by predecessor-LCP front-coding (≈5 MB at D1); alt bigram
   sharing is only 1.19×. Patterns do not share large sub-lists the way pangenome color sets do.
2. **Family site-set factoring — NEGATIVE, and diagnostic.** If patterns of one paralog family
   were all difference-sets of a single absolute site set, storing the family once would collapse
   the dictionary. Measured: 4.51M mm instances use 773k distinct relative patterns but
   **998k distinct absolute site-sets (0.77 patterns per set)** — anchors shift within families,
   so absolute sets *fragment*. Conclusion: the anchor-relative pattern interning we already do
   IS the right quotient; the meta-color move inverts here. (Third and final pattern-restructure
   negative: alt-vocab ×2, now family sets.)
3. **Cell-of-class as locus-conditioned colors ("expression modules") — weak at chunk
   granularity:** Σ H(cell | chunk) = 94% of global H(cell) (25.2 vs 26.7 MB at D0). The module
   structure of the 1,225 real cells is drowned by 324k ambient barcodes whose draw is
   locus-independent. ~6% headroom, not a design driver.

**The honest generalization:** the archive's existing quotients — UMI classes over values,
anchor-relative patterns over placements, junction chains over shapes, cell-of-class over
molecules — already are its "meta colors." Each further two-level indirection we can construct
has now been priced, and none pays more than single-digit MB. This is a claimable result, not a
disappointment: the representation is measured to sit at its own structural quotient.

## The rANS stage (v1.4) — coder wins taken, and what they proved

Per the D26 audit, zstd sat above the order-0 value bound on exactly five streams (class, weight,
rep.pos, mm.pos, mm.weight). v1.4 codes those with a static-table byte-wise rANS (direct symbols
< 128 + bit-length escape classes with raw extra bits; 12-bit frequencies; **global tables in a
`rans.tables` dictionary section, so chunk decode stays fully self-contained** — the cell-scoped
backref recoding that would beat it by another ~40% is rejected for breaking chunk independence,
per explicit decision). Context-structured streams (anchor, layout, shapes, patterns) stay
varint+zstd where zstd beats memoryless.

Measured: **D0 114.07 → 112.89 MB (−1.0%, 12.44× vs CRAM); D1 554.51 → 554.41 MB (−0.02%,
11.76×).** Regression byte-identical; queries identical (open still 0.00–0.01 s); D1 ingest
609.7 → 443.8 s (MT-BGZF extraction landing in the ingest path). The D1 null is itself the
finding: the bit-length approximation and incompressible escape bits consume the theoretical gap
at D1's value distributions. **The "entropy layer is done" audit claim is now confirmed
in-format, not just on paper.** Kept: zero cost, small D0 win, and the guard/table machinery is
the vehicle for any future per-stream model.

## Standing size verdict

D0 112.89 MB / 41.1 → 40.8 bits per molecule; D1 554.41 MB. Remaining measured headroom anywhere
in the format: ≈1–2% total. The size axis is closed as an engineering frontier; further gains are
modeling research (cell-locality models for coc, chain-conditioned shape streams) with measured
upper bounds in the single digits. Framing (near-optimality claim, fidelity ladder, Malva
capability-per-byte, cross-sample dictionary amortization) is where size wins from here.
