# `.aie` v0 — the archive exists

**Date:** 2026-08-29 · **Host:** analysis-host · **Spec:** `docs/format-spec.md` (v0, frozen on Phase-0)
**Code:** `rows.rs` (shared extraction/replay core), `evidence-io/format.rs` (container),
`aie ingest-archive` / `aie replay-rows` · **Data:** `runs/archive/d0.aie`, `logs/archive-full-v2.log`

## Acceptance — all three tests pass at full scale

1. **Byte-identical regression.** Archive-sourced and BAM-sourced replay share one row abstraction
   and produce identical matrices, verified with `cmp` at dev and full scale, through both the
   class-renumbering and class-scope changes.
2. **Fidelity vs fresh STARsolo: 0.244% of UMI mass, 0 of 1,225 cells changing >1%** — statistically
   identical to the 0.248% two-rep emulation. The value-free collapse (classes + edges, class-id
   tie order) costs nothing measurable.
3. **Honest Gate F, from the real artifact:** replay = **25.3 s** (0.8 s archive load + 2.0 s GTF
   compile + 19.6 s classify/collapse) against ~16 min fresh STARsolo → **~38×**, and against
   126.4 s for the same replay fed from the BAM (80.4 s of which is BAM decode + barcode pass the
   archive paid once at ingest).

## The artifact

`d0.aie` = **151.1 MB** for 37,630,563 rows (32.1 bits/row): full D0, with *complete* replay
capability — two span-extreme representatives per junction chain, exact multimapper evidence via
the paralog-pattern dictionary, UMI adjacency instead of values, pseudocount-corrected barcodes.

| | size | ratio |
|---|---:|---:|
| vs CRAM 3.1, all tags | 1.405 GB | **9.30×** |
| vs ingest BAM | 4.209 GB | 27.9× |
| vs compressed FASTQ | 5.039 GB | 33.4× |
| vs Malva index (same data) | 393 MB | 2.60× |

**Placement on the capability/size frontier, stated honestly:** the frozen Gate G-E criterion
(≥10× vs CRAM) was defined and PASSED on the E1 molecule archive (119.5 MB / 11.76× with UMIs;
58.4 MB / 24.1× with the graph). The full `.aie` buys strictly more — exact multimapper mass
(D16's 2.2%), exact re-collapse, span-extreme fidelity — and lands at 9.30×, nominally in the old
MARGINAL band. These are different objects; the paper should present the ladder with a price at
each rung, not quote one number.

Per-stream accounting (the embedded Gate-0 measurement):

| stream | zstd MB | bits/row | note |
|---|---:|---:|---|
| rows.cell | 46.5 | 9.9 | as measured in Phase-0 |
| rows.class | 46.0 | 9.8 | **headroom**: plain delta of renumbered ids; dev compressed to 1.6 bits/row but full-scale saturation interleaves repeats — a token scheme (new-class bit + backref) should recover most of this |
| rows.pos | 13.2 | 2.8 | |
| rows.weight | 13.1 | 2.8 | headroom: geometric model |
| patterns + rows.pattern | 12.3 | 2.6 | exact E4 multimapper evidence |
| rows.shape + dict | 11.7 | 2.5 | |
| edges | 5.5 | 1.2 | cell-scoped, no window (D20) |
| flags, cells, misc | 2.8 | 0.6 | |

## Two bugs the acceptance discipline caught, worth keeping

- **Class-id delta coding** initially cost 23 of 49 bits/row because ids were assigned in
  cell-major order against genome-major rows — the packed-barcode mistake in a new costume.
  Renumbering by row order fixed dev (9.8 MB → 0.66 MB); full scale shows further token-coding
  headroom (above).
- **Class scope (D20):** per-(cell,chrom,strand) classes over-counted 2.9% because a multimapper's
  anchor chromosome can differ from its counted gene's (D16). The first full-scale fidelity read
  1.68% and the size table plus fidelity check caught it same-day; classes are now global per
  (cell, value), and the positional edge window was dropped rather than assumed safe.

## Still open (deliberately)

Per-bin chunking, the query indexes (junction postings, cell postings, site table), and the two
stream-coding improvements above — all v0.1 work, none a format break (sections are versioned and
skippable). D1 re-verification of every number quoted here is mandatory before any appears in a
paper table.


---

# v0.2 addendum — molecule-major layout (+ cell frequency renumbering)

**Date:** 2026-08-29 · same acceptance tests, re-passed: **regression byte-identical, fidelity
0.244% / 0 of 1,225 cells, identical replay output (10,393,179 UMIs).**

| | v0 (flat rows) | v0.2 (molecule-major) |
|---|---:|---:|
| size | 151.1 MB | **117.4 MB** (−22.3%) |
| vs CRAM 3.1 | 9.30× | **11.96×** |
| vs ingest BAM / FASTQ / Malva | 27.9× / 33.4× / 2.60× | **35.8× / 42.9× / 3.34×** |
| replay (archive) | 25.3 s | 47.4 s (load 1.5 s)* |

*Replay-time variance across runs is host-load noise (the BAM reference took 126–140 s across the
same runs); the stable claim is the load step: 1–2 s from the archive vs 80–91 s of BAM decode.

The diagnosis held exactly: the win had shrunk because flat rows re-paid cell/class/weight per row
(1.77 rows/molecule at full saturation), not because the data demanded it. Molecule-major storage
with a class token stream (fresh = 0, else backref) and frequency-renumbered cells recovered it:
cell 46.5 → 28.7 MB, class 46.0 → 8.8 MB, weight 13.1 → 8.1 MB. Dev-scale predicted almost none of
this (−5%) because dev molecules average 1.08 rows — the fourth dev≠full instance, favourable
direction for once.

**Largest streams now:** mol.rep 32.8 MB (per-representative offsets + shape ids — the 2-rep
payload's real price, above the ≤15 MB spec estimate), mol.cell 28.7 MB, patterns 9.2 MB.

**#3 measured (adoption deferred to v0.1 with chunking):** 844,565 patterns, 2.70M alt entries;
**78.7% of alternatives share the anchor's shape** (a same-shape flag removes most per-alt shape
ids, est. 2–4 MB); 1.16M distinct (chrom, offset) hops at 2.33× sharing (hop-vocabulary
factorization est. 2–3 MB). Real but modest; the same-shape flag is clearly worth taking, the hop
vocabulary is optional-with-story.

Remaining honest gap to the 16–19× projection: mol.rep — the projection priced second
representatives at ≤15 MB without per-rep shape ids. Candidate v0.1 coding: same-shape-as-first-rep
flag plus rep-offset models. Estimated landing with those + #3: ~105–110 MB (~13×).


---

# v1.0 addendum — chunking, indexes, and queries

**Date:** 2026-08-29 · regression **byte-identical** through every change; replay output unchanged.

**Container v1:** seekable directory; sections readable individually. **Chunks:** 4 Mb genomic
bins, self-contained decode state (per-chunk class base, delta resets). **Indexes written at
ingest:** chunk range index, junction catalogue + postings, cell postings. **Same-shape pattern
flag adopted** (with anchor-relative shapes merging patterns 845k → 773k).

| | unchunked v0.2 | v1.0 @ 4 Mb | v1.0 @ 16 Mb |
|---|---:|---:|---:|
| size | 117.4 MB | **120.8 MB** (11.63× vs CRAM) | 119.4 MB |

Chunking + all three indexes + pattern savings net out to **+2.9%** — random access costs almost
nothing.

**Queries (full D0 archive):**

| query | `aie` result | latency | samtools on 1.4 GB CRAM |
|---|---|---:|---|
| region (EEF1A1, 36 kb) | 16,288 UMIs across 2,413 cells, per-cell counts | 0.14 s | 0.074 s — but returns 432,969 raw reads, no UMI dedup, no cells |
| junction (top SJ) | 60,156 UMIs across 4,511 cells supporting the exact junction | 0.10 s | 0.038 s — 170,344 raw region reads, no junction filter, no dedup |

The honest comparison: raw latencies are the same order (both index-seek + decode). The difference
is what comes back — the archive answers the biological question (per-cell molecule counts for a
feature) directly; the CRAM returns raw records that still need tag parsing, junction filtering
and UMI collapse downstream. That, plus 11.6× less storage, is the claim — not "faster than
samtools at fetching bytes."

Semantics note (documented in `querycmd.rs`): region listings are anchor-based; a multimapper
molecule anchored elsewhere is reachable through junction/pattern paths, not region listings.
