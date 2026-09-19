# Memo — codec choice, stream split, and index enrichments (v1.1)

**Date:** 2026-08-29 · **Scale:** full D0 · **Regression:** archive replay ≡ BAM replay,
**byte-identical** after every change below. Archive: `runs/archive/d0.aie` (= the z19 build);
sweep artifacts `runs/archive/d0-z{12,19,22}.aie`.

## Codec sweep — keep zstd, keep level 19, reject lz4 with numbers

| Codec/level | Size | vs CRAM 3.1 (1,404.9 MB) | Ingest wall | Region query | Junction query |
|---|---:|---:|---:|---:|---:|
| zstd-12 | 128.76 MB | 10.91× | 107 s | 0.12 s | 0.09 s |
| **zstd-19 (kept)** | **121.41 MB** | **11.57×** | 143 s | 0.14 s | 0.09 s |
| zstd-22 | 121.32 MB | 11.58× | 198 s | 0.10 s | 0.09 s |

The decisive fact: **query latency is flat across compression levels.** A gene-scale query decodes
1–2 chunks of ~4 Mb; the zstd decode of those chunks is milliseconds and invisible next to the
~0.07–0.08 s open cost, which is *dictionary* decompression (cells/shapes/patterns/junctions +
indexes), not chunk decode. Therefore:

- **lz4 is rejected**: it would trade size (lz4 typically 1.5–2× worse than zstd-19 on these
  streams) for decode speed we cannot observe. The premise of the question — that chunk decode
  costs latency — is measurably false in this format.
- **zstd-22 is rejected**: 0.07% smaller for +38% ingest time.
- **zstd-12 is a documented knob**, not the default: 25% faster ingest for +6.1% size. Ingest is
  paid once per dataset; size is paid forever.
- The real latency lever, if 0.07 s ever matters, is **lazy dictionary loading** (decode postings
  and only the dictionaries a query touches) — a reader change, not a codec change. Backlogged.

## Stream split — a falsified estimate, recorded

Hypothesis: interleaving mol.rep position/shape (and the four mm columns) inside one stream was
costing zstd context, estimated ~10 MB. Implemented: the 7 per-chunk streams became 11
(rep.pos/rep.shape split out; mm split into pos/shape/pattern/weight). Measured: **120.78 →
121.41 MB (+0.5%)**, and that delta also *includes* the new index payload below. At chunk
granularity zstd was already exploiting the locality within each frame; the estimate was wrong.
The split is kept anyway — it costs nothing, isolates each column for future per-stream coding
experiments (mol.rep remains the largest stream and the headroom target), and the layout guard
below makes it safe. Filed with the site-ordering flip as instance #6 of "estimates lose to
measurements."

## Index enrichments (implemented, in the 121.41 MB)

1. **Per-junction support totals in `index.jpost`**: genome-wide junction support (e.g. 115,997
   supporting children for the test junction) now answered from the index alone, zero chunk
   decodes. Postings rows gained a leading varint; the junction query handler skips it.
2. **Per-chunk summaries in `index.chunks`**: `max_anchor` (true span coverage — chunks are
   binned by anchor start, so a range query previously had to over-fetch neighbors) and `n_cells`
   per chunk, enabling coverage-style summaries from the index alone.
3. **`chunk_streams` in `meta`**: the reader (`read_dicts`) refuses an archive whose per-chunk
   stream count it does not understand (absent field = 11, for the transition builds). This exists
   because of a real incident: an 11-stream binary silently mis-decoded a 7-stream archive during
   the sweep and produced plausible-looking garbage (0/206 in the refined APA analysis) — caught
   only because the number was implausibly zero. Layout changes are now loud, never silent.

## Decision

`.aie` v1.1 = v1.0 + split streams + enriched indexes + layout guard, zstd-19 per section/chunk.
121.41 MB, **11.57× vs CRAM 3.1**, 3.24× vs Malva, regression byte-identical, queries ~0.1 s.
The codec question is closed with measurements; remaining size headroom is in stream *coding*
(mol.rep models, weight models, hop factorization), not in the entropy coder.
