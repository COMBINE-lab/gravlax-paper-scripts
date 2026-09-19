# Memo — multithreading, lazy loading, and the D1 size breakdown (v1.3)

**Date:** 2026-08-29 · **Regression:** byte-identical on BOTH datasets after every change.
Driver: `logs/v13-pipeline.{sh,log}`. Archives rebuilt: `runs/archive/{d0,d1}.aie`.

## 1. Were we multithreaded? Partially — and the serial parts were the wall-clock.

`replay_rows` already ran its two heavy stages (per-row gene assignment, per-(cell,gene)
collapse) under rayon. Everything between and around them was serial, and at D1 scale that serial
share dominated: the hashmap aggregation between the parallel stages, the flatten sort of ~190M
rows, the zstd-19 compression of every chunk at write, the chunk decompress+decode at load, and
the whole BAM extraction. Changes (all under a rayon pool now capped at **24 threads** by default
— shared-host discipline; `RAYON_NUM_THREADS` overrides):

- **Write:** sections are built serially, then compressed in one rayon pass and written in order
  (`SectionWriter::section_precompressed`).
- **Load:** chunk payloads are read serially (I/O) and decompressed+decoded in parallel; coc
  blocks decode in parallel.
- **Replay:** the aggregation → best-gene filter → tie-merge pipeline is sharded by cell (64
  shards, everything downstream of assignment is cell-keyed) and each shard runs in parallel;
  `flatten` uses a parallel sort.

Measured (D0 / D1):

| | v1.2 serial | v1.3 parallel |
|---|---:|---:|
| replay total | 49.9 s / 181.9 s | **6.4 s / 25.3 s** |
| …replay phase alone | 43.6 / 161.0 | 2.4 / 13.5 |
| ingest wall | ~150 s / 795.6 s | 103.7 s / 609.7 s |
| archive load | 1.5 s / 7.8 s | 0.9 s / 6.2 s |

Replay vs fresh STARsolo: **~150× (D0), 48× (D1)** — the D28 "6.5×" caveat is retired.
The ingest bottleneck is now BAM extraction (~540 s of D1's 610 s, serial noodles decode);
a multithreaded bgzf reader is the identified next lever, not done here.

## 2. Lazy loading — implemented as blocked cell-of-class (format v1.3)

Eager dictionary decode cost 0.69 s at D1 open, 76% of it the cellofclass table. v1.3 stores
cell-of-class in 65,536-class blocks, **each its own zstd section** (`coc.<b>`, first value
absolute, rest deltas; `meta.coc_block` guarded like `chunk_streams`). Queries open with meta +
chrom names + chunk index only; cells/shapes dictionaries load on first use; coc blocks
decompress individually as classes are touched — and because classes are numbered in genome
order, a region query touches a handful of blocks. Full loads (replay) decode all blocks in
parallel, so nothing is lost there.

Measured: **open 0.17 → 0.00 s (D0), 0.69 → 0.01 s (D1)**; region/junction totals 0.09–0.50 s.
Price: +0.34 MB (D0) / +0.60 MB (D1), ≈0.1–0.3% — independent frames compress slightly worse
than one. Open cost is now flat in archive size; the scaling concern from D28 is closed.

## 3. D1 size breakdown (554.5 MB) and where squeeze remains

Sections: chunks 337.5 MB, cellofclass 137.8, patterns 29.9, edges 23.8, shapes 11.0, cells 5.1,
indexes ≈9.3. Value-entropy bounds inside chunks (order-0, memoryless): rep.shape 108.7 MB,
rep.pos 88.2, class 53.8, weight 40.6, layout 33.1, anchor 33.3, mm.* ≈75 raw / ~30 bound.

- **cellofclass (24.8% of file) is the biggest single object** and the main modeling target:
  91.3M entries, delta-coded over a 1.27M-barcode space. Headroom would come from a cell-locality
  model (same cell's molecules cluster at loci → delta-0 runs) or rank/frequency coding with a
  real entropy coder — needs its own measurement pass before believing any number.
- **Dictionary factorizations improve with scale**: junction-chain sharing 6.89× (D0) → **8.20×**
  (D1); hop vocabulary 8.46× → **12.23×** (330k chains for 2.7M shapes — sublinear growth). Still
  only ~2–4 MB at D1, but the trend says revisit at D2+; the chains are also a queryable
  isoform-skeleton object, which may matter more than the bytes.
- **Pattern-alt vocabulary stays negative at scale** (2.19× sharing, 3.85M unique offsets) —
  closed twice now.
- Elision health at D1: only 0.7% residual zeros in rep.pos (duplicate-position chains); class
  tokens 84% fresh.

## Known wart

`aie query … | head` panicked with a broken-pipe message after `head` closed stdout — cosmetic;
FIXED same day (see addendum).

## Addendum (same day) — SIGPIPE, MT BGZF, fairness

- **Broken-pipe wart fixed**: `main()` restores SIGPIPE to default (libc), so `aie query | head`
  now exits silently like any Unix tool. Verified.
- **noodles**: already latest (bam 0.95.0, bgzf 0.51.0). Both extraction passes now read through
  `bgzf::io::MultithreadedReader` (≤8 workers — beyond that the serial record-parsing consumer is
  the limiter). Measured: extraction 93 → 65 s (D0), 520 → 389 s (D1); from-bam replay total
  722 → 448 s (D1). Regression byte-identical on both datasets and both paths — the MT reader is
  transparent to record order, and the matrices prove it.
- **Fairness of the STARsolo comparison**: D1 is a same-budget comparison — STARsolo at
  `--runThreadN 24` (1,225 s) vs replay under a 24-thread rayon cap (25.3–30.5 s wall), same idle
  host. D0 is conservative in STAR's favor (oracle ran at 32 threads vs our 24). Two caveats to
  carry into the paper: STARsolo's wall includes writing the 20 GB sorted BAM we requested (a
  matrix-only run would be faster), and STARsolo is doing more work by construction (alignment
  from FASTQ) — which is the claim, not a confound. Publish both measured walls, not one ratio.
- **Parallel emit/flatten**: matrix formatting and row-building were parallelized; both are
  within measurement noise at current scale (the replay phase was already collapse-dominated).
  Kept — they are the pieces that would surface next at D2 scale.
- **Remaining serial hotspot, recorded**: the extraction consumer (record field scan, interning,
  per-locus clustering) — ~389 s of D1's ingest. The fix is a batch pipeline with shard-local
  interning and a merge step; real surgery on extract_rows, deferred until ingest time matters
  (it is paid once per dataset). Also available if wanted: parallel multi-chunk decode in
  discover/apa full scans.
