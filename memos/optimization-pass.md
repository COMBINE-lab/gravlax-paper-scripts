# Performance-optimization pass over the `aie` workspace (2026-08-30)

Scope: `replay-rows`, `em` (all modes + `--emit` layer), `federate`, `query region|junction`,
under the absolute contract that every checked output stays **byte-identical** (matrices,
`em.mtx`, EM stdout numbers, query/federate listings). All changes live in `$P/src` (not
committed). Binary of record before any edit: `runs/opt/aie-baseline` (copy of the pristine
`cargo build --release` binary at the pinned tree).

Machine: shared 256-core box; rayon capped at 24 threads (the binary's own default). All inputs
page-cache-warm. `/usr/bin/time -v`, wall + max RSS.

## Baseline (pristine binary)

3 reps per command except the two ~25–30-minute EM d1 runs (1 rep each; see note under the
final table). Values are medians (rep range in the raw `runs/opt/base-*.time*` files).

| command | baseline wall (median) | baseline peak RSS |
|---|---|---|
| replay-rows d0 (v49) | 8.39 s (7.45–8.45) | 7.43 GB |
| replay-rows d1 (v49) | 34.25 s (31.44–34.94) | 35.0 GB |
| em d0 (mask .2, 4 modes) | 284.7 s | 9.00 GB |
| em d1 (mask .2, 4 modes) | 1811.5 s (1 rep) | 42.7 GB |
| em d0 --mask 0 --emit | 232.4 s | 8.98 GB |
| em d1 --mask 0 --emit | 1409.3 s (1 rep) | 42.4 GB |
| federate 6 archives (junction) | 4.24 s | 1.17 GB |
| query region d1 | 0.82 s | 0.28 GB |
| query junction d1 | 1.19 s (1.19–1.27) | 0.70 GB |

Internal stage split, baseline replay-rows d1: load 7.9 s / GTF 2.1 s / replay 13.1 s /
emit ≈ 2.9 s / process teardown ≈ 5.4 s.

## Where the time actually went (perf, dwarf call graphs)

1. **EM blend mode recomputed a global sum inside the innermost closure.**
   `unique_counts_global.values().sum::<f64>()` sat inside the per-target per-candidate weight
   lambda: O(targets × cands × genes) full-map traversals per iteration. This alone made
   `em d1` a 30-minute command. (The map never changes after construction, so the sum is a
   constant.)
2. **Replay aggregation lived in nested hash maps.** After the first fixes, perf on replay d0
   showed >50% of samples in map machinery: 18% `malloc` (inner-map allocations), 18% in the
   aggregation closure (≈10% hashing), 6% `or_default`, 6% map iteration, 4%+4% rehash — vs
   only ~9% in the actual classification (`overlapping`/`align_vs_transcript`).
3. **Per-row allocation churn in the classifiers.** Every one of 108M rows allocated a
   `Placement` (2 Vecs), an overlap hit list, and a gene list; mm rows repeated this per
   pattern alternative.
4. **`align_vs_transcript` was quadratic in exon count** (per-segment linear rescan of the exon
   array from index 0) and `overlapping()` chased a pointer into each candidate `Transcript`
   just to read its span.
5. **Archive load:** serial section reads, then a serial 9 GB `extend` concat of decoded chunk
   vectors; `MolRec`/`MolChain` carried two heap `Vec`s per molecule (~200M small allocations
   at D1 scale — also the multi-second teardown).
6. **Queries:** `region` decoded every chunk of the chromosome left of the window before
   discarding the irrelevant ones; `cell_of` decoded 64Ki-entry coc blocks serially per lookup
   (62.7% of a junction query's samples in `rans::decode`); `federate` visited archives
   serially and deep-cloned the shape dictionary per archive.
7. **EM E-step collected a fresh `Vec<f64>` per target per iteration** (~35% of EM wall in
   malloc/free/consolidate after the blend fix), and the GTF parse itself was a serial 1.5–2 s
   (String-per-line, char-level splitting, three key `String` allocations per exon line).

## Changes (all output-byte-identical; each kept only after re-measuring)

* **rows.rs / EM**
  - Hoisted the blend-mode normalizer to a single pre-computed constant (identical f64: same
    map, same deterministic FxHashMap iteration order, never mutated).
  - E-steps parallelized over targets and written in place into `resp` (no cross-target
    accumulation happens in an E-step; per-target arithmetic order unchanged). M-step stays
    serial: its float accumulation order is part of the frozen output.
  - `pi_*.clone_from(&base)` instead of fresh clones per iteration.
  - Per-row candidate lists as `SmallVec<[u32; 4]>` (1–2 genes inline for almost all rows).
* **rows.rs / replay**
  - Fused flatten+classify+shard: rows are classified straight off the molecule records into
    64 cell-shards (per-worker piece lists merged per shard in parallel). The 108M-row flat
    table, its canonical sort, and the per-row `Option` vector are gone from this path — row
    order cannot reach the matrix because every shard fully sorts. (`flatten()` — sorted — is
    unchanged and still used by every EM/velocity path, where construction order feeds f64
    accumulation and IS part of the frozen output.)
  - Aggregation/collapse rewritten hash-free: shard sorted once; (cell,class) gene sums,
    MultiGeneUMI_CR best-gene, and the per-(cell,gene) tie-merge collapse all read contiguous
    runs; 1MM adjacency as CSR arrays (offsets+neighbors) instead of a HashMap of Vecs; rank
    lookup by binary search in the (small) group. Identical semantics, no hashing at all.
  - One-entry memo in the per-worker scratch: (chrom, pos, shape, pattern, strand) fully
    determines a row's gene set, and rows arrive position-clustered.
* **anno**
  - GTF parse: two-phase parallel — parallel `pread` of the whole file, parallel byte-level
    line/field/attr scanning into per-exon records, then a serial in-file-order intern pass
    (identical first-occurrence id assignment). No per-line String, no UTF-8 validation of
    3.3 GB, keys validated only on first intern.
  - `align_vs_transcript`: junction concordance and the block walk now binary-search the
    (disjoint, sorted) exon array — O(log E) per segment instead of O(E), same first-hit
    semantics (unit tests unchanged and passing).
  - `overlapping_into` (allocation-free) + flat `starts/ends` arrays in the chrom index (no
    per-candidate `Transcript` deref); `concordant_genes_into` skips classifying further
    transcripts of an already-admitted gene (same set, same order).
* **evidence-io**
  - `SectionReader::read_compressed_at`: `pread`-based `&self` section reads → chunk / coc
    reads parallelize end-to-end (read + decompress + decode fused per worker).
  - rANS decode: escape bit-gather byte-at-a-time instead of bit-at-a-time (bit-identical).
* **archivecmd**
  - `MolRec.chains`/`.mms`/`MolChain.reps` are `SmallVec`s sized for the dominant shapes
    (1 chain / ≤2 reps / ≤1 mm) — the ~200M per-molecule heap Vecs at D1 scale are gone.
  - `read_archive`: fused parallel pread+decompress+decode per chunk; parallel move-concat of
    decoded chunks into the final Vec (disjoint ranges; sources freed without element drops).
  - `emit_matrix`/em emit: digits pushed directly (byte-identical to `format!` for u32; floats
    keep the std formatter), lines formatted in parallel blocks.
  - `replay-rows`/`em` exit via `exit_without_teardown()` (flush stdio, `process::exit(0)`) —
    serially freeing the multi-GB molecule table cost whole seconds at D1 scale.
  - `LazyArchive`: shape dict behind `Arc` (no per-query deep clone), `prefetch_coc()` decodes
    all touched cell-of-class blocks in parallel before per-class lookups.
* **querycmd**
  - `region`: chunk selection via the index's `max_anchor` — the exact test the old code ran
    on decoded molecules (same selected set, same printed chunk count), then fused parallel
    read+decode; coc prefetch before counting.
  - `junction`/`federate`: fused parallel read+decode+junction-scan per postings chunk; coc
    prefetch; federate queries its archives in parallel and prints per-archive results in
    argument order (identical bytes).

## hashbrown / foldhash status

`hashbrown = "0.17"` (default foldhash) is a workspace dependency, used where iteration order
provably cannot reach an output: all `anno` lookup/index maps, `LazyArchive::coc_cache`, and
the emit-side barcode→column maps. Two deliberate carve-outs, both mandated by the byte-identity
contract or by measurement:

1. **Determinism.** foldhash seeds per process, so map iteration order differs run to run.
   Everywhere iteration order leaks into an observable — EM's `classes` map (its iteration
   order defines target order, and f64 accumulation is not associative, i.e. it reaches
   `em.mtx` and the printed EM numbers) and the ingest-side interning that fixes archive bytes
   — the maps stay rustc-hash `FxHashMap`, whose order is deterministic. Putting foldhash there
   would not just change outputs, it would make them nondeterministic. (Verified the other
   direction: `em.mtx` and EM stdout are byte-identical to the pre-change references.)
2. **Measurement.** The first cut moved the replay aggregation onto hashbrown/foldhash nested
   maps: replay d0 got SLOWER (replay stage 2.5 s → 4.5 s; perf: ~10% of run in
   `folded_multiply` alone — foldhash loses to fxhash on tiny integer keys). Per the mandate's
   "keep what wins", that stage is now hash-free (sorted runs + CSR), which beats both hashers.

## Benchmarks — baseline vs optimized (medians of 3 reps, wall / peak RSS)

| command | baseline | optimized | speedup | RSS baseline → optimized |
|---|---|---|---|---|
| replay-rows d0 | 8.39 s | 2.62 s | **3.2×** | 7.43 → 6.97 GB |
| replay-rows d1 | 34.25 s | 8.08 s | **4.2×** | 35.0 → 18.5 GB |
| em d0 | 284.7 s | 12.93 s | **22.0×** | 9.00 → 8.14 GB |
| em d1 | 1811.5 s¹ | 64.1 s | **28.3×** | 42.7 → 38.1 GB |
| em d0 --emit | 232.4 s | 13.27 s | **17.5×** | 8.98 → 8.16 GB |
| em d1 --emit | 1409.3 s¹ | 73.3 s² | **19.2×** | 42.4 → 38.1 GB |
| federate ×6 (junction) | 4.24 s | 0.50 s | **8.5×** | 1.17 → 1.19 GB |
| query region d1 | 0.82 s | 0.04 s | **20.5×** | 0.28 → 0.21 GB |
| query junction d1 | 1.19 s | 0.25 s | **4.8×** | 0.70 → 0.50 GB |

¹ Baseline EM d1 runs measured once each (30.2 and 23.5 minutes per rep); at a ~20–28×
separation, rep variance is immaterial. ² em d1 --emit reps spanned 62.3–76.2 s (shared-box
noise); the em d1 reps were 63.5/64.1/76.1 s.

Shared-box caveat: the machine hosts other users' work; baseline reps 2–3 of the fast commands
overlapped the tail of the single-threaded baseline em d1 --emit run (≤1 busy core of 256), which
is within the same noise band as the numbers show (e.g. replay-d0 7.45 s quiet vs 8.4 s then).

## Byte-identity verification (transcript)

`$AIE` = optimized `src/target/release/aie`; `$BASE` = pristine `runs/opt/aie-baseline`.
References are either the pre-existing `runs/replay/*` outputs (produced by the original binary
in earlier campaigns) or fresh `$BASE` outputs generated in this pass (`runs/opt/base-*`).

```
# replay-rows d0 — vs pristine-binary output
$AIE replay-rows runs/archive/d0.aie --gtf annotations/gencode.v49.annotation.gtf \
    --barcodes runs/oracle/full/v49/.../raw/barcodes.tsv --out-dir runs/opt/new-d0-replay
cmp runs/opt/base-d0-replay/matrix.mtx  runs/opt/new-d0-replay/matrix.mtx    -> identical
cmp .../features.tsv, .../barcodes.tsv                                       -> identical

# replay-rows d1 — vs BOTH the pristine-binary output and the original archived reference
cmp runs/opt/base-d1-replay/matrix.mtx  runs/opt/new-d1-replay/matrix.mtx    -> identical
cmp runs/replay/d1-aie/matrix.mtx       runs/opt/new-d1-replay/matrix.mtx    -> identical

# EM stdout (mask 0.2; per-cell + pooled + blend numbers, calibration curve)
diff runs/opt/base-em-d0.out runs/opt/new-em-d0.out                          -> identical
diff runs/opt/base-em-d1.out runs/opt/new-em-d1.out                          -> identical

# EM additive layer (mask 0 --emit) — vs the ORIGINAL runs/replay references
cmp runs/replay/d0-em-layer/em.mtx      runs/opt/new-d0-em/em.mtx            -> identical
cmp runs/replay/d1-em-layer/em.mtx      runs/opt/new-d1-em/em.mtx            -> identical
cmp features.tsv / barcodes.tsv for both                                     -> identical

# EM-0 STAR port (--star) — vs runs/replay/d0-em-star
cmp UniqueAndMult-EM.mtx / features.tsv / barcodes.tsv                       -> identical
  (this check caught a real bug in the teardown-skip: a BufWriter still in scope at
   process::exit lost its buffered tail — fixed with explicit flushes, re-verified)

# federate (6 archives) and query region/junction on d1 — $BASE vs $AIE stdout
diff modulo the elapsed-seconds fields (which vary run-to-run for any binary) -> identical
  counts checked explicitly: region 158036 mols / 148652 UMIs / 14747 cells / 2 chunks;
  junction 315157 UMIs / 29789 cells / support 711336 / 1 chunk;
  federate 1066128 UMIs / 93910 cells across 6 archives — all unchanged.

# beyond the mandate (same absolute contract):
query apa d0 (--tsv)         -> identical to $BASE
query discover d0 (--tsv)    -> identical to $BASE
replay-rows d0 --velocity    -> identical to runs/replay/d0-velocity (all 6 files)

cargo test --release: 71 passed, 0 failed (all crates).
```

## Remaining opportunities not taken

* **Parallel EM M-step / classes build.** The M-step's float accumulation order and the
  `classes` map iteration order are part of the frozen numbers; parallelizing them needs an
  ordered reduction (or a one-time blessed output change). Left serial (~half of the remaining
  `em d1` minute).
* **GTF parse**: still ~0.8 s (dominated by the 3.3 GB page-cache read + serial intern pass).
  An mmap or a persistent compiled-annotation cache would remove it; not taken (I/O-layer
  change, marginal for the mandated commands).
* **`emit` velocity path** (three per-line `writeln!` matrices) and the `--star` EM emit keep
  their serial formatting except for shared helpers — not in the mandated benchmark set.
* **`extract_rows` / ingest** untouched beyond the type changes: not a mandated workload, and
  its interning order fixes archive bytes.
* **Junction catalogue walk** in `query junction` is still a serial varint scan of the whole
  catalogue up to the hit (~0.15 s on d1); a binary-searchable layout needs a format change
  (frozen spec).
* **RSS of `em d1`** (~38 GB) could drop by sharding `per_row`→`classes` construction, but the
  map's iteration order is output-bearing (see above).
