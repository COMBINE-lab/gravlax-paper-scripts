# `.aie` archive format — specification draft

**Status: v1.5 IMPLEMENTED and accepted 2026-08-31** — v1.4 plus a per-block cell-of-class codec choice (delta-varints vs rANS over absolute ids via a 6th global table; one codec byte per block; `meta.codec = "rans2"`). Byte-identical on four datasets; sizes 111.1 / 544.0 / 439.7 / 595.2 MB (11.3–18.1 bits/read; 9.0–12.7× vs CRAM). Previously — **v1.4 accepted 2026-08-30** — v1.3 plus a static rANS stage on the five
memoryless streams (class, weight, rep.pos, mm.pos, mm.weight): global tables in `rans.tables`
(chunk decode stays self-contained), `meta.codec = "rans1"` guarded. D0 112.89 MB (12.44× vs
CRAM), D1 554.41 MB (11.76×); the D1 null result confirms the entropy layer is at its floor
(memos/structure-rans-v14.md — structural meta-color factorings also measured and closed there).
Previously — **v1.3 accepted 2026-08-29** — v1.2 plus (D29): cell-of-class stored in
65,536-class blocks, each its own zstd section (`coc.<b>`; first value absolute varint, rest
svarint deltas; `meta.coc_block` MUST equal the reader's block size — guarded like
`chunk_streams`). Queries open lazily (meta + chroms + chunk index; dictionaries on first use;
coc blocks decompressed individually as classes are touched) — open 0.00–0.01 s, flat in archive
size; full loads decode chunks and coc blocks in parallel. Writer compresses all sections in one
rayon pass (24-thread default cap). Sizes: D0 114.07 MB (12.32× vs CRAM), D1 554.51 MB (11.76×) —
+0.1–0.3% vs v1.2 for the independent coc frames. Replay 6.4 s (D0) / 25.3 s (D1).
Previously — **v1.2**: the D26/D27 factorizations:
(a) cell is a pure function of the UMI class (global (cell,value) classes, D20; verified 0
violations), so it lives in a global `cellofclass` section (svarint deltas, one per class in
fresh-introduction order) instead of per molecule in the chunks — chunk decode resolves cells
through this table (backrefs routinely cross chunks: same cell+UMI at another locus);
(b) each chain stores its span-minimum rep first and chains are position-sorted per molecule, so
the first stored rep IS the anchor for every molecule with chains — its rep.pos offset is elided
with no flag (implied by n_chains > 0; writer enforces the invariant). Chunks carry **10 streams**
(anchor, class, layout, weight, rep.pos, rep.shape, mm.pos, mm.shape, mm.pattern, mm.weight);
`chunk_streams` in `meta` guards the layout — readers MUST refuse a count they do not understand
(absent = pre-guard 11-stream file, refused). **113.73 MB, 12.35× vs CRAM 3.1**, regression
byte-identical; open cost 0.17 s (eager cellofclass decode — lazy dictionary loading is the
recorded lever). A per-molecule flag variant of (b) was measured and rejected: the flag bit
destroyed the layout stream's run structure and cost what the zeros cost (D27). Earlier v1.1: 11
split streams, per-junction support totals in `index.jpost`, `max_anchor`/`n_cells` per chunk. Codec CLOSED by full-scale sweep
(`memos/codec-indexing.md`): zstd-19 per section/chunk (z12: −25% ingest/+6.1% size documented
knob; z22 and lz4 rejected — query latency is flat across levels because open cost is dictionary
decompression, not chunk decode). **121.41 MB, 11.57× vs CRAM 3.1**, queries ~0.1 s, regression
byte-identical. Previously v1.0: seekable container, 4 Mb genomic chunks, junction/cell/range
indexes, same-shape patterns; 120.8 MB (11.63× vs CRAM). Previously: v0.2 (`memos/archive-v0.md`): molecule-major layout, class tokens, frequency-renumbered cells — **117.4 MB, 11.96× vs CRAM 3.1**, byte-identical regression, 0.244% fidelity, archive load 1.5 s. Decisions frozen 2026-08-28 on Phase-0 measurements (`results/sigstats-full-v2.json`,
`results/chainrep-full-v49.json`). Structure was pre-registered first; every decision below cites
its measurement. Remaining open: block size and graph sharding (freeze on writer benchmarks).

## ⚖ RESOLVED decision points (Phase-0, full D0: 21,236,782 molecules, 52,283,887 signatures)

| Decision | Verdict | Evidence |
|---|---|---|
| Signature payload | **FROZEN: two representatives per junction chain — the most-contained and most-extended reads (span extremes) — plus the chain read count** | The measured frontier (full D0, vs 0.217% read-level ceiling): 1-rep 0.619% (305 cells >1%); **2-rep 0.248% (0 cells >1%)**; exact ≈ read-level at +16 bits/mol. The extremes bracket every intermediate read's containment behaviour, so 2-rep recovers all but 0.031pp at the cost of one extra placement only where a chain's spans actually differ. Dev-scale said 1-rep cost 0.036pp and did not extrapolate (saturation, third occurrence of the standing warning) — which is why every point was measured at full scale before this froze |
| Multimapper evidence | **Paralog pattern dictionary** (offset vectors from anchor), interned | 5.11M mm signatures share 651,082 patterns (**7.85×**; top-100 cover 34.4%); dict 5.85 MB + id stream 5.33 MB = **4.21 bits/mol** for exact E4 union-rule evidence. Future shave: frequency-modeled ids |
| Position coding | **Flat svarint deltas.** Site two-level coding REJECTED for the codec | With dense cell interning: flat pos+cell = **14.61 bits/mol** vs site-ordered 15.03. The 46% cell-stream cut (11.11 → 6.01) loses to the position-monotonicity cost (3.49 → 9.02). First run of this experiment used packed barcodes and pointed the other way — the interning fix flipped the verdict, which is exactly why nothing froze early |
| Cell stream | (pos, cell) sort, dense interned ids, svarint deltas | same measurement |
| 3′ site table | Kept as the **APA/discovery query index only** (1.78M sites, 2.84 MB anchors), decoupled from the position codec | site table is cheap and query-motivated regardless of codec verdict |

**v0 size projection (full D0), 2-rep payload:** E1 base 54.93 MB + adjacency graph 3.46 MB +
paralog patterns 11.18 MB + chain counts ≈ 1 MB + second representatives (only where span extremes
differ; ≤ ~15 MB upper bound, exact figure from the writer) ≈ **~72–86 MB → ~16–19× vs CRAM 3.1**
at **0.248% fidelity with full E4 multimapper evidence**. The three measured points (1-rep 70.6 MB /
0.619%, 2-rep / 0.248%, exact / 0.217%) are the paper's archive-richness-vs-replay-accuracy figure —
the source plan's §21 curve with real prices on both axes.

## Design principles (fixed)

1. **Fidelity by construction, not approximation.** The molecule payload is the multiset of
   distinct *read signatures* (a signature = the complete placement set of one read, with a read
   count). Replay over signatures reproduces read-level STARsolo semantics exactly — including the
   multimapper unique-gene union rule (D16). No fidelity is negotiated away by the format; only
   bytes are.
2. **Compression comes from structure the data actually has**, each shared with a query index:
   - geometry → **shape dictionary** (measured: 385k shapes / 19.9M molecules, 2.4 bits/mol);
   - junctions → **junction catalogue** ids + inverted postings (the junction-query index);
   - UMI values → **adjacency graph** (measured: 1.20 bits/mol vs 27; and replay *improves*);
   - multimapping → ⚖ **paralog pattern dictionary** (offset vectors between a read's alternative
     placements, conserved across a repeat family — freeze on sig-stats pattern-sharing numbers);
   - positions → ⚖ **3′-site two-level coding** (site anchors + offsets; site table doubles as the
     APA query index — freeze on sig-stats site-vs-flat numbers);
   - cells → ⚖ ordering experiment ((pos,cell) vs (site,cell) runs — same measurement).
3. **Access is genome-major.** Replay is a sequential genome-ordered scan; junction/region/
   discovery queries are range reads; per-cell access goes through postings, never a scan.
4. **Nothing is irreversible.** Versioned sections with offsets in the header; optional sections
   (E2 edits, E3 residual sequence) append without breaking readers; unknown sections are skipped.

## Layout

```
header         magic "AIE0", format version, section directory (offset, length, codec, version)
dictionaries   cells (packed u32 whitelist barcodes, insertion-ordered)
               shapes (block-length/junction-offset vectors)
               junction catalogue (chrom, donor, acceptor, strand)
               ⚖ paralog patterns (vec of (chrom, offset, flip) relative to anchor)
blocks         per genomic bin (~2 Mb target): columnar streams, zstd-19 per stream
                 anchor positions (⚖ flat delta vs site+offset)
                 shape ids · flags · cell ids (⚖ ordering) · n_reads
                 signature payloads: per molecule, (n_sigs; per sig: count, kind,
                   unique → anchor-relative placement; mm → pattern id)
graph          UMI adjacency per (cell, chrom, strand): same-value links, 1MM edges (3 Mb window)
indexes        genomic range → block · junction id → block postings · cell id → block postings
footer         per-stream byte accounting (the Gate-0 measurement, embedded), checksums
```

## ⚖ Decision points and their gating measurements

| Decision | Options | Freezes on |
|---|---|---|
| Signature payload cap | store all distinct signatures vs cap at k | sig-stats signature-multiplicity histogram; cap only if the ≥4 tail is both large and cheap to truncate at measured fidelity cost |
| Paralog patterns | interned pattern dictionary vs per-molecule alt lists | pattern sharing ratio and top-100 coverage; dictionary wins if sharing ≫ 1 |
| Position coding | flat svarint deltas vs site anchors + offsets | measured bits/mol, and whether the site table is wanted for APA regardless |
| Cell stream | (pos, cell) sort vs (site, cell) sort | combined pos+cell bits/mol — the trade must be judged jointly, not per stream |
| Block size | **FROZEN: 4 Mb** (parameter retained) | measured: 4 Mb = 120.78 MB, 16 Mb = 119.43 MB (+1.1% for 4× finer granularity); gene-scale queries decode 1–2 chunks either way |
| Graph placement | **FROZEN: one global edges section** (5.5 MB, loaded once) | classes are cell-global (D20), so shards would split edges across chunks for no access benefit |

## Acceptance tests (fixed)

- Round-trip: write → read → identical molecule set, bit for bit.
- **Replay regression: replay-from-archive ≡ replay-from-BAM, byte-identical matrices** on dev and
  full D0, for v49, v32, and w1.
- Honest Gate F: replay wall time measured from the `.aie` file, reported next to file size and the
  per-stream footer accounting.
- Re-verify on D1 before any number is quoted (project rule: single-dataset numbers don't
  extrapolate — Gate G-D's dev/full reversal is the standing warning).
