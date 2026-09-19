# Memo — query-surface richness: assessment and redesign proposals

**Date:** 2026-08-30 · **Scope:** `aie query` / `aie federate` (querycmd.rs), read-only design
study. No source changes made; no format changes required for the "must-have" tier. Hands-on
verification against `runs/archive/d0.aie` (111 MB, v1.5).

---

## 1. The current query surface, precisely

Four `query` subcommands plus `federate`, all in
`src/crates/aie/src/querycmd.rs` over `LazyArchive` (archivecmd.rs):

| command | question answered | returns | measured latency (d0) |
|---|---|---|---|
| `query A region chr:s-e [--top N]` | which cells have molecules **anchored** in the window | summary line + top-N `barcode\tUMIs` | 0.16 s |
| `query A junction chr:d-a [--top N]` | per-cell class-deduped UMI counts supporting one **exact** junction | summary (total UMIs, cells, index support, chunks) + top-N `barcode\tUMIs` | 0.20 s |
| `query A apa chr:s-e [--site-gap 24] [--strand +/-] [--tsv] [--groups f.tsv]` | clustered 3′-end sites in a window, strand-aware; optional per-group UMI counts | summary or TSV `chrom start end strand umis cells [group…]` | 0.22 s |
| `query A discover --gtf g.gtf [--merge-gap 1000] [--min-umis 10] [--tsv] [--emit-gtf out]` | candidate unannotated loci (molecules unclaimed by the GTF, clustered) | summary or TSV `chrom s e strand umis cells`; optional AIENOVEL GTF | genome scan (~5–10 s class) |
| `federate A B… chr:d-a [--top 5]` | one junction across N archives | one line per archive: UMIs / cells / index support / top cells | 1.34 s for d0+d1 |

Adjacent but not "queries": `replay-rows` (any GTF → full MTX matrix, `--velocity`,
`--audit-multigene`), `em --emit/--star` (fractional multimapper layer), `debug` (byte
accounting).

**What the index actually holds** (format-spec v1.5) vs what queries touch:

- **Used:** chunk index (`index.chunks`: chrom, bin_start, n_mols, class_base, max_anchor,
  n_cells), junction catalogue + postings (`index.junctions`/`index.jpost`), coc blocks
  (cell-of-class), shapes, cells dict, per-molecule records (chains with span-extreme reps +
  read counts, mm anchors + pattern ids, strand).
- **Written but never read by any query:** `index.cellpost` (cell → chunk postings — a whole
  per-cell access path, unused; grep confirms zero readers), `edges` (UMI adjacency graph),
  `patterns` dictionary (only replay uses it), chunk-index `max_anchor`/`n_cells` (region query
  ignores `max_anchor` and instead decodes every chunk left of the window, skipping after
  decode), per-chunk `n_cells`.
- **In the spec narrative but not in the file:** the 3′-site table ("kept as the APA/discovery
  query index") was never persisted in v1.x — `apa` recomputes 3′ ends from molecule decode.

### Semantics fine print (matters for the comparison)

- **Region = anchor containment, not overlap.** A molecule anchored 1 bp left of the window
  that spans across it is invisible. samtools/tabix semantics are "any overlap".
- **Coordinates are 0-based**, `chrom:start-end`, comma-tolerant; junction locus is
  (exclusive end of upstream exon, start of downstream exon) — neither STAR SJ.out (1-based
  intron) nor BED; nowhere documented on the CLI.
- One region per invocation; no BED file of regions; no whole-chromosome form; no gene names.
- Junction lookup requires **exact** coordinates known in advance, and there is **no way to
  list junctions** (my first junction query failed until I found coordinates in a memo).
- Catalogue/postings walk is a linear prefix scan per query (fine at 170k–672k junctions, but
  it decodes all postings lists up to the target id).
- `junction_counts` decodes chunks even when only the support total (already in `index.jpost`)
  is wanted.

### Output & composability

- Human summary + top-N (default 20) to **stdout**, timings embedded in the summary line; no
  way to get *all* per-cell counts (`--top` only); no header rows except `apa --groups`.
- `--tsv` exists only on `apa`/`discover`; no JSON anywhere; no MTX from queries (only
  replay/em emit MTX); no `-o` (stdout only — fine) but summary lines pollute TSV pipes
  (region/junction have no TSV mode at all).
- SIGPIPE handled (pipes to `head` work). Exit code is nonzero on absent junction (good), but
  `federate` prints "junction absent" and exits 0 (inconsistent).
- `federate` runs archives serially; per-archive barcode namespaces (correct — no false
  cross-sample cell identity), but no set semantics across archives beyond a grand total line.

---

## 2. What users of comparable tools expect

**samtools/tabix (region access):** 1-based `chr`, `chr:beg-end`, *overlap* semantics,
multiple regions per call, `-R regions.bed`, streaming record-level output that pipes into
bedtools/awk, machine-readable by default. Gravex has none of these conventions and no
record-level output at all — the molecule is the format's whole point, yet you cannot export
one.

**Snaptron (the bulk-coverage analogue, per its docs: `regions` + `rfilter` + `sfilter`):**
- Region queries by coordinate **or gene name** (annotation-resolved), with `exact` /
  `contains` / `either=1|2` coordinate modifiers for junction matching.
- Interval queries return **all junctions in the region** as TSV with rich columns:
  samples_count, coverage_sum/avg/median, strand, **annotated / left_annotated /
  right_annotated** status.
- **Range filters** (`rfilter=samples_count>:5`, `coverage_sum>:10`) and **sample metadata
  filters** (`sfilter=tissue:cortex`).
- Client-side **high-level queries**: PSI (percent-spliced-in), JIR (junction inclusion
  ratio), TS (tissue specificity), SSC (shared sample count), and **set operations
  (union/intersection) across queries** in a batch file.

Gravex's pitch is "Snaptron at molecule/cell resolution" (querycmd.rs says so), but today it
matches only Snaptron's *point* junction query. No listing, no filters, no annotated status,
no PSI/JIR analogues, no set ops — even though per-cell resolution makes all of these
*stronger* here than in Snaptron (per-cell PSI is a thing Snaptron cannot do at all).

**10x ecosystem:** filtered barcode lists (`barcodes.tsv`) as the lingua franca for "these
cells"; MEX (MTX + barcodes + features) as the matrix interchange; `aggr`-style multi-sample
composition. Gravex accepts a barcode list only in `apa --groups` and in replay's column
ordering; queries cannot be scoped to the filtered cell set, and no query emits MEX.

---

## 3. Proposals — no format change required

Everything below is servable from streams already in the file. Latency classes:
**I** = index/dictionary only, <0.1 s, flat in archive size · **C** = targeted chunk decode,
0.1–1 s · **P** = postings walk over many chunks, 0.2–2 s/archive · **G** = full genome scan,
~5–10 s (d0) / ~25 s (d1). Effort: S < 1 day, M ≈ 1–3 days, L > 3 days.

### MUST-HAVE

**Q1. `query junctions <region>` — junction listing with filters (the Snaptron interval
query).** *The single biggest gap: you cannot ask "what junctions are here?".*
- User story: "Show me every junction in the CD44 locus with ≥20 supporting reads in ≥10
  cells, and whether it's annotated."
- CLI: `aie query A junctions chr11:35138870-35232402 --min-support 20 --min-cells 10
  [--gtf ref.gtf] [--contains|--either d|a] --tsv`
  → `chrom donor acceptor support n_cells n_umis annotated left_anno right_anno`
- Serves from: `index.junctions` is already **coordinate-sorted** (chrom, donor, acceptor) —
  binary-search/scan the window; `index.jpost` support totals answer `--min-support` with
  **zero chunk decodes** (class I). `n_cells/n_umis` and `--min-cells` need the postings walk
  (class P). `--gtf` marks donor/acceptor membership via the existing `anno` crate.
  `--contains d` (all junctions sharing a donor = alternative acceptors) is a catalogue scan.
- Latency: I without cell columns, P with. Effort: **M**.

**Q2. Uniform machine-readable output: `--tsv` / `--json` on every subcommand, `--top 0` =
all rows, headers, summary/timings to stderr.**
- User story: "Pipe per-cell junction counts into pandas without regex-scraping a summary line."
- CLI: `aie query A junction chr16:89562391-89562883 --tsv > counts.tsv` →
  `barcode\tumis` for *all* cells; `--json` for the summary object.
- Serves from: nothing new; pure presentation. Effort: **S**. (Also: exit 0-with-empty-output
  vs error made consistent; document 0-based coordinates and junction convention in `--help`.)

**Q3. Cell scoping and grouping everywhere: `--cells barcodes.txt` (filter) and
`--groups two-col.tsv` (stratify) on region/junction/junctions/apa, with
`--agg cell|group|bulk`.**
- User story: "Count this junction only in my filtered/annotated T cells, and give me
  T-vs-monocyte pseudobulk in one call" — the APA `--groups` machinery generalized (the code
  in `What::Apa` is already exactly this; it just isn't shared).
- CLI: `aie query A junction chr16:89562391-89562883 --groups celltypes.tsv --agg group --tsv`
  → `group n_cells n_umis`
- Serves from: cells dict (packed-barcode → id) + `cell_of` coc blocks. Effort: **S–M**.

**Q4. `query molecules <region>` — molecule-level export (the `samtools view` / BED
equivalent).**
- User story: "Give me the raw evidence at this locus as records I can bedtools-intersect,
  sort, and eyeball: span, strand, cell, reads, junction chain, multimapper status."
- CLI: `aie query A molecules chr1:155234000-155236000 [--bed|--tsv|--json] [--cells f]
  [--spliced|--unspliced] [--multi-only|--unique-only] [--min-reads N] [--min-span/--max-span]`
  → TSV: `chrom span_start span_end strand barcode umi_class n_reads n_chains n_mm
  junctions(d1-a1;d2-a2…) pattern_id` · `--bed`: BED12-ish with junction chain as blocks.
- Serves from: chunk decode + shapes (span/junctions reconstructed exactly as `discover` and
  `three_prime` already do); `patterns` dict for `--multi` detail. Spliced = any shape with
  ≥2 blocks; span filters from rep extremes. Overlap semantics: decode uses span, not anchor
  (fixes the region-semantics gap for free).
- Latency: C. Effort: **M**. *This one line item unlocks the entire unix-toolchain
  composability story and is the honest demo of "molecule resolution".*

**Q5. Junction set composition per cell: `query jset` (per-cell PSI / A-AND-NOT-B).**
- User story: "Which cells support the inclusion junctions of this cassette exon but not the
  skipping junction?" / "Per-cell PSI for exon X" — Snaptron's JIR at single-cell resolution;
  no other tool can answer this from a 111 MB file.
- CLI: `aie query A jset --incl chr1:100-200 --incl chr1:300-400 --excl chr1:100-400
  [--groups f] --tsv` → `barcode incl_umis excl_umis psi`; `--expr 'J1 & !J2'` as the general
  form later.
- Serves from: per-junction postings intersected at chunk level (usually the same 1–2
  chunks), then one decode pass evaluates the predicate per (cell, class); the class-dedup
  in `junction_counts` already gives correct molecule semantics.
- Latency: P (≈ one junction query, chunks shared). Effort: **M**.

**Q6. Region/locus conventions + scan fix.**
- Multiple loci per call and `--regions-file f.bed`; whole-chromosome `chr1`; `--gene NAME
  --gtf ref.gtf` resolves a footprint (annotation *projected at query time* — the project's
  own thesis); `--overlap|--anchor` semantics flag, overlap default.
- Fix chunk selection to use the already-stored `max_anchor`: today `region`/`apa` decode
  every chunk left of the window on the chromosome and discard (querycmd.rs:110–147, 196–216)
  — `c.max_anchor >= start && c.bin_start < end` prunes without decode.
- Serves from: chunk index + anno crate. Latency: C. Effort: **S**.

**Q7. Federation upgrades: k-of-N and sample-set semantics.**
- User story: "Which junctions in this region are present in ≥3 of my 6 samples? Which are
  private to the tumor archive?" — FED-3 already measured 44.2% junction recurrence d0→d1;
  turn that measurement into a query.
- CLI: `aie federate A B C… junctions chr17:7565097-7590856 --min-samples 3 --tsv` →
  `chrom donor acceptor n_samples support_per_sample…`;
  `aie federate … junction chr:d-a --tsv` → long-form `sample barcode umis`;
  `--only-in A` / `--absent-in B` for sample-specific junctions.
- Serves from: per-archive catalogues alone for presence/support (class **I per archive** —
  no chunk decode needed when cell columns are off); parallel archive opening (rayon, already
  a dependency). Latency: I–P. Effort: **M**.

### NICE-TO-HAVE

**Q8. `query cell <barcode>` — per-cell evidence profile (finally use `index.cellpost`).**
- User story: "This barcode looks like a doublet/novel state — show me all its molecules /
  its junction usage genome-wide."
- CLI: `aie query A cell AAACCTGAGCGATATA --tsv [--summary]` → molecule rows (Q4 schema) or
  per-chromosome summary.
- Serves from: `index.cellpost` (written today, read by nothing) → decode only that cell's
  chunks. Latency: P (spread cells touch many chunks; still ≪ scan). Effort: **S–M**.

**Q9. MEX/MTX matrix output from queries: `--mtx out/` on `junctions` (junction × cell) and
`apa` (site × cell).**
- User story: "Load junction-level or 3′-site-level counts for this gene panel into
  Seurat/scanpy next to the gene matrix" — the alternative-splicing-in-scRNA workflow
  (comparable to STARsolo SJ.out, but queryable per locus after the fact).
- CLI: `aie query A junctions --regions-file panel.bed --mtx out/ --barcodes barcodes.tsv`
  → `matrix.mtx` + `features.tsv` (junction ids) + `barcodes.tsv`; emit_matrix in
  archivecmd.rs is 90% of the writer already.
- Latency: P × regions. Effort: **M**.

**Q10. Paralog-pattern queries.**
- User story: "Is the signal at this locus shared with a paralog? Which loci co-mumap with
  it, with what read mass?" — evidence no count matrix or BAM slice shows compactly.
- CLI: `aie query A ambiguity chr1:155234000-155236000 --tsv` → for mm molecules in the
  window: `pattern_id n_mols n_umis alt_chrom alt_offset strand_flip` (the cross-locus
  ambiguity map); `aie query A pattern <id>` dumps one pattern's placements (class **I** —
  dictionary only).
- Serves from: `patterns` dict + window chunk decode. Reverse queries ("all molecules whose
  pattern touches locus X") need either a genome scan (G) or a pattern→chunk postings index
  (format addition, §5). Latency: C forward / G reverse. Effort: **M**.

**Q11. Spliced/unspliced/span filters as global flags** (subsumed into Q4's filter set, but
also on `region`/`junctions` aggregations): `--spliced-only`, `--min-span`, `--min-reads`,
`--strand` — shape block-count, rep extremes, chain weights all in the decoded record.
Effort: **S** once Q4 exists.

**Q12. `apa` genome-wide mode + built-in differential statistic.**
- `aie query A apa --all --gtf ref.gtf --groups f --tsv` → per-gene 3′-site table with the
  93_apa_diff.py shift statistic inline (median-shift + null option). Latency: G. Effort: **M**.

**Q13. Arrow/Parquet output (`--format arrow`).** Worthwhile only after Q2/Q4 settle the
schemas; adds a heavy dependency (arrow-rs) to an otherwise lean binary. Effort: **M**,
deferred — TSV/JSON/MTX cover the immediate audiences.

### Priority order

1. **Q1 junctions listing** — removes the "must already know the answer" wall; mostly index-only.
2. **Q2 uniform TSV/JSON + `--top 0`** — cheapest credibility fix; everything else inherits it.
3. **Q4 molecule export** — the format's thesis made visible/composable; enables Q11.
4. **Q3 cell scoping/grouping + aggregation** — makes every query single-cell-native.
5. **Q7 federate k-of-N** — the atlas story with set semantics, index-only fast path.
6. Q5 jset/PSI · 7. Q6 conventions/scan fix (S, can ride along early) · 8. Q9 MTX ·
9. Q8 cell profile · 10. Q10 patterns · 11. Q12 apa-wide · 12. Q13 Arrow.

---

## 4. Sketch of the composed surface (after must-haves)

```
aie query d0.aie junctions chr11:35138870-35232402 --min-support 20 --gtf ref.gtf --tsv
aie query d0.aie junction  chr16:89562391-89562883 --groups celltypes.tsv --agg group --tsv
aie query d0.aie molecules --gene CD44 --gtf ref.gtf --spliced-only --bed | bedtools ...
aie query d0.aie jset --incl chrX:a-b --excl chrX:a-c --groups celltypes.tsv --tsv   # per-cell PSI
aie federate runs/archive/d*.aie junctions chr17:7565097-7590856 --min-samples 3 --tsv
```

---

## 5. Would need a format addition (recorded, not proposed for v1.x)

| addition | enables | cost estimate |
|---|---|---|
| strand byte in `index.junctions` (spec layout lists it; writer never stored it) | strand column in Q1 without decoding supporting molecules | ~1 bit/junction |
| pattern → chunk postings | reverse paralog queries (Q10) at class P instead of G | small (651k patterns, sparse) |
| persisted 3′-site table (the v0 decision that lapsed) | genome-wide APA (Q12) at class I/C instead of G | 2.84 MB measured in Phase-0 |
| per-chunk per-cell UMI histograms (or per-(chunk,cell) counts beside `index.cellpost`) | pseudobulk region counts and cell profiles without chunk decode | few MB; measure first |
| sample-metadata section (free-form JSON) | Snaptron-style `sfilter` at federate level | ~KB |
| per-junction per-cell postings (junction → (cell, count)) | class-I per-cell junction counts | large — likely rejected on bytes; measure |

All of §3 works without any of these; they only move latency classes leftward.

---

## 6. Bottom line

The current surface proves four capabilities (region, exact junction, APA, discover) with
excellent open latency (lazy open ≈ 0.00 s, flat in archive size), but it is a demo surface,
not a query language: point lookups with human-formatted top-N output, no junction discovery,
no filters, no set operations, no record export, no cell scoping, and a per-cell index
(`index.cellpost`) that nothing reads. The must-have tier is presentation and composition
work over streams the format already pays for — the format itself is *richer than its query
surface*, which is the best possible position to be in before a systems paper: every proposal
above is an eval-section table, not a format renegotiation.
