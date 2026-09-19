# Decision log

Chronological record of choices that would be expensive or misleading to reconstruct later.
Detailed records live in `docs/decisions/`.

## 2026-08-28

**D1 — Primary oracle is STARsolo, not alevin-fry/simpleaf.** *(user-confirmed)*
The archive is genome-coordinate based; STARsolo is genome-based and therefore structurally
matched. alevin-fry maps to a transcriptome/splici index and so commits to an annotation at
mapping time, which cannot be replayed from genomic coordinates even in principle. simpleaf stays
as a second oracle family after Gate -1, where evidence sufficiency becomes genuinely hard
(selective-alignment scores depend on transcript sequence, so E2/E3 should start to matter).

**D2 — Gate order deviates from the source plan: storage is front-loaded.** *(user-confirmed)*
The source doc puts the compression audit at M4. Fidelity is very likely to pass (see D4), while
storage is what killed the predecessor project — and killed it on day ~12. G-E now runs in week 1.

**D3 — Rust from the start.** *(user-confirmed)*
Costs some up-front I/O plumbing during the phase where the design changes fastest, but the
archive writer is reusable without a rewrite. Python (the predecessor's `env/py312`) is for
figures and report tables only.

**D4 — E1 is a sufficient statistic for STARsolo gene assignment, by construction.**
Read from STAR's source on this host. `Transcriptome_classifyAlign.cpp::alignToTranscript` and
`Transcriptome_geneCountsAddAlign.cpp` consult only aligned blocks, splice junctions, strand and
`NH`. No sequence, no quality, no alignment score. Consequences: Gate -1A is near-tautological;
the E0→E4 ablation has a largely predictable answer; the real risk is elsewhere (G-C, G-D, G-E).

**D5 — UMI collapsing is annotation-dependent, and this is the project's crux.**
`SoloFeature_collapseUMIall.cpp` collapses in a `for (iG=0; iG<nGenes; iG++)` loop, i.e. per
(cell, gene). The storage win comes from collapsing reads into molecules, but the fidelity target
needs a collapse that consults the annotation. Gate G-D exists to bound the resulting error.

**D6 — Blocks split at CIGAR `N` only.**
STAR walks across indel separators when matching exons (`alignToTranscriptMinOverlap` expands the
block while `canonSJ` is negative and > -3), and advances the exon cursor only at `canonSJ >= 0`.
So a `D` extends a block rather than splitting it, and an `I` is skipped. Encoded and unit-tested
in `ingest::cigar`; getting this wrong would make replay demand exon boundaries that STAR does not.

**D7 — Contributions (1) and (3) demoted on prior-art grounds.**
`docs/prior-art.md`. The barcoded BAM already is a molecule-indexed annotation-independent evidence
store that people re-quantify today; recount3/Monorail owns the annotation-free-archive-plus-replay
framing for bulk; UMI-tools/BKC/UMI-nea own annotation-free UMI collapsing. The surviving spine is
the **replay-fidelity taxonomy** plus the **compact format with indexed compatibility queries**.

**D8 — Per-annotation STAR indices pre-built rather than inserting junctions per run.**
`scripts/05_build_anno_indices.sh`. Junction insertion is per-run overhead on many runs; paying it
once per annotation is strictly cheaper. The annotation-free base index (`ref/star-base`) remains
separate and is the only thing the archive itself ingests from.

**D9 — G2 grouping must read `UR`, not `UB`.**
STARsolo's UMI *correction* is also performed per-gene, so grouping annotation-independent
molecules on the corrected `UB` tag would smuggle the annotation back in and overstate our
fidelity. The oracle therefore emits `CR UR` alongside `CB UB`, and Gate G-D groups on `UR`.
Caught before the full-scale runs; the dev-scale oracle was re-run to add the tags.

**D10 — Dev-scale validation before every full-scale run.**
The 4M-read `dev` subsample caught two parameter faults within minutes that would each have cost
hours at full scale: STAR refuses to emit `CB`/`UB` into an unsorted BAM, and the `CR`/`UR` tags
of D9 were initially missing. Kept as standing practice.

**D11 — An annotation-free ingest cannot get STARsolo's corrected barcodes.**
STARsolo refuses to emit `CB`/`UB` unless a gene feature (`Gene`/`GeneFull`/...) is requested, and
those require a GTF. So `scripts/40_ingest_align.sh` emits only raw `CR`/`UR` and the archive must
correct barcodes itself against the whitelist. This is the honest arrangement in any case: barcode
correction is whitelist-based and genuinely annotation-independent, whereas STARsolo's UMI
correction is per-gene and must not leak into our grouping.

**D12 — Gate G-C BAMs must be name-sorted before joining; STAR's unsorted order is not stable.**
Two STAR runs over identical input FASTQs do **not** emit records in the same order — the streamed
join detected a read-name mismatch at record 168,642 between the two-pass and sjdb_v49 runs. The
comparator's order guard caught it rather than silently reporting nonsense, which is why the guard
exists. Fix: `samtools sort -n` every G-C BAM before comparison.

**D13 — Gate G-D's operating point is 50 kb, not the 10 kb chosen at dev scale.**
The 4M-read subsample sits at 17.2% sequencing saturation; full D0 is at 70.0%. That reversed two
dev-scale conclusions: UMI error correction went from negligible (0.13 pp) to essential (7.14% ->
1.21%), and the error curve does not flatten until ~50 kb rather than ~10 kb. The headline error
roughly doubled, 0.49% -> 0.94%. **Dev-scale gate results are for catching bugs, never for
predicting full-scale numbers** — recorded in `memos/gate-D-decision.md` rather than quietly fixed.

**D14 — bioconda's STAR installs as `STAR-avx2`; never health-check it with `pgrep -x STAR`.**
The conda package dispatches through a wrapper, so the running process is named `STAR-avx2`.
Exact-match checks (`pgrep -x STAR`, `ps -C STAR`) return nothing for a perfectly healthy run.
This produced a false "process died" reading on the full-scale ingest, and acting on it launched a
*second* STAR against the same `--outFileNamePrefix` while the first was still running — two writers
on one output directory. Both were killed and the run restarted clean; no results were derived from
the overlapping period. Substring matching (`ps aux | grep '[S]TAR'`) was correct all along, which
is why earlier polls looked fine and the later exact-match ones did not.

**Lesson recorded rather than silently fixed:** a liveness check that is wrong in the *negative*
direction is worse than none, because the recovery action it triggers (relaunch) actively corrupts
a healthy run. Monitors must be validated against a known-running process before being trusted.

**D15 — the host is shared and genuinely memory-contended.**
`dmesg` shows a global OOM in which another user's python process reached 1.4 TB RSS. Our jobs cap
threads at 24-32 and `--limitBAMsortRAM` at 32 GB (reduced from 80 GB) to keep our own footprint
modest and our exposure low. Long runs should be launched with `setsid nohup` and an EXITCODE file
so their fate is unambiguous after the fact.

**D16 — STARsolo's `Gene` counts unique-GENE reads, not unique-ALIGNMENT reads, and the `GX` tag
hides it.** A read mapping to several loci (NH>1) is still counted when the concordant genes of
ALL its alignments reduce to one gene; the primary record's GX often stays "-" while secondaries
carry the evidence. Found by entry-level matrix-vs-tags diffing: the dev matrix exceeds every
tag-derived count by exactly 37,590 UMIs (2.2%), concentrated on multimapping genes; verified on
EEF1A1, where one cell shows a matrix count of 102 against primary-GX tags of 4, with 118
secondary alignment records supplying the missing evidence. Consequences: (a) replay must union
concordant genes across a read's alternative placements (`--read-level` replay now does); (b) the
archive MUST retain alternative placements for multimapped molecules — the plan's §14.4/E4 field
is not optional, it carries 2.2% of matrix mass; (c) the tag-derived "G1" reference used by
gate-d and the graph experiment undercounts the true matrix by the same 2.2%, so their absolute
errors were measured against a slightly wrong target (relative policy conclusions unaffected;
re-basing needed before any paper number is quoted).

**D17 — the Gate -1A error decomposition, dev scale, after D16.** With multimapper union counting:
assignment+collapse port vs oracle on identical alignments+barcodes = 0.032% of UMI mass (zero
cells off by >1%); adding our own barcode correction costs +0.68pp; adding the annotation-free
alignment costs +0.15pp; end-to-end annotation-free replay = 0.869%. Barcode correction is now the
LARGEST fidelity component — our exact-or-unique-1MM rule drops reads that STARsolo's
1MM_multi_Nbase_pseudocounts resolves probabilistically (whitelist frequency + base quality), which
is implementable annotation-free and is the next fidelity lever, ahead of the junction catalogue.

**D18 — Phase-0 format measurements resolved four spec decision points (2026-08-28).**
(1) Exact read signatures are unaffordable: fragmentation gives 2.46 distinct placements per
molecule, so "fidelity by construction" at molecule prices was wrong; adopted instead: one
most-contained representative per junction chain (78.8% of molecules have exactly one chain;
1.004 chains/mol) + per-chain read counts, fidelity cost measured at 0.036pp (dev). (2) Multimapper
evidence via an interned paralog-pattern dictionary: 7.85x sharing, 4.21 bits/mol for exact E4
union-rule evidence. (3) Site-based position coding REJECTED by measurement: flat 14.61 vs
site-ordered 15.03 bits/mol (pos+cell, dense interning) — the first run of this experiment used
packed barcodes and pointed the OPPOSITE way; the interning fix flipped it. The 3' site table
survives as the APA query index only. (4) v0 projection: ~70.6 MB ≈ 19.9x vs CRAM 3.1 with full
multimapper evidence (vs 24.06x without). Spec: docs/format-spec.md; data:
results/sigstats-full-v2.json.

**D19 — signature payload frozen at TWO representatives per junction chain (span extremes).**
Measured frontier on full D0 against the 0.217% read-level ceiling: 1-rep 0.619% (305/1225 cells
>1%), 2-rep 0.248% (0/1225), exact +16 bits/mol. The two extremes bracket every intermediate
read's containment behaviour. Dev-scale priced 1-rep at 0.036pp and was wrong by 11x at full scale
(saturation; third occurrence) — all frontier points were therefore measured at full scale before
freezing. `results/{chainrep,tworep}-full-v49.json`.

**D20 — UMI classes must be scoped (cell, value) globally, not (cell, chrom, strand).**
First full-scale archive replay over-counted by 2.9% (fidelity 1.68% vs the emulation's 0.248%).
Cause: a multimapper's anchor chromosome can differ from the chromosome of the gene it counts for
(the D16 union rule), so per-(cell,chrom,strand) classes split the same (cell, UMI) across two
classes inside one gene — counted twice where value semantics count once. Classes are now global
per (cell, value), and 1MM edges are cell-scoped with no positional gate, exactly mirroring the
value-based collapse (an edge that never co-occurs within one gene is never consulted). The
positional edge window survives only as a future storage optimization to be re-validated against
fidelity, not assumed.

**D21 — `.aie` v0 accepted (2026-08-29).** Full D0: 151.1 MB (9.30x vs CRAM 3.1, 2.60x vs Malva),
fidelity 0.244% / 0 cells >1% (= the emulation, so the value-free class collapse costs nothing),
replay 25.3 s vs ~16 min fresh (~38x), archive load 0.8 s vs 80.4 s BAM decode. Regression:
archive-sourced and BAM-sourced matrices byte-identical. Honest placement: the frozen G-E >=10x bar
was defined for the E1 molecule archive (11.76x/24.1x); the full-capability .aie lands at 9.30x —
different objects, present the ladder with prices, not one number. Known coding headroom:
rows.class token scheme, rows.weight model. `memos/archive-v0.md`.

**D22 — molecule-major layout accepted (v0.2, 2026-08-29).** 151.1 -> 117.4 MB (11.96x vs CRAM,
back above the 10x bar with full capability), regression byte-identical, fidelity/output unchanged.
Flat rows had re-paid cell/class/weight per row (1.77 rows/molecule); molecule-major + class tokens
+ frequency-renumbered cells recovered it (class 46.0->8.8 MB, cell 46.5->28.7 MB). Dev-scale
showed only -5% (1.08 rows/molecule at 17% saturation) — fourth dev!=full instance. #3 measured:
78.7% same-shape alts (adopt in v0.1), 2.33x hop sharing (optional). Largest stream now mol.rep
(32.8 MB) — the true price of the 2-rep payload; coding headroom identified. `memos/archive-v0.md`.

**D23 — chunking & indexing campaign complete (v1.0, 2026-08-29).** Seekable v1 container; 4 Mb
chunks FROZEN (120.78 vs 119.43 MB at 16 Mb — +1.1% buys 4x granularity); global edges section
FROZEN (cell-global classes make shards pointless); junction catalogue+postings, cell postings,
chunk range index written at ingest; same-shape flag adopted (patterns 845k->773k). Net cost of all
of it: +2.9% (117.4 -> 120.8 MB, 11.63x vs CRAM). Queries live: region and junction with per-cell
UMI counts in ~0.1 s from a 120 MB file; the claim vs samtools/CRAM is semantics returned (counts
vs raw reads), not raw fetch latency. Regression byte-identical throughout.

**D24 — APA/DISCOVER gates demonstrated; two frozen criteria recorded as mis-specified
(2026-08-29).** `aie query apa|discover` implemented and evaluated against the withheld panel at
full scale. DISC-1: 100% recall (45/45) where the withheld gene overlaps no residual w1 gene;
82.7% where it does (claimed-by-design under the frozen claim-generously rule) — 86.2% overall,
MARGINAL as frozen. DISC-2: the frozen comparator was wrong (discover counts all molecules =
GeneFull semantics; first cut vs Gene read 76%); vs GeneFull: median rel err 26.3%, corr 0.9944 —
MARGINAL as frozen, boundary effects. APA-2: 72.0% in-window (MARGINAL); the aggregate-mass
criterion compared two different denominators (in-window mass vs net oracle deficit) and was not
scored — refined per-gene analysis: 87.9% of deficit≥5 genes show in-window mass, corr 0.328.
DISC-3: 33.9% of candidate mass hits withheld genes, 7.1% other v49; of the 59% "novel" mass, 57%
sits ≤5 kb downstream of a v49 3′ end (readthrough — the APA phenomenon again); 3,151 truly
intergenic loci remain. APA-1/APA-3 not yet run. Mis-specifications recorded, not silently
replaced. `memos/gate-apa-discover.md`.

**D25 — codec closed: zstd-19 kept, lz4 rejected with measurements; v1.1 accepted (2026-08-29).**
Sweep z12/z19/z22 at full scale: 128.76/121.41/121.32 MB, ingest 107/143/198 s, query latency
FLAT (~0.1 s at every level — open cost is dictionary decompression, chunk decode is invisible),
so lz4 would buy nothing observable and cost ~1.5–2× size; z22 buys 0.07% for +38% ingest. Stream
split (7→11 per chunk): the ~10 MB interleaving estimate was falsified (+0.5% incl. new indexes);
kept for column isolation, estimate-vs-measurement instance #6. Index enrichments: per-junction
genome-wide support totals in jpost (zero-decode junction stats), per-chunk max_anchor/n_cells
summaries, and a chunk_streams layout guard in meta — added after an 11-stream binary silently
mis-decoded a 7-stream archive mid-sweep (the invalidated first APA rerun). v1.1 = 121.41 MB,
11.57× vs CRAM 3.1, regression byte-identical. Real latency lever if ever needed: lazy dictionary
loading (backlogged). `memos/codec-indexing.md`.

**D26 — coding-stack audit: entropy layer exonerated; two exact factorizations found
(2026-08-29).** New read-only `aie debug` instrument (per-section/stream accounting, value-entropy
bounds, factorization sharing ratios). Findings: zstd-19 sits at or below the order-0 value bound
on every chunk stream (below on 5 of 11 — it exploits cross-value structure); chunk framing costs
+1.9% vs whole-stream compression; a rANS swap would harvest <1 MB — the entropy coder is DONE,
remaining bytes are modeling questions (strengthens D25). Factorizations measured: (a) **cell is a
pure function of the UMI class — 0 violations over 22.16M molecules** — so only the 87.0%
fresh-class molecules need cell entries: ≈3.5 MB exact win on the largest stream (28.8 MB);
(b) 58.0% of rep.pos offsets are 0 (first rep = anchor): ≈1.5–2 MB behind a flag bit in the unused
strand byte; (c) shapes collapse 6.89x into 117k junction chains but factored encoding nets only
~0.5 MB raw (semantic value > storage value); (d) pattern-alt vocabulary NEGATIVE (2.14x sharing,
1.14M unique offsets — break-even at best); (e) edges near-incompressible by construction.
Recommended: implement (a)+(b) under byte-identity regression (→ ~116 MB, ≈12.1x vs CRAM), then a
class-RLE experiment. `memos/coding-stack-audit.md`.

**D27 — cell-of-class + implied first-rep elision implemented; v1.2 accepted (2026-08-29).**
The two D26 factorizations, with one design correction forced by measurement. (1) Cell moved out
of the chunks into a global `cellofclass` section (one entry per class, fresh order — cell is a
pure function of the class): 28.8 MB in-chunk stream → 26.0 MB section. (2) First-rep elision v1
(flag bit in the layout byte) was a WASH: rep.pos −2.2 MB but the flag destroyed layout's run
structure (+2.6 MB) — a per-molecule flag costs what the zeros cost; elision only wins if the flag
is IMPLIED. v2: reps now store the span-minimum first (rows.rs; all consumers iterate reps
order-independently) and chains are position-sorted per molecule, making "first stored rep ==
anchor" an enforced invariant for every chain molecule — elided with no flag (rule: n_chains>0).
The writer bails if the invariant breaks (and did catch the contained-first rep order on the first
attempt). Chunks: 10 streams, layout back to pure strand bytes. Class RLE measured and REJECTED:
8.89 vs 9.09 MB (−0.20). Result: **113.73 MB, 12.35× vs CRAM 3.1** (was 121.41 / 11.57×),
regression byte-identical, query outputs identical. Cost: open 0.07→0.17 s (eager cellofclass
decompress) — lazy dictionary loading is the standing answer if it ever matters. Old 11-stream
archives (d0-z12/z19/z22 sweep artifacts) are refused loudly by the layout guard.

**D28 — D1 re-verification PASS: headline numbers reproduce at 4.9x scale (2026-08-29).**
pbmc_5k_v3, 383.9M reads, 5,038 cells, full pipeline (annotation-free ingest + v49 oracle in
parallel, CRAM 3.1 baseline, .aie v1.2 build, replay, byte-identity, fidelity). Results: 108.7M
molecules in 553.9 MB — **40.8 bits/molecule vs D0's 41.1 (format cost is scale-stable)**; 11.78x
vs CRAM 3.1 / 36.3x vs BAM / 48.1x vs FASTQ; regression BYTE-IDENTICAL; fidelity 0.307% mass
moved (D0 0.217%; deeper per-cell saturation, ~8x under the annotation signal), 10/5,038 cells
>1%; replay 187 s vs 1,225 s fresh (6.5x; D0's ~19x reflects thread counts — publish both pairs,
not one ratio). Measured cost that scales: open 0.17 -> 0.69 s (eager dict decode; cellofclass =
24.8% of file) — lazy dictionary loading is now a measured need, not a hunch. No number near STOP
territory. `memos/d1-reverification.md`, `results/gate1a-d1-v49.json`.

**D29 — multithreading + lazy loading accepted; v1.3 (2026-08-29).** rayon pool capped at 24
threads by default (shared host; RAYON_NUM_THREADS overrides). Parallelized: section compression
at write, chunk decompress+decode at load, replay's cell-sharded aggregate/filter/collapse (the
serial hashmap aggregation was the actual replay wall-clock), parallel flatten sort. Lazy loading
implemented as the v1.3 format change: cell-of-class in 65,536-class blocks, each its own zstd
section (coc.<b>, meta.coc_block guarded), queries open with meta+chroms+chunk-index only and
decompress just the coc blocks their classes touch; full loads decode blocks in parallel.
Measured: replay 49.9->6.4 s (D0), 181.9->25.3 s (D1) — **48x vs fresh STARsolo at D1, ~150x at
D0** (retires D28's thread-count caveat); ingest 795.6->609.7 s (D1; extraction now the
bottleneck — MT bgzf reader is the next lever); query open 0.17->0.00 s / 0.69->0.01 s, flat in
archive size. Price: +0.1-0.3% file size (independent coc frames). Regression byte-identical on
BOTH datasets. D1 audit: cellofclass (24.8%) is the top modeling target; chain/hop factoring
sharing IMPROVES with scale (8.2x/12.2x — revisit at D2+); pattern vocab negative again at scale.
Wart: query|head SIGPIPE panic, cosmetic, recorded. `memos/parallel-lazy-v13.md`.

**D30 — SIGPIPE fix, MT BGZF extraction, fairness audit (2026-08-29).** (1) SIGPIPE restored to
SIG_DFL in main() — `aie query | head` exits silently now. (2) noodles confirmed latest (bam
0.95.0 / bgzf 0.51.0); both extract_rows passes read via bgzf::io::MultithreadedReader (<=8
workers): extraction 93->65 s (D0), 520->389 s (D1); from-bam total 722->448 s; byte-identical
on both datasets and both paths. The remaining extraction time is the serial record-parse
consumer — the recorded fix is a batch pipeline with shard-local interning + merge, deferred
(ingest is paid once). (3) Fairness audit of the speed claim: D1 is same-budget (STAR 24 threads
vs rayon 24, same idle host) — 1,225 s vs 25.3-30.5 s; D0 is conservative in STAR's favor (32 vs
24). Caveats to publish with the walls: STARsolo's time includes writing the sorted BAM we
requested; and doing more work (alignment) is the claim, not a confound. (4) Parallel emit/
flatten: kept, within noise at current scale. `memos/parallel-lazy-v13.md` addendum.

**D31 — structural probes measured; rANS stage accepted (v1.4, 2026-08-30).** Meta-color-style
two-level factorings (Mac-dBG partial/meta colors) priced on real data: (a) pattern sub-set
sharing weak (20.4% LCP, 1.19x bigrams, ~5 MB at D1); (b) family absolute-site-set factoring
NEGATIVE (998k site-sets for 773k patterns, 0.77/set — anchor-relative interning is already the
right quotient; third pattern restructure closed); (c) locus-conditioned cell coding weak (Σ
H(cell|chunk) = 94% of H(cell); ambient barcodes drown the module structure). Generalization now
claimable: classes/patterns/chains/cell-of-class already ARE the archive's meta-colors; every
further indirection is priced at single-digit MB. rANS stage (5 noisy streams, static global
tables in `rans.tables`, chunk decode stays self-contained; meta.codec="rans1" guarded): D0
114.07 -> 112.89 MB (12.44x vs CRAM), D1 554.51 -> 554.41 (−0.02% — the null CONFIRMS the D26
"entropy layer is done" claim in-format). Cell-scoped backref recoding (−40% on backrefs)
REJECTED by explicit user decision: breaks independent chunk decode. Regression byte-identical
both datasets; D1 ingest 609.7 -> 443.8 s. Size axis closed as engineering frontier (~1–2%
measured headroom anywhere); remaining size strategy = framing (near-optimality, fidelity ladder,
capability-per-byte vs Malva, cross-sample dict amortization). `memos/structure-rans-v14.md`.

**D32 — capability campaign: velocity, discover→replay, APA-diff, federation (2026-08-30).**
Gates frozen first (docs/decisions/2026-08-30-capability-gates.md). VELO: full port of STARsolo
Velocyto semantics (alignToTranscriptMinOverlap TOL=6 + cell-scoped UMI intersection; 8 unit
tests); four measured iterations (no-merge 7.2/4.3/10.0 -> over-merge broke unspliced 10.0 ->
Gene-exact canon map -> empty-set transparency) landing at **1.81/0.52(PASS)/3.06% (D0),
2.67/0.86(PASS)/6.13% (D1) — MARGINAL as frozen**; residual = S<->A churn = the measured price of
2-rep evidence under velocity semantics (the D19 frontier under a harder statistic). DISCR:
**PASS/PASS — 45/45 withheld genes through discover→GTF→replay, median err 7.6% vs GeneFull**
(raw discover was 26.3%: the loop IMPROVES novel-locus quantification). APA-D: `apa --groups`
demo T-vs-mono with frozen null — **35/239 genes (15%) beyond null P95 (5% expected)**, top hits
PTPRC/CD44/LCP1 (correct biology). FED: `aie federate` — FED-1 PASS (byte-equal to single
queries), FED-2 PASS (0.66 s, 2 archives), FED-3: **44.2% of D0's junctions recur in D1** — the
atlas amortization argument with numbers. EM information-sharing: design recorded, multi-gene
discarded-mass pre-measurement owed before design freeze. `memos/capability-campaign.md`.

**D33 — EM upside bound + VELO-2 attribution (2026-08-30).** Multigene audit (new
`--audit-multigene`): 3.82% of read mass discarded for multi-gene ambiguity; **691k classes
(+6.25% of counted UMIs) are gene-informative yet fully discarded** — the EM design's funded
upside, plus 0.91M mixed classes counted with evidence ignored; 40.2% of rows match no v49 gene
(kept by the archive by design). VELO-2 (entry completeness strata via entries.complete.tsv):
2-rep substitution CONFIRMED as the driver of ambiguous-component disagreement (3.25% vs 1.16%)
but REFUTED for spliced (1.78 vs 1.41% — completeness-independent; suspects: class-level vs
read-level UMI correction, min-overlap edges) and inverted for unspliced (small single-read
molecules). D19's payload survives velocity better than feared. `memos/capability-campaign.md`
addendum.

**D34 — EM multimapper recovery: PASS both datasets; pooled sharing adopted (2026-08-30).**
Masked-evidence protocol (labeled truth from mixed classes, no simulation; gates frozen first):
uniform 39/35% -> cell-EM 94.4/93.4% -> **pooled cross-cell 98.2/98.3%** top-1 (D0 n=173k / D1
n=904k) — sharing cuts per-cell EM error ~75%; alpha- and seed-stable; pooled adopted. EM-2
impact: **+6.28% of the counted matrix recovered (688,558 classes), 85% at r>0.8**, recipients =
ribosomal/mito/histone paralog families (correct biology). EM counts are an additive opt-in
layer; byte-exact replay matrices untouched. EM-0 (STARsolo --soloMultiMappers EM replication)
deferred. `memos/em-recovery.md`, `aie em`.

**D35 — EM-0 PASS: STARsolo --soloMultiMappers EM replicated at 0.581% (2026-08-30).** Port of
STAR's per-cell multimapper EM (collapseUMIall: intersection candidate sets, skip-if-unique-seen
= our mixed/multi-only split, raw UMIs, 0.01 zero/convergence constants) over the byte-exact
Gene replay as the unique side; fresh v49mm oracle; filtered-cell UniqueAndMult-EM rel L1 =
0.581% [PASS <=1%]. Combined with D34: archive reproduces STARsolo's EM AND improves on it
(pooled sharing +4 pts top-1). Products: `aie em --star`, `aie em --emit` (594k-entry additive
layer). EM gates doc updated. D2 campaign launched (Parent_SC3v3_Human_Glioblastoma, 19.2 GB,
R1=28bp v3 verified; one path-nesting relaunch).

**D36 — D2 (glioblastoma) campaign complete (2026-08-30).** 250.7M reads, 5,573 cells, third
dataset / first non-blood tissue. REG byte-identical; FID PASS 0.426% (but 5.1% of cells >1% —
tissue-dependent concentration, un-decomposed); VELO **PASS** 1.59/0.45/1.97 (best yet — shallow
saturation = fewer lost middles, consistent with VELO-2); EM PASS pooled 95.7 vs cell 90.0.
SIZE: per-MOLECULE metric wrong invariant — 27.2 b/mol at 1.89 reads/mol vs 40.8 at 3.07; **per
read 14.4 vs 13.4 (+8%) — the stable invariant**; CRAM ratio 9.3x (∝ reads/molecule). SIGNAL
MARGINAL: tumor mass churn 2.12% < PBMC 2.43% (gene churn matches); the annotation-sensitive
premise needs a single-NUCLEUS brain dataset (D2' recorded); tumor cells bought discovery
richness, not churn — honest miss. D2-DISC not run, owed. `memos/d2-glioblastoma.md`.

**D37 — D2′ (brain nuclei) campaign complete: the premise gate PASSES (2026-08-31).** Brain_3p,
263.4M reads, 6,460 nuclei. **D2′-SIGNAL PASS: v32→v49 churn 4.90% = 2× PBMC — the
annotation-sensitivity premise demonstrated**; Gene counting uses only 34% of nuclei molecules
(the motivation in one number). REG byte-identical (4th dataset); FID PASS 0.454%; VELO PASS
under the frozen ≤5% bars (2.79/0.57/3.95 — nuclei invert to 48.1M unspliced vs 13.7M spliced
and unspliced is our best component); EM MARGINAL at the 95-bar but the finding is the gap:
**pooled 90.9% vs per-cell 59.4% (+31.5 pts)** — sharing matters most where cells are sparsest.
Size: 629.0 MB = 19.1 bits/read, 8.5x vs CRAM — per-read cost tracks complexity; coc audit:
7.7% ambient (aggregation tier minor), delta codec 15% above iid bound (inverts D0 verdict) →
per-block codec choice (~−25-30 MB) + transcript-relative second reps (~−10-15 MB) queued as
v1.5 AFTER campaigns idle; pattern restructure 4th negative. `memos/d2p-brain-nuclei.md`.

**D38 — v1.5: per-block cell-of-class codec choice; transcript-relative reps killed at design
(2026-08-31).** Each coc block now encoded both as delta-varints and as rANS over absolute cell
ids (6th global table; one codec byte per block; block independence and lazy queries preserved;
meta codec "rans2" guarded). Wins EVERYWHERE, byte-identical on all four datasets: D0 112.89 ->
111.08 MB, D1 554.41 -> 544.01, D2 449.42 -> 439.70, **D2' 628.99 -> 595.24 (−5.4%)** — the
per-dataset delta-vs-frequency verdict is now decided per block by measurement at write time.
Ranges: 11.3–18.1 bits/read; 9.0–12.7x vs CRAM. The transcript-relative second-rep lever was
KILLED at design: junction-free chains group same-UMI reads with potentially disjoint spans (no
shared transcript coordinate), and those are exactly the large rep.pos offsets; spliced chains'
second offsets are already bounded by the first block. UMI value population measured en route
(D0): 62.9% of 4^12 seen globally (skewed usage), 0.18% per-cell occupancy, <1% within-block
value collisions — MST-over-values evaluated and rejected: the value SET is ~free (2 MB bitmap);
the ASSIGNMENT is the cost and nothing downstream consumes it. `memos/structure-rans-v14.md`
lineage continues.

**D39 — D1 ambiguous-component deviation decomposed (2026-08-31).** The one velocity number
above 5% (D1 ambiguous, 6.13%) is fully attributed: entry-completeness strata on D1 show
ambiguous deviation of **1.49% on complete-evidence entries vs 6.50% on incomplete ones** (89%
of D1's counted mass is incomplete at 3.5 reads/molecule); spliced shows only a mild gradient
(2.04 vs 2.64 — the class-level-correction residual), unspliced inverts as at D0 (single-read
intronic molecules). Confirms the VELO-2 attribution at the depth where it matters: the residual
is the 2-rep payload's saturation-scaled price, purchasable via extra representatives for deep
chains (a frontier extension, not a fix). Paper's velocity section rewritten around this
decomposition; process language (gates/frozen/verdicts) removed from the paper per direction.

**D40 — application-support evaluations: prospective discovery, EM external concordance, EM
calibration (2026-08-31).** (1) PROSPECTIVE DISCOVERY (external truth = annotation history):
discover on real GENCODE v32 recovers **100% (52/52)** of expressed v49-added genes whose loci
are free of v32 features, median quant err 58% (~1.6x; genuinely novel lncRNA structures) —
overlapping additions absorbed by the conservative claim rule (38.4% of 315). (2) EM x
ALEVIN-FRY (disjoint information: sequence vs pooled patterns): top-254 EM-recipient genes,
unique-only median |log2| = 1.34 (2.5x) from alevin-fry -> **0.149 (11%) with the EM layer**,
corr 0.79->0.87, 79% of genes moved closer. (3) EM CALIBRATION: r>0.9 bucket = 97% of targets
at 99% empirical accuracy; mid-range (~2%) modestly overconfident, reported. All three folded
into the paper's application sections. Remaining queue: APA T-vs-mono replication in D1 +
annotated-3'-end concordance; D0∩D1 intergenic-locus replication; EM recovered-fraction D0<->D1
replication. Scripts 98/99. Paper title finalized: "Gravex: align once and query forever — a
minimal and sufficient index for annotation-independent single-cell RNA-seq".

**D41 — replication pack: all three applications replicate (2026-08-31).** (1) APA: independent
T-vs-mono contrast on D1 → **26/27 D0-significant genes replicate (chance: 5.3)**; 3'-end
concordance: 30/35 significants have both dominant sites within 200bp of annotated termini, with
differential usage manifesting as redistribution among shared ends (3/30 switch dominant end) —
the proportional mode of 3' regulation, stated as such. (2) Discovery: **79.4% of D0's 3,164
intergenic loci (89.8% of mass) recur in D1** (D1 at depth holds 14,851). (3) EM layer: per-gene
recovered fractions D0-vs-D1 **Pearson 0.992**, median |Δ| 0.0037 (n=2,823). Every application
now carries an external-truth + replication pair. Scripts 100/101; d1-em-layer emitted (3.15M
entries). Folded into the paper.

**D42 — Pool et al. (2023) reference-optimization delivered as a replay (2026-08-31).** Built a
Pool-style 3'-extended v49 (507,335 terminal exons, <=3 kb, neighbor-clipped) and replayed:
7 s (D0) / 35 s (D2'). **D0: +16.4% counted UMIs, 7,920 genes >10% (MBNL1, FOSB, DPYD). D2'
nuclei: +125.9% — counted molecules more than double — led by the canonical neuronal
long-3'UTR genes NRXN1 (40k->378k), NRXN3, RBFOX1 (5k->313k), KCNIP4, CADM2, GRIK2, RORA** —
exactly the biology Pool et al. document as lost. Caveat recorded: blanket per-isoform extension
also captures intronic mass near internal isoform ends on nuclei (blends their extension +
pre-mRNA fixes). Citation added; new capability subsection "Reference optimization as a replay";
their reanalysis call = our thesis, quoted in intro + discussion. Federation-extension downloads
(pbmc_10k_v3, 500_PBMC LT) in flight for the 5-archive atlas chapter.

**D43 — EM-DE probe (honest modest result) + PTPRC vignette (2026-08-30).** EM layer's impact on
T-vs-mono pseudobulk DE: 22/197 multimapper-heavy genes shift log2FC >0.5; aggregate agreement
with alevin-fry FCs unchanged (0.93 vs 0.96 median |Δ|; 54% closer = chance). Read as a SAFETY
property (consistent with r=0.992 gene-intrinsic fractions): the layer adds sensitivity without
reshaping contrasts — one sentence in the paper, not a capability claim. PTPRC vignette figure
built from both PBMC indexes (per-population 3'-site usage, replicated) and added to the APA
section. Atlas build chain (d3=pbmc_10k, d4=500_LT; whitelist match 95.5/93.1%, R1=28) running.

**D44 — six-index atlas built; support-weighted amortization measured (2026-08-30).** Two more
public PBMC samples ingested annotation-free: d3 = pbmc_10k_v3 (638.9M reads, 228.3M molecules,
1.02 GB = 12.8 b/read) and d4 = 500_PBMC LT Chromium X (75.7M reads, 99.0 MB = **10.5 b/read —
new low**; whitelist verified 93-96%). `aie federate` across all six indexes (3 tissues, v3→LT
chemistries, 20x cell range): **3.8 s, ~94k cells, byte-equal per-index**. Amortization: raw
junction-set overlap is depth-confounded (deep catalogues = singleton tails; 18-55% set-level),
but **support-weighted sharing is 89.2–96.5% per PBMC sample** — the shared core carries the
evidence; recurrence doubles as discovery confidence. ~120 GB FASTQ ↔ 2.8 GB of indexes.
Federation subsection rewritten with these numbers.

**D45 — genome signature + evidence-gated extension + APA statistics (2026-08-30).** User
directives: implement Pool-style extension; improve APA with optional genome access; stamp a
genome hash in the archive (blake3 per user). Shipped in COMBINE-lab/gravex:
(1) `evidence-io::genome`: per-contig BLAKE3 over uppercased bases (wrap/gzip-invariant) +
combined digest; stored in meta `genome_sig` (~8 KB, backward-compatible); per-contig lazy
verification. `ingest-archive --genome` (hashing overlaps extraction), `stamp-genome` (rewrites
meta only; all streams copied compressed byte-for-byte; d0 stamp = 10 s, replay stays
BYTE-IDENTICAL). Wrong genome refused (negative control verified).
(2) `query apa`: +cleavage column, internal-priming filter (>=12A/20nt or 8A-run/140nt
downstream, transcript-oriented), site x group G-test + permutation p. PTPRC: 100/1290 sites
flagged (one 625-UMI intronic site had a 20-A tract), p=4e-129.
(3) `query apa-test`: genome-wide per-gene multinomial G-test + BH-FDR; d0 = 11,166 genes in
14 s, 3,355 q<0.05 (T vs mono). Replaces the half-split P95 heuristic (memo P1).
(4) `extend`: evidence-gated per-gene 3' extension (corridor <=10 kb, clip at ANY same-strand
gene occupying the corridor — first cut clipped only at downstream *starts* and LOST 3% of UMIs
to induced multi-gene ambiguity; fixed — coverage-gap walk <=2 kb, >=5 UMIs/>=3 cells, IP
filter). D0: 5,048 genes, +2.3% UMIs. D2': 22,904 genes, **+20.2%** in 35 s.
KEY FINDING: gated extension shows the blanket-ext3p NRXN1-class explosions involve ~no movement
of true 3' ends (NRXN1 gated 1.00x) — that mass is intronic pre-mRNA captured by terminal-exon
inflation, i.e. the component of Pool's correction that GeneFull semantics addresses in a
principled way. Decomposition written into the paper.
(5) PolyASite 2.0 external validation: clean vs IP-flagged sites hit rate 16.9%/3.6% (all >=5
UMIs), 59.9%/11.4% at >=100 UMIs. herrmann2020polyasite added to refs.
Artifacts: runs/apa2/ (stamped d0/d2p copies, extgated GTFs, reports, sites-all.tsv, polyasite
bed). Canonical archives untouched (optimization agent benchmarks in flight); stamp them after
the optimization merge.

**D46 — Gravlax (2026-08-30).** "Gravex" was a typo: the name is **Gravlax** — cured salmon,
i.e. evidence preserved once and served forever; continues the lab fish lineage. Repo renamed
github.com/COMBINE-lab/gravlax (GitHub redirects the old URL); docs moved to
combine-lab.github.io/gravlax; logo wordmark, README, docs pages, paper title/text all updated.
Local clone now ${GRAVLAX_SRC}. Binary remains `aie`.

**D47 — visualization layer shipped (2026-08-30).** Per the viz memo's build order: (0) IGV
exports (`region --export-prefix`: per-strand molecule-coverage bedGraph + BED12 junctions);
(1) `region --plot` sashimi (strand-split coverage + support-weighted junction arcs + --gtf
underlay; arcs filtered to >=2 mols, top 80); (2) `apa --plot` lollipops (log stems, IP sites
hollow, per-group usage band; top 250 sites >=5 UMIs); (3) `em --plot` reliability diagram
(cal deciles per mode, exposed from em_experiment). All render via a zero-dep internal SVG
module (viz.rs/plots.rs, paper palette). Crate survey redux (user asked): adopted **resvg 0.48**
so `--plot x.png` rasterizes in-tool (DejaVuSans.ttf embedded — headless hosts have no fonts;
that bug ate the first PNGs). Rejected: kuva (no genome tracks), plotters (nothing we need),
charming/plotly/plotlars (JS-runtime/HTML interactive — the interactive path is igv.js around
our exports, backlogged). Examples in runs/apa2/viz/.

**D48 — optimization pass merged; archives stamped (2026-08-30).** Agent's performance pass
merged with D45/D47 features via git three-way (one conflict hunk). Headline (medians of 3):
replay-rows d0 8.4->2.6 s, d1 34.3->8.1 s (RSS 35->18.5 GB); em d0 285->12.9 s (22x), d1
30.2 min->64 s (28x); federate x6 4.24->0.48 s; query region d1 0.82->0.04 s. Byte-identity
re-verified post-merge (replay/EM/federate/apa/extend; cargo test 77/77). hashbrown 0.17
(foldhash) only where iteration order cannot reach output; EM/ingest stay rustc-hash for
determinism; replay aggregation rewritten hash-free (sorted runs + CSR). All six canonical
archives stamped with the genome signature (replay unchanged). Paper speed claims updated
(150-370x replay, 0.5 s federation, EM walls added). Full details:
memos/optimization-pass.md. INCIDENT: opt-branch commit swept docs/node_modules (56 MB
pagefind binary) into history via the baseline .gitignore; untracked in a follow-up commit —
history still carries the blobs; a force-push rewrite needs Rob's call.

**D49 — recovery EM packed and cell-sharded; memory gate PASS (2026-08-31).** The 38.1 GB D1
peak was representational rather than intrinsic: eager molecules + flat rows + per-row candidate
vectors + hash-owned class vectors + per-target responsibility vectors were simultaneously live.
Default `aie em` now streams batches into 8-byte support records, reduces 64 deterministic cell
shards, and iterates flat `u32` labels with `u64` CSR offsets; `--eager` preserves v1 as the
reference. Three-replicate medians: D0 3.10 s / 1.221 GB (masked), D1 10.06 s / 4.559 GB —
4.17--6.37x faster and 6.55--8.35x lower RSS than optimized v1. Impact-layer values, mass, and
all nonzero coordinates equal the frozen reference at emitted precision (relative L1 and moved
mass zero); 9/3 additional explicit `0.0000` coordinates reflect sub-rounding accumulation order.
Keyed `(seed,cell,class)` masking removes hash-layout dependence; pooled/blend scores remain within
0.1 point of v1 and pass the original gate. Code `3d7dcd5`; protocol/results scripts 133--135.

**D50 — bounded residual-site discovery selected; query/discovery optimization retained
(2026-08-31).** The conservative span rule's complete-denominator recall is 164/367 (44.7%).
The locked `residual-sites` channel retains those calls and adds transcript-incompatible evidence
clustered at transcript-oriented terminal bases. The predeclared 25/50/75/100-UMI grid selected
75: **266/367 (72.5%)**, 81.5% UMI-weighted recall, 26,787 D0 candidates (2.16x baseline), and
99.0% same-strand D1 recurrence. Recall is 51/52 clean, 52/60 antisense-only, 163/255
same-strand. All 14,398 added site intervals are <=1,001 bp; later-annotation confirmation is a
recall endpoint, not biological precision. Immutable query structures now avoid deep copies
(D1 chromosome-wide junctions: 0.51->0.44 s, 1.19->1.05 GB). Density-capped two-unit decode
makes default discovery 31.5%/29.2% faster on D0/D1 for +63/+81 MiB median RSS, with a hard
500,000-molecule pair cap and byte-identical output. Two replay rearrangements failed the locked
10% gate and were reverted. Code `c1ef9c4`; scripts 140--143.

**D51 — hierarchical EM exploration closed at untouched D4 (2026-08-31).** The registered
hard-group additive model first failed because whole-transcriptome normalization diluted its
80/20 priors; an exploratory 12,800/12,800 audit diagnosed rather than rescued that result.
Candidate-normalized convex pooling repaired the scale on D0, a posterior-mean Dirichlet proxy
found evidence-dependent gains only above 16 fitted unique-candidate units, and the monotone
depth hybrid `s(d)=d^2/(d^2+64^2)` repaired the shallow-depth loss. Its D0 verdict remained
MARGINAL because one registered ratio conflated model-form and group-information effects.
Before touching D4, commit `8b11be1` locked those as two orthogonal contrasts. Fresh v49
STARsolo fixed 565 called cells; replay-only, target-gene-excluded grouping selected 20 PCs,
resolution 0.70, seven groups (minimum 20 cells), and median pairwise Leiden-seed ARI 1.0.
On D4, hybrid-real beat proxy-real on all seeds by **0.460% mean NLL**, and real groups beat
shuffled groups on all seeds by **1.017%**. All top-1, depth, and memory gates passed; top-1
improved 0.018 percentage point and peak RSS was 1.09 GiB. The sole failure was mean Brier:
0.03021193 versus 0.03019876, an absolute regression of **0.00001318** (0.044% relative),
localized to the 4--16 and 16+ strata. Verdict: **MARGINAL**, without relabeling or tolerance
revision. The preregistration therefore does not open held-back D3. Real group signal has not
replicated, so the full hierarchical Dirichlet model is not licensed despite the structured
depth residual. Production emission remains pooled; the hierarchy branch is closed.

**D52 — shared cell/group query scopes and exact junction sets shipped (2026-08-31).** Before
implementation, the manifest fixed one strict scope contract for region, point-junction,
interval-junction, and batch queries, plus a new `jset` primitive. Cell and group maps reject
malformed, duplicate, empty, or unknown-barcode input. Junction sets place each UMI class in
exactly one of inclusion-only, exclusion-only, or shared support; conservative usage excludes
the explicit shared category. Synthetic coverage includes all three classes, absence, and strict
failures. On D0, point support equals exclusive plus shared support for every barcode, group rows
sum cell rows and bulk totals, and the same reductions reconcile on D1. Unscoped batch, region,
point, and enumeration outputs remain byte-identical to the pre-scope binary. The D0 union query
takes 0.05 s median versus 0.09 s for two point-query processes; grouped execution also takes
0.05 s and peaks near 73 MiB. The locked loci are sparse (11 informative D0 molecules; the
exclusion is absent in D1), so they are a correctness/performance fixture rather than a
differential-splicing result. Code `f1f99ec`; reproduction scripts 175--176.

**D53 — event-engine semantics and untouched descriptive gate locked (2026-08-31).** Before
implementation, `experiments/specs/event-engine-v1.md` fixed coordinate-only discovery for
alternative acceptors, alternative donors, and cassette structures; deterministic side/ID
rules; a single multi-event union decode; exact class-level include/exclude/both reduction;
strict overflow; annotation-neutral labels; and named-sample cohort federation without
pseudoreplicate p-values. The implementation/performance gates use D0/D1. D4 event output has
not been inspected: the predeclared AKR1A1 coordinate-defined cassette must have >=100
informative molecules and >=0.05 usage range among groups with >=20 informative molecules.
That is a descriptive structure gate, not replication of the post-hoc T-versus-monocyte effect.
Implementation `cc28ab4` now closes the decision as **PASS engine / FAIL held-out biology**. D0
chromosome 1 yields 1,905 coordinate candidates and 1,211 retained events in 0.201 s median,
reducing 4,263 posting decodes to 57; the first 50 exactly equal standalone junction-set totals
and group rows. D1 yields 17,951/9,242 in 1.64 s. D0/D1 cohort federation returns 1,291 shared
events in 0.80 s, exactly matching local reductions. The locked D4 AKR1A1 cassette is present
with 110 informative molecules, but eligible-group usage spans 0.04198 rather than the frozen
0.05. The earlier D0 junction set also lacks one cassette flank. The capability is published;
no biological-generalization claim is made. Scripts 177--178.

**D54 — federated event biology screen recovers external positive controls; FNBP1 confirmation
locked (2026-08-31).** After the D53 negative gate, a separately labeled post-hoc screen ran the
coordinate-defined event federation genome-wide. Four PBMC archives used the frozen marker rule
T = CD3E>0/LYZ=0 and M = LYZ>3/CD3E=0 for cell scope only. Of 20,028 events, 784 were testable
with conservative denominator >=10 per group in every dataset and had a strictly concordant
nonzero direction; 44 retained minimum |delta usage| >=0.10. FYB1 was the dominant locus: its
independently RT-PCR-validated cassette was 5.7--8.8% included in T cells and 98.3--99.1% in
monocytes across D0/D1/D3/D4; CD47 was another high-effect literature-positive locus. A
six-shard called-cell scan recovered MYL6, RPS24, and PPP1R12A compartment programs and exposed
an annotated in-frame 183-nt FNBP1 cassette at 94.1--100% inclusion in D2/D2p versus 4.0--13.1%
in four PBMCs. Exact STAR SJ aggregation agrees in direction but is a reducer audit of the same
alignments, not independent biological evidence. The scans took 43.78 s / 2.32 GiB and 73.71 s /
2.57 GiB as sequential chromosome processes. No confirmatory p-values are assigned. Before
opening another dataset, `experiments/manifests/federated-biological-screen.yaml` locked FNBP1,
the public 10x 2026 normal-brain Multiome GEX dataset, thresholds, and failure reporting. Scripts
179--181; compact results `post-v1-federated-biological-*`.

**D55 — prospectively locked FNBP1 normal-brain confirmation PASS (2026-08-31).** After D54
fixed the event, public dataset, thresholds, and failure rule, the 10x 2026 GEM-X Epi Multiome
normal-brain files were opened. Its filtered matrix defines 10,261 called nuclei. The fixed event
has 72 inclusion-only and 2 exclusion-only molecules (74 informative, 97.30% inclusion), passing
the >=50 and >=0.75 gates. Original barcode/UMI tags and cDNA sequence from 46,108 primary reads
at the locked locus were reconstructed and rerun through annotation-free two-pass STAR; the
resulting archive is exactly 72/2, while the independent STAR SJ matrix is 61 inclusion-right
versus 2 skip UMIs (96.83%), passing the >=0.50 and direction gates. Direct corrected-UMI
aggregation from the public Cell Ranger BAM is 62/2 (96.88%). All five locked criteria pass.
Qualification: selecting reads at the fixed window before STAR validates alignment direction and
reducer semantics but not whole-transcriptome mapping sensitivity. The supported statement is
prospective replication of the FNBP1 normal-brain cassette candidate, not novelty, cell-type
specificity, or mechanism. Scripts 182--185; result and paper figure
`results/post-v1-fnbp1-normal-brain-confirmation.json`.

**D56 — broad-cell localization is descriptive; federated execution is optimized
(2026-09-01).** A frozen marker-only rule assigned D5 nuclei before the fixed FNBP1 event was
inspected. The event reproduces 72 inclusion / 2 exclusion molecules. Both exclusions localize
to the four-molecule microglia/immune stratum, but one donor and four informative molecules do
not license a differential claim. Independently, exact row-denominator pushdown, packed
group-only reduction, bounded chromosome concurrency, and streaming summarization preserve
locked rows while reducing the legacy 394-MiB summarizer peak from 826,756 to 32,336 KiB.

**D57 — prospective eight-donor FNBP1 panel is underpowered (2026-09-01).** The gate and
event-blind feasibility record were remotely committed before acquiring 80 FASTQs from
GSE234790/PRJNA983239. Eight donor-isolated annotation-free STAR/Gravlax workflows processed
977,312,622 read pairs. Donors A--H contribute 4, 0, 1, 9, 6, 4, 13, and 2 informative
molecules; all 39 are inclusion-only and STAR-SJ agrees exactly. Only G reaches the registered
ten-molecule primary threshold and no donor has microglial informative evidence. The binding
verdict is **STOP_UNDERPOWERED** primary and **STOP_DESCRIPTIVE** microglia; thresholds are not
relaxed after inspection.

**D58 — exact molecular splice-graph primitive passes (2026-09-01).** Before implementation,
the v1 gate fixed strand-aware junction edges, exact UMI-class junction-set path fragments,
strict scopes, edge/path conservation, hard overflow, and explicit lower-bound/non-transcript
semantics. Gravlax `2d4d51e` passes all synthetic gates. Donor G yields 37 nodes, 22 edges, and
24 paths; three two-edge reverse-strand fragments carry six UMIs. Twenty runs take 0.004785 s
median and 13,204 KiB maximum RSS. This is a one-donor capability observation, not population
inference.

**D59 — replicate-aware graph layer passes engineering gate; FNBP1 STOP is preserved
(2026-09-01).** Registration commit `8d01b5a` fixed a strict unique-archive sample design, common
edge recurrence, complete zero-explicit sample matrices, sample-depth eligibility, and a
path-versus-rest beta-binomial contrast before implementation. Gravlax `173d925` passes balanced
null/signal controls and strict malformed-input, reference, low-depth, and conservation gates;
all 124 workspace tests and the documentation build pass. The eight-donor counts-only FNBP1 run
returns eight recurrent coordinates, ten strand edges, and 12 paths (86 fragments, including 11
on two-edge patterns) in 0.006347 s median and 14,844 KiB peak RSS. It emits no test. Only D, E,
and G reach ten reverse-strand fragments, so D57 remains binding. A population claim now requires
a prospectively locked, genuinely replicated two-condition dataset and sample-label calibration;
the full joint Dirichlet model remains reserved.
