# Memo — Visualization plan: optional `--plot` outputs and IGV export

**Date:** 2026-08-30 · **Status:** design only, no code changed. Scope: which `aie` commands
should grow an optional plot flag, what each plot shows, whether to use the `kuva` crate, and
whether emitting IGV-loadable tracks beats drawing our own.

---

## 1. The `kuva` crate — assessment

**What it is.** `kuva` 0.5.0 (crates.io, released 2026-08-07; MIT; repo
github.com/Psy-Fer/kuva — the author is a nanopore-bioinformatics developer, which shows in the
plot roster). A "scientific plotting library in Rust": builder-pattern plot structs → `Layout` →
`Scene` → backend. ~865 stars, 646 commits, ~14k downloads, actively developed through Aug 2026;
MSRV 1.87 (1.92 for PDF) — we are on 1.98, fine.

**Chart types.** 60+ plot structs, including several bioinformatics-native ones: `ManhattanPlot`,
`UpSetPlot`, `VennPlot`, `VolcanoPlot`, `PhyloTree`, `SyntenyPlot`, `Clustermap` +
`AnnotationTrack`, plus the full statistical standard kit (scatter/line/bar/box/violin/ridgeline/
heatmap/ECDF/QQ/ROC-PR/histogram/lollipop…). Also a CLI (`kuva scatter data.tsv -o plot.svg`)
that reads TSV/CSV/Parquet — worth knowing about independently of the library question: any
`--tsv` output we already emit can be piped into it today.

**Output.** SVG is the default, always-on backend; PNG and PDF behind `png` / `pdf` features;
a terminal backend for quick looks; optional embedded DejaVu font, LaTeX-ish math in labels.
SVG-first is exactly the shape we want (the paper pipeline is SVG-based).

**Dependency weight.** Required runtime deps: `chrono`, `colorous`, `image ^0.24`, `rand`,
`rand_distr`, `ryu`, `ttf-parser`. `image` as a *required* dep is the heavy one (pulls a pile of
codec crates); PNG/PDF features add `png`, `krilla`, `usvg`. Not enormous, but non-trivial for a
workspace that currently has ~10 direct deps and a codec-auditing ethos
(memos/coding-stack-audit.md).

**The decisive gap.** kuva has **no genome-track primitives**: no coverage track on a genomic
axis, no junction/sashimi arcs, no gene-model glyphs, no stacked per-locus tracks. Every plot we
actually most want (§2) is a genome-axis composition of rectangles, arcs, and rugs — precisely
what kuva does not provide and what is trivial to draw directly (the paper's figures are already
generated as styled SVG via `paper/figs/gen_figures.py`; the aesthetic is defined).

**Alternatives.**
- **`plotters`** — the mature general-purpose choice; SVG backend exists; verbose API, raster
  bias, and again zero genomics primitives. No advantage over kuva for us.
- **`charming`** — ECharts wrapper; static rendering drags in a JS engine. Wrong shape for a CLI.
- **Direct SVG generation** — a small internal helper (`viz.rs`: viewbox, path, rect, text,
  cubic-Bezier arc, axis ticks ≈ 250–350 lines, zero deps). Full control of the genome-axis
  layout, styles matched to the paper palette, output diffs cleanly in git.

**Recommendation: direct SVG for the genome-native plots; no kuva dependency now.** The 2–4
plots worth building first are all either genome-track compositions (kuva can't) or trivial
bars/scatter (kuva overkill). Revisit kuva as an optional `viz-kuva` cargo feature only if we
later want an UpSet panel for federation (§2.4) — its `UpSetPlot` is the one type that would
save real work. Interim: kuva's *CLI* on our existing `--tsv` outputs covers exploratory needs
with zero code.

---

## 2. Per-command plot designs

Flag convention: `--plot <out.svg>` everywhere; plot generation must never change stdout output
or counting semantics (plots are a side product of the same decode pass). All plots: white
background, paper palette (`#1a6faf` blue / `#c98f1a` orange / grays), genomic x-axis with
comma-grouped tick labels and a kb scale bar.

### 2.1 `query region --plot` — junction-arc / sashimi-style locus view  ★ flagship

**Data already in hand** (querycmd.rs `Region` arm): the decode pass touches every molecule in
the window with anchor, strand, cell, UMI class, junction chains (rep positions + shape ids →
exon-block intervals and donor/acceptor pairs via `has_junction`-style block walking), and mm
placements. Today we only accumulate per-cell counts; coverage and junction support are one
extra accumulation in the same loop:
- per-bp span coverage: for each class-deduplicated molecule, add its rep block intervals into a
  binned depth array (~1–2k bins across the window);
- junction support: `FxHashMap<(donor,acceptor,strand), FxHashSet<(cell,class)>>`.

**Layout.**
- Track 1 (top, ~40% height): coverage area chart, + strand filled upward in blue, − strand
  mirrored downward in orange (or overlaid at 50% alpha when one strand dominates); y = molecule
  depth, linear with auto-scale, y-axis label "molecules".
- Track 2 (middle, ~45%): junction arcs — cubic Bezier from donor to acceptor, arc height ∝
  genomic span (rank-scaled so nested introns stay separable), stroke width ∝ log2(UMI support),
  + strand arcs above the baseline, − below; UMI count labeled at each apex for the top ~15 arcs,
  faint gray for the singleton tail.
- Track 3 (bottom): axis + optional gene models if a `--gtf` is passed (thin exon boxes — makes
  the "annotation drawn *under* annotation-free evidence" point visually).
- `--groups <tsv>` (add, mirroring `apa`): small multiples — one coverage+arc row per cell
  group, shared x-axis, per-group UMI-normalized. This is the per-cell-resolution picture no
  coverage-domain resource (Snaptron/recount) can draw.

**Effort:** medium — 2–3 days. Junction extraction + coverage accumulation ~80 lines in the
existing loop; SVG track code ~250 lines, reusable by 2.3.

### 2.2 `query apa --plot` — 3′-site usage on a genomic axis  ★ cheapest win

**Data already in hand:** the `out: Vec<Site>` vector — (start, end, strand, UMIs, cells,
per-group UMI counts). Nothing new to compute; the plot is a renderer over `out`.

**Layout.**
- Genomic x-axis over the window. Each site: a lollipop at the site midpoint — stem height ∝
  UMIs (log scale), head colored by strand (blue +, orange −), thin underline marking the site's
  clustered extent (start–end), top ~8 sites labeled with coordinate and UMI count.
- With `--groups`: replace each lollipop head with a compact per-group bar pair (counts
  normalized per group's window total), and add a bottom "proportion band": for each site a
  stacked 100% bar of group shares — differential 3′ usage between cell populations at a glance
  (the PTPRC/fig_ptprc story as a per-locus product, straight from the archive).

**Effort:** small — 1 day. Pure rendering; data structure is final.

### 2.3 `query discover --plot-dir <dir>` — locus portraits for candidate loci

**Data:** candidates (chrom, start, end, strand, UMIs, cells) exist; the *member molecules*
(spans) live transiently in `runbuf` during clustering. For a portrait we additionally want the
members' junctions and 3′ ends — a small extension of the buffered tuple (or a second targeted
decode of each candidate's window reusing the 2.1 machinery, which is simpler and touches only
`n_candidates` chunks).

**Layout** (one SVG per top-k candidate, `portrait_<chrom>_<start>.svg`, default k=12):
- Header line: locus, strand, UMIs, cells (the evidence summary already printed).
- Track 1: unclaimed-molecule span coverage over locus ± 2 kb flank.
- Track 2: junction arcs among member molecules — spliced structure is the strongest "this is a
  real gene, not noise" evidence for an unannotated locus.
- Track 3: 3′-end rug + 24-bp-clustered site ticks (a focused terminator = real transcription
  unit; a smear = artifact).
- Footer: nearest annotated gene and distance (annotation available in this command via `--gtf`).
This is the reviewer-facing artifact of the discover→replay loop: "here is the evidence" as one
picture per AIENOVEL id.

**Effort:** medium — 1–2 days *after* 2.1 exists (reuses all three track renderers).

### 2.4 `federate --plot` — cross-sample junction recurrence

**Data today:** per-archive (UMIs, cells, index support) for **one** junction — enough for a
per-sample bar panel, not for an UpSet. The paper's recurrence claim (89–97% of support mass on
shared junctions) came from script-side catalogue sweeps.
- **Minimal (small):** two aligned horizontal bar rows per archive — UMIs and cells — plus index
  support as a faint background bar; archives sorted by UMIs. Honest and cheap.
- **Full (medium-large):** a `federate --region chrom:s-e` mode enumerating all catalogued
  junctions in a window per archive → junction × sample support **heatmap** (rows = junctions
  sorted by position, columns = samples, cell shade = log UMIs) with a presence-pattern UpSet
  strip on top. This needs a new range scan over the junction catalogue per archive — a real
  query-path addition, and the one place kuva's `UpSetPlot` would pay for itself.

**Effort:** bars small (0.5 day); region-recurrence heatmap medium (2–3 days incl. the new
query path). Defer the full version.

### 2.5 `em --plot` — calibration diagnostics

**Data already computed** in the EM-1 scoring pass: per masked target, max responsibility and
correct/incorrect; per mode, top-1/expected accuracy (the memo table). The paper already quotes
responsibility buckets ("r > 0.9 bucket 99% correct") — this plot makes the instrument
self-reporting.

**Layout** (two panels side by side, ~900×360):
- Left — reliability diagram: bin masked targets by max responsibility (10 bins), x = mean
  predicted r, y = empirical top-1 accuracy, y=x reference diagonal, point area ∝ log bucket
  count (the r>0.9 bucket dominates; the overconfident mid-range buckets become visible instead
  of a prose caveat). One series per mode (pooled solid, cell-only faint).
- Right — grouped bars: top-1 and expected accuracy for uniform / cell-only / pooled / blend,
  the memo's table as a figure. With `--mask 0 --emit`, substitute a histogram of
  responsibilities on the real recovered layer (the "85% at r>0.8" claim).

**Effort:** small — 1 day. Plain scatter + bars; no genome axis.

### 2.6 `ingest-archive --plot` — knee / QC at ingest

**Data:** the writer sees every (cell, class); a per-cell class-count accumulator is one map.
Caveat: ingest is already whitelist-filtered, so the classic ambient-vs-cell knee is truncated —
diagnostic value is moderate. The per-stream byte accounting (footer; also `aie debug`) is the
more distinctive QC surface.
**Layout:** (a) barcode-rank knee, log-log, rank vs class-deduplicated UMIs; (b) per-stream
bytes as a horizontal bar stack (anchor/class/layout/… + dictionaries + indexes) — the Gate-0
picture, self-drawn. **Effort:** small, but **lowest priority** — the paper's fig_bits already
covers (b), and (a) is weakened by the whitelist.

---

## 3. Export for existing browsers (IGV) — complement, and partly done

Sometimes emitting tracks beats drawing: IGV gives free zoom/pan, gene annotation context, and
BAM juxtaposition. The archive's decode passes already have everything needed, and one export is
**already implemented**: `discover --emit-gtf` loads in IGV as a gene track today.

Proposed `--export-prefix <p>` on `query region` (same decode pass as `--plot`):
- `<p>.cov.plus.bedGraph` / `<p>.cov.minus.bedGraph` — span coverage per strand; users can run
  `bedGraphToBigWig` themselves (do not add a bigWig writer dep; bedGraph loads fine at locus
  scale).
- `<p>.junctions.bed` — TopHat/RegTools-style junctions BED12 (score = UMI support, blocks =
  1-bp anchors at donor/acceptor): IGV renders these natively as sashimi-style arcs.
- With `--groups`: one bedGraph *set per group* — per-population tracks stacked in IGV, which is
  exactly the comparison IGV cannot make from a BAM.
- `query apa --bed` — BED6 of sites, score = UMIs, name = `site_<n>_<umis>u_<cells>c` (the
  `--tsv` output is 90% of the way there).

**Effort:** ~0.5 day total; highest value-per-line in this memo. Recommendation: ship the
exports *with* the first `--plot` flags — they de-risk the plots (any doubt about our drawing is
checkable against IGV over the identical data) and serve users who live in IGV anyway. Our own
SVGs remain necessary where IGV stops: per-cell-group small multiples, per-sample federation
panels, calibration diagnostics, and publication-styled locus portraits.

---

## 4. Prioritized recommendation — build these first

| # | What | Why first | Effort |
|---|---|---|---|
| 0 | `query region --export-prefix` + `apa --bed` (IGV tracks) | near-free; validates plots against IGV; immediately useful | 0.5 d |
| 1 | `query region --plot` sashimi (coverage + arcs + optional GTF underlay, `--groups` small multiples) | the flagship image of annotation-independence; renderer reused by discover portraits | 2–3 d |
| 2 | `query apa --plot` per-group 3′ usage | unique capability ("the matrix cannot represent this") made visible; data fully computed | 1 d |
| 3 | `em --plot` reliability diagram + mode bars | turns the calibration *claim* into a self-generated instrument readout | 1 d |

Then, in order: discover locus portraits (reuse #1's tracks), federate per-sample bars, and only
after a decided need: federate region-recurrence heatmap/UpSet (where kuva as an optional
feature becomes worth revisiting). Skip `ingest-archive` knee for now.

Implementation shape: one new `viz.rs` in `crates/aie` (or a tiny `viz` crate), zero new
dependencies, SVG only, styles copied from `paper/figs/gen_figures.py`'s palette so CLI output
and paper figures are visually one family.
