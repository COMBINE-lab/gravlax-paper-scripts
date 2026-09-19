# Gravlax — experimental & provenance manifest

Written 2026-08-30 on host **analysis-host**. This document is self-contained: it is
written for a reader (human or agent) with **no prior context**, and describes every experiment
run for the Gravlax project — where the data live on this host, how the tools were run, the
current results, and the pre-registered gates each result was judged against.

Path convention: `$P` = `${GRAVLAX_PROJECT_ROOT}` (the project
root). All relative paths below are relative to `$P` unless absolute.

---

## 1. What this project is

**Gravlax** builds a compact, molecule-resolved, annotation-independent evidence index (`.aie`)
for 10x Chromium 3′ scRNA-seq, from raw reads plus a genome — no annotation is consulted at
build time. Supplying any GTF *at query time* reproduces, byte-for-byte at matrix level, the
cell×gene matrices (Gene / GeneFull / Velocyto semantics) that fresh STARsolo processing of the
reads would produce under that annotation. The index also supports indexed region/junction/
3′-site queries, unannotated-locus discovery with re-quantification, evidence-gated annotation
extension, pooled cross-cell EM multimapper recovery, and federated queries across N indexes.

- Paper (typst): repo `github.com/COMBINE-lab/moleceular-evidence-store-paper`,
  local clone `${GRAVLAX_PAPER_SRC}`, entry `main.typ`.
  Title: *"Gravlax: align once and query forever — a minimal and sufficient index for
  annotation-independent single-cell RNA-seq."*
- Code (Rust workspace, binary `aie`): repo `github.com/COMBINE-lab/gravlax`
  (renamed from `gravex` 2026-08-30; GitHub redirects), local clone `${GRAVLAX_SRC}`
  at commit `d9257ca`. `$P/src` mirrors the repo (synced after each merge).
- Docs site: `https://combine-lab.github.io/gravlax/` (Astro Starlight in `docs/` of the code
  repo, deployed by `.github/workflows/docs.yml`).
- Predecessor project (STOPped; conventions reused): `${GRAVLAX_PYTHON_ROOT}`.

**Authoritative running records (read these for full detail):**

| Record | Path | Content |
|---|---|---|
| Decision log | `docs/decision_log.md` | Numbered decisions **D1–D48**, chronological, every result and incident |
| Frozen gates | `docs/decisions/*.md` | Pre-registered pass/stop thresholds, dated, written **before** measurement |
| Gate & experiment memos | `memos/*.md` | One memo per gate/campaign with verdict tables |
| Machine manifest (week 1) | `experiments/manifests/gates.yaml` | Hardware, toolchain, references, D0, baselines (partially superseded by this file) |
| Format spec | `docs/format-spec.md` | The `.aie` on-disk format |
| Original plan | `docs/project_plan.md`, `docs/implementation_plan.md` | The pre-registered plan this all executed |

---

## 2. Host and environment

- Host: analysis-host — 2× AMD EPYC 9575F (256 logical cores), 1.5 TiB RAM, RHEL 8.
  Shared interactive box, **no scheduler**; project policy caps runs at 24–32 threads.
- **Every command below requires `source $P/env/setup.sh` first.** It sets project-local
  `CARGO_HOME`/`RUSTUP_HOME`/`TMPDIR`/`CONDA_PKGS_DIRS` (the home NFS is nearly full and
  `~/.cargo`, `~/.rustup` are dangling symlinks on this host), puts `env/sc/bin` on PATH, and
  exports `AIE_THREADS` (default 32). `AIE_CHROM=chr19` restricts the whole pipeline to one
  gene-dense chromosome for minute-scale iteration; unset for whole-genome.
- Toolchain: rustc/cargo 1.98.0; STAR **2.7.11b** in `env/sc` (conda; the host-wide STAR 2.7.2b
  lacks GeneFull/Velocyto); samtools 1.23.1; gffread 0.12.9; python via
  `${GRAVLAX_PYTHON_ROOT}/env/py312` (3.12; numpy/scipy/polars/pysam/
  matplotlib — analysis and figures only, never in the measurement path); typst vendored at
  `<paper-repo>/env/bin/typst` and `$P/env/bin/typst`.
- Node 22.18.0 for the docs site at `${GRAVLAX_TOOLS}/node-v22.18.0-linux-x64/`.
- Build the tool: `cd $P/src && cargo build --release` → `$P/src/target/release/aie`
  (identically buildable from `${GRAVLAX_SRC}`). `cargo test --release`: 77 tests.

## 3. Reference data (all digests in `data/DIGESTS.sha256` unless noted)

| Item | Path | Provenance |
|---|---|---|
| Genome | `data/GRCh38.primary_assembly.genome.fa.gz` | GENCODE/EBI, sha256 `3023c070…5452`; also at `/shared-resources/references/GRCh38.gencode.v49/` |
| Genome blake3 signature | stamped inside every archive (`meta.genome_sig`) | combined digest `2817ea0bc4a2b919…`, algo `aie-genome-blake3-v1`, 194 contigs |
| 10x v3 whitelist | `data/3M-february-2018.txt` | sha256 `843a6f70…819b` |
| GENCODE v49 (A1, baseline) | `annotations/gencode.v49.annotation.gtf[.gz]` | EBI release_49; 78,691 genes / 507,365 transcripts |
| GENCODE v32 (A2-far) | `annotations/gencode.v32.annotation.gtf[.gz]` | headline cross-annotation case (CellRanger 2020-A era) |
| GENCODE v48 (A2-near) | `annotations/gencode.v48.annotation.gtf[.gz]` | adjacent-release hard case |
| Withheld panel (w1) | `annotations/withheld/w1.gtf` + `w1.manifest.tsv` | v49 minus 900 real expressed features (300 genes / 300 isoforms / 300 3′-trims), built by `scripts/70_withhold_annotation.py`, seeded |
| Pool-style 3′-extension GTF | `annotations/gencode.v49.ext3p.gtf` | every terminal exon +≤3 kb, neighbor-clipped; built by `scripts/102_pool_extension.py` |
| Evidence-gated extension GTFs | `runs/apa2/gencode.v49.extgated*.gtf` | built by `aie extend` (§8.11) |
| STAR indexes | `runs/star-index-base` (annotation-free base; sjdb inserted at map time), `runs/star-index-{32,48,49}` | see `scripts/05_build_anno_indices.sh` |
| PolyASite 2.0 atlas | `runs/apa2/polyasite.bed.gz` | polyasite.unibas.ch, GRCh38.96 clusters (downloaded 2026-08-30) |

## 4. Datasets

All are public 10x Genomics demonstration datasets (3′ chemistry; R1 = 16 bp CB + 12 bp UMI).
FASTQs live under `data/`; identity is by the 10x dataset name (the original tarballs are kept
alongside the unpacked FASTQs; `DL_STATUS` files record successful download+untar).

| ID | 10x dataset | Reads | Cells (called) | FASTQ location |
|---|---|---|---|---|
| D0 | `pbmc_1k_v3` | 66.6 M | 1,225 | `data/fastq/` (sha256s in DIGESTS) |
| D1 | `5k_pbmc_v3` | 383.9 M | 5,038 | `data/d1/5k_pbmc_v3_fastqs/` |
| D2 | `Parent_SC3v3_Human_Glioblastoma` | 250.7 M | 5,573 | `data/d2/Parent_SC3v3_Human_Glioblastoma_fastqs/` |
| D2′ | `Brain_3p` (single-**nucleus** brain) | 263.4 M | 6,460 | `data/d2p/Brain_3p_fastqs/` |
| D3 | `pbmc_10k_v3` | 638.9 M | not called² | `data/d3/pbmc_10k_v3_fastqs/` (atlas demo only) |
| D4 | `500_PBMC_3p_LT_Chromium_X` | 75.7 M | not called² | `data/d4/500_PBMC_3p_LT_Chromium_X_fastqs/` |

D0–D2′ carry the full validation battery; D3/D4 were ingested (annotation-free) for the
six-index federation/atlas demonstration only. ² No STARsolo oracle / cell calling was run for
D3/D4 — their archives carry whitelist-corrected barcodes without an EmptyDrops call.

## 5. Oracle: STARsolo, and how it was run

The fidelity oracle is STARsolo 2.7.11b. Two alignment configurations exist:

1. **Oracle runs** (`scripts/10_oracle_starsolo.sh`): one shared annotation-free base index;
   the annotation enters only via `--sjdbGTFfile <gtf>` at mapping time. Key flags:
   `--soloType CB_UMI_Simple --soloFeatures Gene GeneFull SJ Velocyto
   --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR
   --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter EmptyDrops_CR`.
   Outputs: `runs/oracle/{dev,full,d1,d2,d2p}/<annotation>/Solo.out/...` (`full` = whole-genome
   D0; annotations v32/v48/v49/w1, plus `v49mm` for the multimapper study).
2. **Ingest runs** (`scripts/40_ingest_align.sh`): annotation-free two-pass alignment
   (`--twopassMode Basic`, no GTF, `--soloFeatures SJ`, raw `CR/UR` tags only, secondaries
   kept). These BAMs (under `runs/ingest/`) are the sole input to `aie ingest-archive` —
   the tool performs its own whitelist correction (STARsolo's exact
   1MM_multi_Nbase_pseudocounts algorithm, reimplemented) and its own UMI handling.

Independent sanity oracle: a validated simpleaf/alevin-fry run of D0 at
`/fs/cbcb-lab/rob/pbmc_1k_v3_simpleaf/af_quant/alevin/` (15,519,436 UMIs / 71,123 barcodes),
also used for the EM concordance study (§8.9). Storage baseline: Malva's index of the same D0
(`${MALVA_ROOT}/indices/pbmc_1k_v3`, 375 MB).

## 6. The `.aie` archives (current, format v1.5, all genome-stamped)

Built by `aie ingest-archive <ingest.bam> --whitelist data/3M-february-2018.txt --out <x>.aie`
(see `scripts/40_ingest_align.sh`, `96_d2_campaign.sh`, `97_d2p_campaign.sh`,
`103_atlas_build.sh` for the per-dataset invocations). Format spec: `docs/format-spec.md`.

| Archive | Path | Bytes | bits/read | vs CRAM 3.1 |
|---|---|---|---|---|
| D0 | `runs/archive/d0.aie` | 111,087,381¹ | 13.3 | 12.7× |
| D1 | `runs/archive/d1.aie` | 544,022,039¹ | 11.3 | 12.0× |
| D2 | `runs/archive/d2.aie` | 439,709,823¹ | 14.0 | 9.5× |
| D2′ | `runs/archive/d2p.aie` | 595,246,459¹ | 18.1 | 9.0× |
| D3 | `runs/archive/d3.aie` | 1,020,831,458¹ | 12.8 | — (atlas only) |
| D4 | `runs/archive/d4.aie` | 99,050,193¹ | 10.5 | — (atlas only) |

¹ Sizes listed are pre-stamp; the genome signature added ~8 KB to each on 2026-08-30 (`aie
stamp-genome`; stamping copies all evidence streams compressed, byte-for-byte, and replayed
matrices were verified unchanged). CRAM baselines (all tags retained): `runs/gatee/`.
Superseded archives kept for the codec sweep record: `d0-z{12,19,22}.aie`, `d0-c{4,16}mb.aie`,
`dev.aie` (chr19).

**The standing correctness contract (held through every format revision v0→v1.5):**
`aie replay-rows <archive> --gtf G ...` and `aie replay-rows <ingest.bam> --from-bam
--whitelist W --gtf G ...` produce **byte-identical** matrices. Reference outputs for the
current format: `runs/replay/{d0,d1,d2,d2p}-v15/` and `*-frombam/`.

## 7. Gate methodology

Every major claim passed through a **pre-registered gate**: thresholds were frozen in a dated
file under `docs/decisions/` *before* the measurement ran, and each gate ended in a memo with a
verdict table (PASS / MARGINAL / STOP). The frozen-gate files are:
`2026-08-28-frozen-gates.md` (G-B/C/D/E), `2026-08-28-archive-entropy-budget.md`,
`2026-08-29-apa-discover-gates.md`, `2026-08-30-capability-gates.md`, `2026-08-30-d2-gates.md`,
`2026-08-30-em-gates.md`. A hard rule learned early (D18): **never trust dev-scale (chr19)
numbers for frontier decisions** — the 17%-saturation subsample flipped two G-D conclusions.

## 8. Experiment catalog

Each entry: what was asked → how it was run → where results live → verdict.

### 8.1 Week-1 gates (D0, whole genome; memos `gate-{B,C,D,E}-decision.md`)
- **G-B, annotation divergence** (`scripts/30_gateb_compare.py`): does annotation churn move the
  matrix enough to be worth replaying? v32→v49 moves **2.43%** of UMI mass but changes **21.8%**
  of expressed genes by >10%; v48→v49 only 0.369%. **GO** (capability matters for *distant*
  pairs).
- **G-C, alignment invariance** (`scripts/20_gatec_align_configs.sh`, `aie gate-c`): is
  annotation-free alignment a fidelity ceiling? Identical evidence tuples for **98.87%** (vs
  v32-sjdb) / **98.70%** (vs v49) of gene-assigned reads; two-pass essential (93.52→98.36%).
  **MARGINAL** (bar was 99.0), residual = splice-junction placement, matrix footprint ~0.15 pp.
- **G-D, annotation-independent UMI grouping** (`aie gate-d`): locus-based grouping error
  0.939% median per-cell L1 — later superseded by the **UMI adjacency-graph representation**
  (memo `umi-adjacency-graph.md`): store classes+1MM edges, not values; 1.20 bits/molecule and
  error **0.232%** (merging deferred to replay). **PASS** (error ≪ 2.43% signal).
- **G-E, storage** (`scripts/50_gatee_baselines.sh`): first real encoder hit 119.5 MB =
  **11.76× smaller than CRAM 3.1** (all tags), 3.29× smaller than Malva. **PASS.**

### 8.2 Replay fidelity, Gates −1A/−1B/−1C (memos `gate-1A-replay.md`, `gate-1BC-decision.md`, `gate-1C-withheld-first-cut.md`)
`aie replay` / `aie replay-rows` reproduce fresh STARsolo to **0.217%** of UMI mass on D0
(0/1,225 cells off by >1%). One archive replayed under v32/v48/v49/w1 without rebuild:
0.207/0.212/0.217/0.246% mass moved. Withheld-annotation recovery (w1 = v49 minus 900 real
features): **99.88% mass recall, 300/300 withheld genes within 10%**. Key mechanism discovered
en route (D16): STARsolo counts unique-gene multimappers via secondary alignments (2.2% of
matrix mass), so alternative placements are mandatory archive evidence. Cross-dataset fidelity
(replay deviation vs v49 oracle): D0 0.217%, D1 0.307%, D2 0.426%, D2′ 0.454% — always ≥5×
under that dataset's annotation signal (2.43/2.43/2.12/4.90%).

### 8.3 Format evolution v0 → v1.5 (memos `archive-v0.md`, `codec-indexing.md`, `coding-stack-audit.md`, `parallel-lazy-v13.md`, `structure-rans-v14.md`; log D18–D38)
Measurement-gated compression ladder; every step kept the byte-identity contract. Landmarks:
2-rep span-extreme chains (1-rep rejected at +0.37 pp fidelity), paralog-pattern dictionary
(7.85× sharing), global (cell,value) UMI classes (per-chrom version over-counted 2.9% — D20
incident), cell=f(class) elision, zstd-19 codec sweep (lz4 rejected with numbers), 10-stream
chunks with meta guards (`chunk_streams=10`, `coc_block=65536`, `codec="rans2"` — guards exist
because an 11-stream binary silently mis-decoded a 7-stream archive, D25 incident), static-table
rANS on the 5 memoryless streams, v1.5 per-block cell-of-class codec choice (delta-varint vs
rANS-absolute per 65,536-class block). Negative results are recorded too: transcript-relative
second reps, MST-over-UMI-values, pattern-alt vocabulary, site-ordered positions, class RLE,
cross-sample junction catalogue (`junction-catalogue-negative.md`: it *lowers* G-C agreement).

### 8.4 Capability gates (frozen `2026-08-30-capability-gates.md`; memo `capability-campaign.md`)
- **Velocity replay** (`scripts/91_velo_gate.py`, `94_velo2_strata.py`; `aie replay-rows
  --velocity`): spliced/unspliced/ambiguous deviations 1.8/0.5/3.1% (D0), 2.7/0.9/6.1% (D1),
  1.6/0.5/2.0% (D2), 2.8/0.6/4.0% (D2′). D1 ambiguous decomposed (D39): complete-evidence
  molecules 1.49% vs 2-rep-incomplete 6.50% — the saturation-scaled price of 2 representatives.
- **Discover → replay** (`aie query discover --emit-gtf` + replay; `scripts/92_discr_gate.py`,
  `81_disc_apa_v2.py`, `82_disc3_classify.py`): 45/45 well-posed withheld genes recovered;
  quantification via the discover→GTF→replay loop at 7.6% error. **Prospective test**
  (`scripts/98_prospective.py`, D40): discovery run against *v32* recovers **100% (52/52)** of
  clean-stratum genes that v49 later added.
- **APA v1** (`scripts/80_apa_gate.py`, `93_apa_diff.py`, `100_apa_replicate.py`,
  `105_ptprc.py`): differential 3′ usage T vs monocytes; 26/27 D0-significant genes replicate
  in D1; 30/35 dominant sites within 200 bp of annotated termini; PTPRC vignette in both
  datasets. Group files: `runs/replay/apa-groups-{d0,d1}.tsv` (+ null files).
- **Federation** (`aie federate`): per-archive byte-equal counts, six archives.

### 8.5 EM multimapper recovery (frozen `2026-08-30-em-gates.md`; memos `em-recovery.md`; log D33–D35)
Masked-evidence protocol (`aie em --mask 0.2`): mixed classes with unique truth are masked and
must be recovered from multi evidence alone. Top-1 recovery, uniform → per-cell → **pooled**:
39.1→94.4→**98.2%** (D0), 35.4→93.4→**98.3%** (D1), 40.0→90.0→95.7% (D2), 33.5→59.4→**90.9%**
(D2′ — sharing matters most where cells are sparsest). Calibration is decile-tracked (see
`aie em --plot`). **EM-0** (`scripts/95_em0_gate.py`): STARsolo `--soloMultiMappers EM`
replicated to 0.581% (`aie em --star`). Impact layer: `aie em --mask 0 --emit` →
`runs/replay/{d0,d1}-em-layer/em.mtx` (+6.28% of matrix recovered). External concordance
(`scripts/99_em_af_concordance.py`): EM layer moves median discrepancy vs alevin-fry from 2.5×
to 11%. Downstream-safety probe (`scripts/104_em_de.py`): T-vs-mono DE on EM-heavy genes — 22
of 197 genes move log2FC >0.5; framed as a safety property (D43).

### 8.6 Dataset campaigns (memos `d1-reverification.md`, `d2-glioblastoma.md`, `d2p-brain-nuclei.md`; scripts `90/96/97`)
Each new dataset re-ran: regression byte-identity, fidelity, velocity, EM, size accounting,
premise (annotation-signal) gate. Headline scale law: **bits/read is the stable invariant**
(10.5–18.1 across six datasets; per-molecule cost varies with saturation). D2′ premise: nuclei
churn 4.90% = 2× PBMC; Gene counting uses only 34% of nuclei molecules.

### 8.7 Six-index atlas & junction amortization (D44; `scripts/103_atlas_build.sh`)
`aie federate` over all six archives, junction `chr16:89562391-89562883`:
**1,066,128 UMIs across 93,910 cells**, byte-equal per index, **0.48 s** post-optimization
(3.8 s pre). Amortization: raw junction-set overlap is depth-confounded (18–55%), but
support-weighted (via `index.jpost` totals) **89.2–96.5% of each PBMC sample's junction support
mass** lies on junctions replicated in the other PBMC samples. (Probe script inline in log
D44; per-archive `jpost` decode.)

### 8.8 Reference optimization as replay: Pool et al. (D42; `scripts/102_pool_extension.py`)
Blanket Pool-style GTF (`annotations/gencode.v49.ext3p.gtf`) replayed:
D0 +16.4% UMIs; D2′ **+126%** led by NRXN1 (43k→396k), RBFOX1, NRXN3, KCNIP4, CADM2, GRIK2,
RORA. Outputs `runs/replay/{d0,d2p}-ext3p/`.

### 8.9 Evidence-gated extension (`aie extend`; D45)
The principled version of 8.8: per-gene 3′ extension gated on clustered site support (≥5 UMIs /
≥3 cells), ≤2 kb coverage-gap reachability, corridor-occupancy clipping, and the
internal-priming filter. D0: 5,048 genes, **+2.3%** UMIs; D2′: 22,904 genes, **+20.2%**.
**Decomposition finding:** gated extension shows the NRXN1-class blanket gains involve ~no
movement of true 3′ ends — that mass is intronic pre-mRNA (GeneFull's job), not annotation
error. Outputs: `runs/apa2/gencode.v49.extgated*.gtf`, `runs/apa2/extend-*-report.tsv`,
replays `runs/apa2/replay-extgated*/`.

### 8.10 Genome signature (D45)
Per-contig BLAKE3 of the reference, stamped in archive meta (`aie stamp-genome`, or
`ingest-archive --genome`); sequence-consulting queries verify the supplied FASTA per contig
and refuse mismatches. Stamping is a meta-only rewrite (evidence streams copied compressed);
replay byte-identity through stamping verified on D0. All six archives stamped 2026-08-30.

### 8.11 APA v2: statistics + internal priming + external validation (D45; memo `apa-literature-survey.md`)
- `aie query <a> apa <window> --groups G --genome FA [--permute N]`: site×group multinomial
  G-test (+ label-permutation p), internal-priming flags (≥12 A/20 nt or 8-A run/140 nt
  downstream, transcript orientation). PTPRC: 100/1,290 sites flagged; G-test p=4.4e-129.
- `aie query <a> apa-test --gtf --groups --genome`: genome-wide per-gene G-test + BH-FDR —
  **11,166 genes in 14 s** on D0; 3,355 at q<0.05 (T vs monocyte). Output
  `runs/apa2/apa-test-d0.tsv`.
- **PolyASite 2.0 validation** of the IP filter (site dump `runs/apa2/sites-all.tsv`): at ≥100
  UMIs, 59.9% of passing sites lie within 50 bp of an annotated cluster vs **11.4%** of flagged
  ones (16.9% vs 3.6% at ≥5 UMIs).

### 8.12 Visualization layer (D47; memo `visualization-plan.md`)
`--plot out.{svg,png}` on: `query region` (sashimi: strand-split molecule coverage + junction
arcs + `--gtf` underlay), `query apa` (site lollipops + per-group usage band), `aie em`
(reliability diagram). `query region --export-prefix` emits IGV tracks (per-strand bedGraph +
BED12 junctions). PNG via resvg with an embedded DejaVu Sans (headless hosts have no fonts).
Examples: `runs/apa2/viz/`.

### 8.13 Performance pass (D48; memo `memos/optimization-pass.md` — full baseline/final tables)
Measure-first optimization, byte-identity re-verified after every change. Medians of 3:
replay-rows D0 8.4→**2.6 s**, D1 34.3→**8.1 s** (RSS 35→18.5 GB); `em` D0 285→**12.9 s**
(22×), D1 30.2 min→**64 s** (28×); federate×6 4.24→**0.48 s**; `query region` D1 0.82→0.04 s.
hashbrown 0.17 (foldhash) only where iteration order cannot reach output; EM/ingest maps stay
rustc-hash for determinism; replay aggregation is hash-free (sorted runs + CSR). Baseline
binary preserved at `runs/opt/aie-baseline`; benchmark outputs under `runs/opt/`.

### 8.14 Downstream stability (memo `memos/post-v1-downstream-multiresolution.md`)
The frozen 30-PC, resolution-0.8 Leiden gate is retained as a high-resolution instability
reference. A subsequent locked post-exploratory analysis selected dimension/resolution from the
oracle arm only, requiring six 80% cell-subsample refits, repeated optimizer seeds, adjacent-grid
stability, nontrivial cluster sizes, and split-half marker reproducibility before replay was
opened. Locked settings yield medoid oracle/replay ARI 0.992 (D0), 0.978 (D1), 0.965 (D2′), and
0.997 (mouse 5′); cross-arm loss beyond mean within-arm variability is -0.001--0.024. All four
pass partition, marker, and corruption-control components. D1 and mouse retain their previously
known pseudobulk-threshold misses, so the deliberately broader composite gate remains FAIL.

### 8.15 Federated splice-event biology (D54; scripts `179`--`181`)
The post-gate exploratory screen applies the coordinate-defined event engine across all autosomes
and sex chromosomes. Four marker-scoped PBMC archives yield 20,028 events; 784 have a strictly
consistent nonzero T-versus-monocyte direction with at least ten conservative-denominator
molecules per group in every dataset, and 44 retain minimum absolute usage shift >=0.10. The
dominant FYB1 locus and high-effect CD47 locus match independently RT-PCR-validated immune-cell
programs checked only after ranking. A six-shard called-cell scan yields 46,767 events and 66
strict descriptive PBMC-versus-brain-derived candidates, including previously validated MYL6,
RPS24, and PPP1R12A programs. Exact STAR-SJ aggregation audits the reducer against another
representation of the same alignments. Sequential chromosome totals are 43.78 s / 2.32 GiB and
73.71 s / 2.57 GiB. FNBP1 is held as a candidate under an untouched normal-brain confirmation
locked in `experiments/manifests/federated-biological-screen.yaml`; no confirmatory p-values are
assigned to the exploratory screen. Compact outputs are
`results/post-v1-federated-biological-{screen.json,pbmc.tsv,tissue.tsv}`.

The predeclared FNBP1 test subsequently passed on the untouched 10x 2026 GEM-X Epi Multiome
normal-brain dataset (10,261 called nuclei): archive 72 inclusion-only / 2 exclusion-only
molecules (97.30%), identical targeted annotation-free STAR archive totals, and STAR-SJ 61
inclusion-right / 2 skip UMIs (96.83%). Scripts `182`--`184` download the exact 26.2-GB public
BAM and 158-MB HDF5 by pinned SHA-256, reconstruct the 46,108 targeted primary reads, run STAR,
build both audit archives, and evaluate the unchanged gate. Because reads were selected at the
locked locus, the STAR arm is not a whole-transcriptome mapping-sensitivity comparison. Result:
`results/post-v1-fnbp1-normal-brain-confirmation.json`; script `185` renders the paper figure.

### 8.16 D5 FNBP1 broad-cell localization (scripts `186`--`189`)

Before inspecting event evidence by cell type, a marker-only rule was frozen in
`experiments/manifests/fnbp1-broad-cell-localization.json`. FNBP1 is absent from the marker modules
and excluded from the expression normalization denominator. The rule assigns five broad
compartments and retains explicit unresolved and ambiguous/doublet strata. Because the targeted
archive contains only 6,584 nuclei represented by local reads, its strict query scope is obtained
by restricting the frozen 10,261-nucleus H5 label map without relabeling; rates retain the full H5
group sizes as denominators.

The unchanged executable exactly reproduces the 72 inclusion-only / 2 exclusion-only / 0 both
bulk totals in 0.08 s at 135,112 KiB peak RSS. Event-bearing cells occur in every broad
compartment. Oligodendrocyte/OPC has 30 event-bearing cells (9.36 per 1,000 full-H5 nuclei), while
both exclusion-only molecules occur among four informative microglia/immune molecules (usage
0.50, Wilson 95% interval 0.15--0.85). This is a single-donor descriptive localization with 3,533
unresolved labels, not evidence of differential or cell-type-specific splicing. Result:
`results/post-v1-fnbp1-broad-cell-localization.{json,tsv}`.

The packed federated reducer, exact row-denominator pushdown, bounded-concurrency chromosome
driver, and streaming summarizer are recorded in
`memos/federated-query-workflow-optimization.md` and
`results/post-v1-federated-query-workflow-optimization.json`.  Raw chromosome JSON and timing
directories remain excluded; the compact result records executable/output digests, exact-equivalence
gates, and resource measurements.

### 8.17 Prospective multi-donor FNBP1 validation (GSE234790; scripts `190`--`193`)

The fixed FNBP1 event and unchanged D5 broad-cell marker rule were registered for eight independent
adult human SEZ donors before raw sequence acquisition or event inspection. The lock was pushed as
commit `95b79ad9190063d3a6702f8bdfab7d3cf9cd1cb6` at
`2026-09-01T11:01:55-04:00`. GSE234790 / PRJNA983239 resolves to eight experiments and 40 paired
ENA runs: 977,312,622 read pairs and 44,692,024,247 compressed FASTQ bytes with complete MD5
identities. Donor is the replicate; age is not tested.

The event-blind feasibility run parses metadata and non-FNBP1 expression only. The FNBP1 counts
row is skipped without tokenizing or retaining its values, and no sequence, alignment, junction,
archive, or event input is opened. All 36,626 nuclei and all locked markers pass. Frozen-rule
microglia/immune labels per donor are A 258, B 116, C 99, D 395, E 71, F 115, G 527, and H 42;
seven donors exceed the predeclared feasibility minimum of 50. The feasibility run therefore
passes in 15.76 s at 102,864 KiB peak RSS.

The registered microglia outcome is inferential only with at least six eligible donors. Four or
five donors are a descriptive STOP, because even uniformly negative signs cannot attain the
registered two-sided exact-sign-test threshold (`p=0.125` and `p=0.0625`, respectively). Compact
result and frozen group maps are in `results/post-v1-fnbp1-multidonor-sez-feasibility/`; the full
protocol is `experiments/specs/fnbp1-multidonor-sez-validation.md`.

Script `191` then acquired all 80 ENA FASTQs by their registered size and MD5 identity. The compact
acquisition manifest has SHA-256
`f3822f7c22880fee1f9a73ed414337fa18b5019f08cbe1572e2de22ac56a7126`. Script `192` aligned all
977,312,622 read pairs without a GTF and streamed the full-genome BAM output through the fixed
FNBP1 interval before Gravlax ingest. Donors were isolated and run eight-way concurrently. Script
`193` applied only the registered fixed event, frozen group maps, donor eligibility thresholds,
and STOP/PASS/FAIL rules.

The primary outcome is `STOP_UNDERPOWERED`: only donor G reaches ten informative molecules, versus
the four eligible donors required. The microglia outcome is `STOP_DESCRIPTIVE`: no donor has a
microglial informative molecule, versus at least six donors meeting the compartment-depth minima
required for inference. Descriptively, donors A--H have 4, 0, 1, 9, 6, 4, 13, and 2 informative
molecules, respectively; all 39 are inclusion-only. STAR-SJ independently reduces the same
annotation-free alignments to identical per-donor totals. These descriptive molecule counts do
not override the donor-level STOP.

The first parallel wrapper completed all donor work in 22:24.87 but then exited 1 at collation due
to a GNU-time label parser bug. Commit `6a784cb` repairs that execution-only parser. The preserved
checkpoint resumed successfully in 6.49 seconds, and direct summarization completed in 6.74 seconds;
the timings are recorded separately. Compact acquisition, outcome, donor-group, resource, and memo
artifacts are `results/post-v1-fnbp1-multidonor-sez-*` and
`memos/fnbp1-multidonor-sez-validation.md`; all large sequence and alignment artifacts remain
excluded.

### 8.18 Molecular splice graph v1 (script `194`)

The semantics and acceptance criteria were pushed before implementation as scripts commit
`1c54ec1aa9bc2ad6bf3eabdc722a80897551aacb`. Gravlax commit
`2d4d51e0ae47bbfa402d3d39ba797f588a067e49` implements a strand-aware graph whose edges are
archived junctions and whose path fragments are the selected junction sets co-supported by one
archive UMI class. The command applies the existing strict cell/group/bulk scopes, unions posting
lists, decodes each selected chunk once, deduplicates classes across representatives and chunks,
and derives edge counts from retained path counts. `--max-paths` is a hard early error. Output
schema `gravlax.query.splice-graph.v1` states that fragments are lower bounds rather than complete
transcripts and that multimappers retain anchor-placement semantics.

Script `194` validates the versioned output, both global and per-group edge/path UMI conservation,
strand separation, a real multi-edge fragment, repeated byte identity, and the frozen semantic
guards. The synthetic integration test also checks repeated representatives, opposite strands,
reverse edge direction, identical selected-cell totals at cell/group/bulk aggregation, malformed
inputs, and overflow; all 119 workspace tests pass.

The real capability check uses the donor-G FNBP1 archive and frozen group map from section 8.17.
In `chr9:129907500-129926000`, 21 independent junction postings collapse to one chunk decode. The
result contains 37 nodes, 22 strand-specific edges, and 24 path patterns representing 51
strand-path UMI counts. Three reverse-strand two-edge patterns carry six UMIs. Twenty warm release
invocations take 0.004785 s median (0.004249--0.005645 s) with 13,204 KiB maximum peak RSS. The
compact result is `results/post-v1-molecular-splice-graph-v1.json`; full graph JSON and compiler
artifacts remain excluded. This is a one-donor capability check, not a population contrast. A
future replicate layer must preserve one row per donor or sample before inference.

### 8.19 Replicate-aware splice graph v1 (script `195`)

Registration commit `8d01b5ac87c327598a9ac1e980f07c767d8fe079` fixed the design contract,
common-graph rule, sample-depth eligibility, beta-binomial model, correction family, failure
guards, and counts-only real check before implementation. Gravlax commit
`173d9256b48eb4d729e73732b8a0b7e0cf01c3ff` adds `aie cohort splice-graph`. The strict TSV has
one unique archive per biological sample; a coordinate recurs when local catalogue support is
at least one in at least two archives. Every archive is then reduced against the common graph,
yielding complete zero-explicit sample × path and sample × edge matrices.

For inference, one sample contributes path count $k_i$ and all same-strand common-graph fragment
count $n_i$. Samples below the ten-fragment default remain reported but are excluded. The null
fits one mean and beta-binomial concentration; the alternative fits condition-specific means and
a shared concentration. Deterministic profile maximum likelihood feeds a one-degree-of-freedom
asymptotic likelihood-ratio test, with BH correction across tested locus paths. Integration tests
cover exact count conservation, low-depth suppression, complete zeros, strong/null controls,
duplicate archives, malformed designs, unknown cells, zero thresholds, overflow, and genome
signature mismatch. All 124 workspace tests pass; the Node-22 documentation build also passes.

The real gate uses the tracked design
`experiments/designs/fnbp1-sez-cohort-splice-graph.tsv` and the eight targeted called-cell
GSE234790 archives from section 8.17 under one neutral condition. It is counts-only by
registration. From 44 union coordinates, eight recur in at least two archives and form ten
strand-aware edges, 15 nodes, and 12 paths. The 96 sample-path rows contain 63 explicit zeros;
the 80 sample-edge rows contain 49. Ten one-edge patterns carry 75 fragments and two two-edge
patterns carry 11. Only donors D, E, and G reach ten reverse-strand fragments; no forward-strand
sample does. Twenty byte-identical warm runs take 0.006347 s median
(0.005930--0.007258 s) and 14,844 KiB maximum RSS.

No donor contrast or p-value was constructed. The registered primary `STOP_UNDERPOWERED` and
microglial `STOP_DESCRIPTIVE` outcomes remain binding. The compact result is
`results/post-v1-replicate-aware-splice-graph-v1.json`; full query JSON and compiler artifacts
remain excluded. Promotion to a biological population claim requires a prospectively locked,
genuinely replicated two-condition dataset, comparison with an established differential-splicing
method, and sample-label null calibration. A joint Dirichlet-multinomial or hierarchical
Dirichlet graph model remains a reserved extension, not a current result.

### 8.20 Prospective protocol-native SEZ transcript-end atlas

Before any whole-genome GSE234790 endpoint inspection, script
`196_build_sez_transcript_end_feasibility.py` verified the published metadata digest and froze
the source paper's cluster identities into one barcode/group map per donor. The primary paired
contrast is published Astro/NSC clusters 0/7/16 versus mature neuronal clusters
2/3/10/12/19. All 36,626 published-QC barcodes are retained; clusters 15/17/18 are labeled
`other`, not discarded. Every donor has at least 164 nuclei in each primary group; the event-blind
feasibility result is `PASS_EVENT_BLIND`. Donors A--D are youth and E--H middle-aged, but age is
secondary and descriptive under the registered four-versus-four design.

The prospective gate is `experiments/manifests/sez-transcript-end-atlas-v1.json`, with full
semantics in `experiments/specs/sez-transcript-end-atlas-v1.md`. It fixes common cross-donor site
construction, unique same-strand gene assignment, high-confidence PolyASite/motif/IP tiers,
paired donor inference, a directional neural-differentiation hypothesis, shuffled-label and
site-gap controls, exact-count checks, and <=90-s/<=8-GiB performance targets. Alignment endpoints
must not be called cleavage sites without tail or orthogonal support. The feasibility JSON digest
is `ad9b0c1e...e351b`; the frozen design digest is `867a3fcb...f0ea5`.

### 8.21 Protocol-aware SEZ PolyASite mixture and NTRK2 result

The production implementation is Gravlax commit
`e378516a0fb6602e63d626d0c49fb11e58f463c8`.

The registered raw-endpoint analysis revealed the expected fragmentation geometry of 10x 3′
libraries and failed its global effect-size threshold. The post-hoc method amendment was remotely
frozen in `experiments/manifests/sez-polyasite-mixture-posthoc-v1.json` before inspecting the final
cross-fitted 12/48-bp and shuffled-label controls. It keeps the frozen donor contrast and
inference, but takes site identity from PolyASite, removes genomic internal-priming candidates,
learns a label-blind fragment-distance kernel from unique-candidate molecules, and fits each donor
with a leave-one-donor-out kernel.

All method-amendment controls pass. The primary run has 45,705 recurrent sites, 35.7 million
assigned expected UMIs, 2,182 tested genes, and 95 reported genes. Every held-out donor gains
0.418--0.598 bits/UMI over a uniform distance model; the shuffled arm reports zero genes; and
93/95 primary gene signs persist at both neighboring merge gaps. Simulation quantifies the
close-site/low-depth resolution limit. Sparse likelihood caching reduces the primary run from
37.16 to 22.02 s while preserving the three numerical TSVs byte-for-byte; this is 5.18× faster
than eight independent legacy scans and uses 3,660,468 KiB maximum RSS.

The strongest exploratory gene is NTRK2. A GENCODE-structure-only audit assigns sites near the
proximal shorter-coding terminal group and the distal group containing the MANE Select transcript.
The distal MANE/full-length-like fraction changes from 20.0% in Astro/NSCs to 78.1% in mature
neurons, with the shift positive in all eight donors (median +56.8 percentage points; exact
sign-flip p=0.0078125). This is terminal-group evidence, not full-molecule isoform phasing. The
failed preregistered global-lengthening gate remains binding. Full records and exact reproduction
paths are in `memos/sez-polyasite-mixture.md` and `results/post-v1-sez-polyasite-mixture*`.

### 8.22 Federated sparse cohort output validation

The final sparse-output run used the same six called-cell chromosome-17 cohort-event command twice,
changing only `--json` to
`--sparse-dir runs/post-v1/federated-sparse-output-r2-final/six-called.chr17`. Inputs were archives
D0/D1/D2/D2p/D3/D4, their fixed called-cell maps, compiled GENCODE v49 labels, minimum support 2,
minimum recurrence 3 samples, minimum total informative count 50, and the every-row minimum 20.
The matched dense JSON is
`runs/post-v1/federated-sparse-output-r2-final/six-called.chr17.json`; GNU-time records and the
sparse tables remain excluded large/rebuildable run artifacts. The run used clean Gravlax commit
`6b92503cffbd48da8312bdd75c1f1f70d5ff2a27` and binary SHA-256
`3c5faaf38d8b9ea14001a8561b78b026452804d454b16a5f2e11dbc089605a4d`.

`scripts/205_validate_federated_sparse_output.py` decompresses the three tables and fails closed on
schema/header drift, malformed unsigned counts, duplicate or reordered facts, unexpected explicit
zero facts, inconsistent derived quantities, metadata/dimension disagreement, non-matched commands,
or any difference from the dense event/sample/group records. Missing presence means catalogue
absence; missing count facts mean exact zeros only over the declared event/sample/scope dimensions.
The frozen output uses the final `annotation_genes_json` event column. The validator retains a
compatibility reader for the early split gene-ID/name pilot layout, while requiring either layout
to reconstruct the dense ordered gene objects exactly.

The gate passes: 543 event definitions, 2,383 positive catalogue-presence facts, and 6,516 nonzero
total/group facts reconstruct all 3,258 dense event–sample observations exactly. The compressed
tables are 72,741 bytes; with 4,229 bytes of metadata the bundle is 76,970 bytes, 45.27× smaller
than the 3,484,594-byte JSON. Sparse and JSON arms took 3.92 and 3.96 seconds with 2,473,920 and
2,452,780 KiB maximum RSS, respectively. Cache order was not counterbalanced, so timing is
descriptive rather than evidence of a speed difference.
The compact record is `results/post-v1-federated-sparse-output-validation.json` and interpretation
is in `memos/post-v1-federated-sparse-output.md`. Reproduce the compact record with:

```sh
python3 gravlax-paper-scripts/scripts/205_validate_federated_sparse_output.py \
  --dense-json runs/post-v1/federated-sparse-output-r2-final/six-called.chr17.json \
  --sparse-dir runs/post-v1/federated-sparse-output-r2-final/six-called.chr17 \
  --dense-time runs/post-v1/federated-sparse-output-r2-final/six-called.chr17-json.time.txt \
  --sparse-time runs/post-v1/federated-sparse-output-r2-final/six-called.chr17-sparse.time.txt \
  --metadata-stdout runs/post-v1/federated-sparse-output-r2-final/six-called.chr17.metadata.stdout.json \
  --binary env/cargo-target-atlas-v2/release/aie \
  --gravlax-commit 6b92503cffbd48da8312bdd75c1f1f70d5ff2a27 \
  --out gravlax-paper-scripts/results/post-v1-federated-sparse-output-validation.json
```

### 8.23 Adaptive annotation-independent evidence frontier (rejected)

This negative pilot is the third track of the atlas/evidence-v2 campaign. The umbrella registration
in `experiments/manifests/atlas-collection-and-evidence-v2.json` established that the archive-v1
two-representative default could not change without measurable value. Before implementation or
frontier measurement, the dedicated specification
`experiments/specs/adaptive-evidence-frontier.md` froze a stricter binding gate: no more than 3.0%
archive growth, no more than 1.15× Gene replay wall time, exact weight/determinism/roundtrip checks,
and either at least 0.015 percentage-point D0 Gene improvement or a 0.50-point single-component /
0.25-point mean Velocity improvement. Its pre-implementation content hash was
`3f0a037a...c88505`. After the all-chain arm failed the resource limits, one and only one rescue was
pre-declared in the same file: apply the identical frontier only to junction-containing chains and
retain two extremes for empty-junction chains. The augmented specification hash, recorded before
the rescue archive, is `ba311385...13eff3`.

The experimental binary was built from Gravlax base
`e378516a0fb6602e63d626d0c49fb11e58f463c8` plus a local, uncommitted prototype touching only
`crates/aie/src/{archivecmd,debugcmd,rows}.rs`. The rejected source was not merged or pushed. For
identity, the 27,866-byte binary diff has SHA-256 `1c527b09...f0e6f`. Rust was 1.98.0; the binary's
standing default selected 24 Rayon threads on the 256-logical-CPU snapshot host. The exact
91,679,920-byte executable has SHA-256 `a73a0024...364de4`. Commands below are the executed command
forms; the opt-in flags exist only in that rejected prototype.

```sh
AIE=$P/env/cargo-target-frontier/release/aie
ROOT=$P/runs/post-v1/adaptive-frontier-pilot-r1

# D0 paired archives. With no option, the final pilot binary reproduced the old bytes exactly.
$AIE ingest-archive $P/runs/ingest/full/Aligned.sortedByCoord.out.bam \
  --whitelist $P/data/3M-february-2018.txt --zstd-level 19 \
  --out $ROOT/d0-control/d0.aie
$AIE ingest-archive $P/runs/ingest/full/Aligned.sortedByCoord.out.bam \
  --whitelist $P/data/3M-february-2018.txt --zstd-level 19 \
  --unique-representatives containment-frontier --out $ROOT/d0-frontier/d0.aie
$AIE ingest-archive $P/runs/ingest/full/Aligned.sortedByCoord.out.bam \
  --whitelist $P/data/3M-february-2018.txt --zstd-level 19 \
  --unique-representatives junction-containment-frontier \
  --out $ROOT/d0-junction-frontier/d0.aie

# Gene replay; GNU time -v wrapped each of three sequential runs per paired arm.
$AIE replay-rows $ROOT/d0-frontier/d0.aie \
  --gtf $P/annotations/gencode.v49.annotation.gtf \
  --barcodes $P/runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv \
  --out-dir $ROOT/d0-replay-frontier

# Velocity used the same archive/GTF/barcode universe with --velocity.
$AIE replay-rows $ROOT/d0-junction-frontier/d0.aie \
  --gtf $P/annotations/gencode.v49.annotation.gtf \
  --barcodes $P/runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv \
  --out-dir $ROOT/d0-velocity-junction-frontier --velocity

# Brain-nuclei confirmation: all-chain only; the rescue stopped at the binding D0 accuracy gate.
$AIE ingest-archive $P/runs/ingest/d2p/Aligned.sortedByCoord.out.bam \
  --whitelist $P/data/3M-february-2018.txt --zstd-level 19 \
  --unique-representatives containment-frontier --out $ROOT/d2p-frontier/d2p.aie
$AIE replay-rows $ROOT/d2p-frontier/d2p.aie \
  --gtf $P/annotations/gencode.v49.annotation.gtf \
  --barcodes $P/runs/oracle/d2p/v49/Solo.out/Gene/raw/barcodes.tsv \
  --out-dir $ROOT/d2p-velocity-frontier --velocity

# Called-cell Gene and Velocity fidelity reductions (run once per replay arm).
UV_CACHE_DIR=$P/env/uv-cache uv run --with numpy --with scipy python \
  $P/gravlax-paper-scripts/scripts/30_gateb_compare.py \
  $P/runs/oracle/full/v49/Solo.out/Gene/raw $ROOT/d0-replay-frontier \
  oracle frontier --cells \
  $P/runs/oracle/full/v49/Solo.out/Gene/filtered/barcodes.tsv \
  --out $ROOT/d0-replay-frontier/gene-vs-oracle.json
UV_CACHE_DIR=$P/env/uv-cache uv run --with numpy --with scipy python \
  $P/gravlax-paper-scripts/scripts/91_velo_gate.py \
  $P/runs/oracle/d2p/v49/Solo.out/Velocyto/raw $ROOT/d2p-velocity-frontier \
  $P/runs/oracle/d2p/v49/Solo.out/Gene/filtered/barcodes.tsv
```

The all-chain D0 arm retained 48,640,255 representatives across 21,601,508 chains while preserving
the exact 56,981,029-read chain weight. It improved called-cell Gene movement from 0.243904% to
0.215574% (0.028330 points) and mean Velocity relative L1 by 0.413260 points. It nevertheless grew
from 111,079,672 to 132,487,429 bytes (**+19.272%**) and changed median Gene replay from 3.13 to
3.92 s (**1.252×**): both resource gates fail. The extra payload localized primarily to
empty-junction evidence: 6,387,206 such chains exceeded two representatives, versus 664,661
junction-containing chains.

The pre-declared junction-only rescue is 114,201,153 bytes (**+2.810%**) and replays in 3.25 s
against its paired 3.17-s control (**1.025×**), so both resource gates pass. But Gene movement is
fractionally worse at 0.244068%, and its largest Velocity improvement is only 0.000130 percentage
points. It therefore fails the binding accuracy gate by orders of magnitude. D2-prime independently
confirms the unfavorable all-chain trade: +7.534% bytes, 1.091× Velocity replay, and only 0.234898
point best-component / 0.126011 point mean Velocity improvement. No D2-prime rescue was run after
the D0 STOP.

The formal guarantee is intentionally narrow. For the interval-containment partial order inside a
fixed absolute junction chain, the union of observed maximal and minimal elements is the unique
smallest observed subset preserving every existential predicate that is upward-closed or
downward-closed. It does not preserve interior endpoint-band predicates, geometry multiplicity,
arbitrary Velocity state changes at exon boundaries, sequence/edit/quality predicates, or complete
isoform phasing; executable tests include the first two counterexamples.

All 143 workspace tests pass, including wide-token roundtrip and malformed-stream rejection. The
final pilot binary's default D0 archive is byte-identical to the paired control; independent rescue
ingests and rescue eager/streaming Gene replays are byte-identical; all exact chain weights are
conserved. The compact record is `results/post-v1-adaptive-evidence-frontier.json`, interpretation
is in `memos/adaptive-evidence-frontier.md`, and every excluded archive, matrix, and timing-log
identity is in `manifests/adaptive-evidence-frontier-artifacts.tsv`. The final verdict is **STOP**:
do not merge either prototype, change the default archive, weaken the gates, or describe this as a
production archive-v2 format.

### 8.24 Bounded-memory PolyASite cohort execution

The fourth campaign track replaced atlas-wide endpoint and nested count structures in the
eight-donor PolyASite mixture analysis with a bounded execution frontier at Gravlax commit
`655a7de923e71d9eaf0ef5b40f2a0afb3c3f3e55`. The executable was built from a clean detached
checkout with Rust 1.98.0; its SHA-256 is
`b7af8ad99d78a260cf74df5c73b211238587a3c9f3a9732dc0a131453302ddb9`.

The exact executed analysis command is the primary gap-24 invocation recorded by
`scripts/203_benchmark_polyasite_mixture_hotpath.sh`, using the eight SEZ called-cell archives,
the frozen donor/group design, PolyASite catalogue, reference, and leave-one-donor-out kernels in
`runs/post-v1/sez-polyasite-mixture-crossfit-controls-r1/primary-gap24`. GNU time recorded the
bounded arm. The fail-closed reduction was:

```sh
python3 gravlax-paper-scripts/scripts/204_summarize_bounded_polyasite.py \
  --frozen runs/post-v1/sez-polyasite-mixture-crossfit-controls-r1/primary-gap24 \
  --bounded runs/post-v1/sez-polyasite-mixture-bounded-r1/primary-gap24 \
  --time runs/post-v1/sez-polyasite-mixture-bounded-r1/primary-gap24.time.txt \
  --out gravlax-paper-scripts/results/post-v1-polyasite-bounded-memory.json
```

The summarizer requires exit status zero, hashes and compares `sites.tsv`, `genes.tsv`, and
`fragment-kernel.tsv` byte for byte, checks a 30.0-s wall limit and a 2.25-GiB (2,359,296-KiB)
maximum-RSS limit, and refuses incomplete timing records. All gates pass. The bounded arm takes
28.21 s at 2,262,980 KiB versus the prior optimized arm's 22.02 s at 3,660,468 KiB: a 38.18% RSS
reduction for a 1.281× wall-time ratio. All three table digests are identical to the frozen
result. The compact record is `results/post-v1-polyasite-bounded-memory.json`; the interpretation
and the deliberately limited memory claim are in `memos/post-v1-polyasite-bounded-memory.md`.

### 8.25 Optional federated collection index

The fifth atlas/evidence-v2 track built an optional derived collection index over the eight called-
cell SEZ archives from §8.20. The locked run is
`runs/post-v1/collection-index-v2-r5`, generated by clean Gravlax commit
`6b92503cffbd48da8312bdd75c1f1f70d5ff2a27` and release binary SHA-256
`3c5faaf38d8b9ea14001a8561b78b026452804d454b16a5f2e11dbc089605a4d`. R1 and r2 were failed
benchmark-driver attempts; r3 and r4 were successful pre-final dry runs. They are execution
history, not scientific arms. R5 is the sole final benchmark and adds a matched timed aggregate of
the eight ordinary archive queries for each query kind.

The exact driver and reduction are:

```sh
cd "$P"
AIE="$P/env/cargo-target-atlas-v2/release/aie" \
BENCH_OUT="$P/runs/post-v1/collection-index-v2-r5" \
  gravlax-paper-scripts/scripts/206_benchmark_collection_index.sh

python3 gravlax-paper-scripts/scripts/207_validate_collection_index.py \
  --run-dir runs/post-v1/collection-index-v2-r5 \
  --binary env/cargo-target-atlas-v2/release/aie \
  --gravlax-commit 6b92503cffbd48da8312bdd75c1f1f70d5ff2a27 \
  --out gravlax-paper-scripts/results/post-v1-collection-index-v2.json
```

The driver supplies and the builder verifies the recorded BLAKE3 digest for each source archive.
The panel totals 1,326,036,691 bytes and 6,848 chunks and shares stamped reference digest
`2817ea0b...47472`. It builds a fresh eight-archive root in A--H and H--A input order, then builds
an A--D base and an E--H extension. It records the base digest before and after extension. The
fresh roots are byte-identical (SHA-256 `9b33caf1...05245`), and the base is unchanged (SHA-256
`25e322d6...163bb`). The 8,931,229-byte root is 0.674% of the archive panel; the 4,751,118-byte
base plus 4,718,653-byte extension is 0.714%. Both clear the 5% gate.

Five fixed root queries and the corresponding two-layer queries are logically identical. Each of
their eight sample rows is also compared to an ordinary source-archive query, for 40 exact
comparisons including explicit absence markers. The fixed totals are 182 UMIs / 180 cells at the
dense FNBP1 junction; 4 UMIs / 4 cells at the sparse neighbor (A and G only); zero at the absent
neighbor; 74,848 molecules / 74,352 UMIs / 20,019 sample-cells across the NTRK2 region; and 182
include-only / 4 exclude-only / 0 both in the joint set. A support threshold of 211 exceeds the
dense catalogue upper bound of 210 and prunes every source without opening an archive.

Root wall times for dense, sparse, absent, region, and joint-set queries are 0.13, 0.12, <0.01,
0.14, and 0.13 seconds. The matched sequential eight-query aggregates are 0.45, 0.24, 0.11, 0.18,
and 0.52 seconds. These are single warm-host observations, not replicated distributions. Peak RSS
is higher than the sequential naive aggregate for dense and joint-set output, so no general memory
advantage is claimed. Source-reader instrumentation records 36,588,673, 14,100,036, zero,
60,464,172, and 36,588,673 archive bytes, respectively. The joint plan decodes eight unique chunks
instead of ten independent component decodes.

`--verify-content` on the dense query accounts for exactly 1,326,036,691 identity bytes—one full
pass over the source panel—and leaves all scientific rows unchanged. Ordinary source guards use
the opened file's size, nanosecond mtime/ctime, device, and inode. The sidecar directory and every
decompressed section are authenticated and full inspection performs a semantic audit. Because v1
archives do not yet expose a reusable archive-root commitment, index construction still reads all
compressed source bytes to verify identity; it does not decompress molecule chunks. The 0.76-s
root build and 0.20-s verified query reflect warm local storage and must not be extrapolated to cold
or remote filesystems.

`scripts/207_validate_collection_index.py` fixes and hashes the complete 114-file r5 layout, checks
all four sidecar and executable SHA-256 records, parses 24 successful GNU-time records, validates
root/chain structure and internal byte/chunk/phase accounting, and emits a canonical digest over
every excluded run file. The compact record is `results/post-v1-collection-index-v2.json`; complete
interpretation and limitations are in `memos/post-v1-collection-index-v2.md`.

### 8.26 Sparse terminal-tail evidence and archive-promotion hard stop

This campaign deliberately separates a successful evidence pilot from an unsuccessful archive
representation. The untouched Gravlax base is
`e378516a0fb6602e63d626d0c49fb11e58f463c8`. The biological gates were frozen in
`ff0cecc75c31aaa9d1c67f47de243cca0ab15500` before cohort enrichment was examined, and the
experimental sidecar/result landed in
`db817cf967c4c6d8df442bb1a73efa81e1152909`. The one-u16-per-ordinal production gate was then
frozen in `50b50a2360d11caf85eb3b1c56c16033f6863546` before its structural audit; the executable
counterexample and hard-stop result are
`d8b7482e1995433d9b017a098ab65b1c700cc2f2`.

The original pilot and promotion specifications have SHA-256 identities
`12de8559...af1ecf` and `25f02890...8d5e1`; the aggregation script is copied byte-for-byte as
`scripts/208_aggregate_sparse_terminal_tail.py` with SHA-256 `c802e5ef...b6d52`. The isolated
pilot code patch from base through `db817cf` hashes to `fb17ba45...abf19`; the later counterexample
patch to `rows.rs` hashes to `31a14af3...00f1a`. Those hashes identify rejected experimental
source without copying it into this repository. The final 95,068,232-byte pilot executable has
SHA-256 `7733f923...d268c`.

The measurements used Rust 1.98.0 (`rustc 88d9e12ae`, Cargo 1.98.0) on Linux
4.18.0-553.156.1.el8_10, a two-socket AMD EPYC 9575F host with 256 logical CPUs. The standalone
tail scanner itself is single-threaded. The code was built and checked in its isolated worktree:

```sh
P=${GRAVLAX_PROJECT_ROOT}
export CARGO_HOME=$P/env/cargo
export RUSTUP_HOME=$P/env/rustup
export CARGO_TARGET_DIR=$P/env/cargo-target-tail-pilot
cd $P/work/gravlax-tail-pilot

# At db817cf: 141 tests, including orientation, packing, roundtrip, AUC, and exact-key dedup.
cargo +1.98.0 test --workspace
cargo +1.98.0 build --release
AIE=$CARGO_TARGET_DIR/release/aie
```

The exact eight-donor execution used the sequence-bearing annotation-free STAR BAM, called-cell
whitelist, matching archive, PolyASite catalogue, and uncompressed GRCh38 reference recorded in
`manifests/sparse-terminal-tail-artifacts.tsv`. Defaults are shown explicitly here; each output
path was absent because both commands refuse overwrite:

```sh
ROOT=$P/runs/post-v1/sparse-terminal-tail-pilot-r2-dedup
for D in A B C D E F G H; do
  mkdir -p "$ROOT/donor-$D"
  "$AIE" tail-evidence-pilot \
    "$P/runs/post-v1/sez-transcript-end-atlas-r1/donor-$D/Aligned.out.bam" \
    --whitelist "$P/runs/post-v1/sez-transcript-end-atlas-r1/donor-$D/called-whitelist.txt" \
    --polyasite "$P/runs/apa2/polyasite.bed.gz" \
    --genome "$P/ref/GRCh38.primary_assembly.genome.fa" \
    --archive "$P/runs/post-v1/sez-transcript-end-atlas-r1/donor-$D/sez.called.aie" \
    --min-clip 6 --min-tail-fraction 0.8 --min-terminal-run 4 \
    --site-window 25 --control-shift 500 --control-exclusion 100 --zstd-level 9 \
    --out "$ROOT/donor-$D/tail-evidence.ait" \
    --json-out "$ROOT/donor-$D/audit.json"
done

python3 "$P/gravlax-paper-scripts/scripts/208_aggregate_sparse_terminal_tail.py" \
  --run-root "$ROOT" \
  --workspace-tests-passed --workspace-test-count 141 \
  --aie-binary "$AIE" \
  --out "$ROOT/cohort-summary-locked.json"
```

The aggregate has 715,967,790 mapped primary reads, 136,763 exact-key tail witnesses, and 53,346
external-site witnesses. External candidates have 53,345 positives among 186,470 terminal-clip
witnesses, versus 763 among 20,386 shifted controls (odds ratio 10.299). The continuous tail score
improves AUC from 0.526843 to 0.654823 (+0.127980). Tail positivity is 9.00% at genome-flagged
internal-priming sites versus 31.39% at non-flagged sites. The eight sidecars occupy 1,671,772
bytes against 1,326,036,691 archive bytes (0.126073%). Every donor independently clears the
enrichment, AUC, storage, and internal-priming gates; the weakest odds ratio and AUC gain are 6.962
and 0.087829.

The release binary was independently rerun on donor A under
`runs/post-v1/sparse-terminal-tail-pilot-r3-release-verify`. Its 73,709-byte sidecar is byte-
identical to the measured r2 sidecar (`dabf7114...b90014`), and the audit is identical after
removing only elapsed time and output path. All relevant inputs and both outputs are hashed in the
artifact manifest.

The secondary NTRK2 query used the nine-site TSV frozen with the pilot and the already frozen
Astro/NSC and mature-neuron cell maps:

```sh
SAMPLE_ARGS=()
for D in A B C D E F G H; do
  SAMPLE_ARGS+=(
    --sample "$D=$ROOT/donor-$D/tail-evidence.ait,$P/gravlax-paper-scripts/results/post-v1-sez-transcript-end-feasibility/groups/donor-$D.tsv"
  )
done
"$AIE" tail-evidence-query "${SAMPLE_ARGS[@]}" \
  --sites "$P/work/gravlax-tail-pilot/pilots/sparse-terminal-tail/ntrk2-sites.tsv" \
  --group-contrast astro_nsc:mature_neuron \
  --numerator-set distal_mane_full_length_like \
  --denominator-set proximal_shorter_coding \
  --window 25 --json-out "$ROOT/ntrk2-query.json"
```

This secondary check fails: only three witnesses occur, all in donor G and all at the proximal
shorter-coding group (one Astro/NSC, two mature-neuron), versus frozen minima of 20 witnesses and
six donors. It is not evidence against the fragment-mixture NTRK2 result, but it cannot confirm
that result and cannot phase coding-to-terminal isoforms.

Finally, the production-layout test was run after the prospective hard-stop gate:

```sh
cd $P/work/gravlax-tail-pilot
cargo +1.98.0 test --workspace
cargo +1.98.0 test -p aie \
  rows::terminal_tail_attachment_tests::middle_tail_read_endpoint_is_not_recoverable_from_emitted_molecule_on_either_strand
```

All 142 tests pass, including the forward/reverse BAM counterexample. On `+`, three same-chain
spans `90..190`, `100..140`, and `110..130` emit one weight-3 molecule whose retained endpoints
are 190 and 130, losing the middle tail-positive anchor 140. On `-`, input starts 290, 300, and 310
retain 290 and 310, losing the positive anchor 300. Thus one ordinal plus one signal cannot encode
either a missing coordinate or multiple distinct anchors summarized by one molecule.

The structural gate superseded size/runtime measurement for that lossy layout. No production
section was written and no experimental Gravlax code was merged. The minimum exact redesign is a
sparse per-selected-molecule event list: delta-coded local ordinals; an event count per ordinal;
and, for each event, a signed varint cleavage-anchor delta from `MolRec::anchor()` followed by the
frozen u16 signal. Global deduplication is by cell, UMI class, chromosome, strand, and exact anchor.
An optional event-coordinate extent/posting index can preserve targeted query skipping.

The compact three-way verdict is `results/post-v1-sparse-terminal-tail.json`; its interpretation
and claim limits are `memos/sparse-terminal-tail-evidence.md`. The large source BAMs, archives,
sidecars, audits, group maps, reference, catalogue, binary, and release check are excluded from Git
and pinned exactly in `manifests/sparse-terminal-tail-artifacts.tsv`.

### 8.27 Content-addressed archive root and local-shape routing

The two ordered gates are frozen in
`experiments/specs/2026-09-01-archive-root-and-jshape-routing.md` and
`experiments/manifests/archive-root-jshape-routing.yaml`. The exact Gate-A source panel, including
byte sizes and full-file digests, is `manifests/archive-root-gate-a-inputs.tsv`.

The complete execution interface and commands are in
`experiments/specs/2026-09-01-archive-root-gate-a-runbook.md`. In brief, script `209` performs
byte-preserving v1-to-v2 sealing, independent root and payload-digest checks, repeated migration,
reversed collection builds, four-arm scientific query comparisons, Gene/Velocity replay, eight
alternating warm blocks, and the adversarial hook. Script `210` freezes the complete artifact set
and is the only component allowed to emit a compact PASS. It enforces the structured collection
build I/O record, including zero rooted identity-content bytes and the exact six metadata/index
sections.

The confirmatory r3 run passes every gate. It is bound to Gravlax commit
`fbd48e536b39729632f189a74e7ece891166ecd1`, reproducibility commit
`df1e363d397c740300c97ec2b1f04fe48cd9eb6b`, and release-binary SHA-256
`f0fb4c0d6a12bc9b3877bbab09168e5754939fe97ad4a363db51fa5c91c06211`. The complete r3 manifest
binds 813 files and 6,896,522,383 bytes with SHA-256
`d64a953320a56d53e0caaa11d81c91bffbb55a2067027466201f161a49def64d`. The compact result is
`results/post-v1-archive-root-gate-a.json`, committed by
`9cba7af125dfce2d0b0c0b14ea56e92e00215a5d` with SHA-256
`3f75fb1f7d5cc8b726405017a320c3698b7384c5da1f833f6d1f0a29453429db`; the corresponding 59-row
excluded-artifact inventory is `manifests/archive-root-gate-a-artifacts.tsv` (SHA-256
`60ac744b559e42918ebab40570bde981600cfd01736e2f744203958f6e87dada`). A second reduction into an
ignored under-project destination was semantically identical after removing only the requested
inventory-output path, and its inventory was byte-identical.

Across archives A--H, sealing adds 355,840 bytes to 1,326,036,691 bytes (0.026835%), exactly 32
bytes per section plus 32 bytes per archive. The rooted collection build reads zero
identity-content bytes and 10,341,234 total source bytes, versus 1,326,036,691 identity-content
bytes and 1,336,022,077 total source bytes for the v1 control: a 129.19-fold total-I/O reduction.
All Gene/Velocity replay outputs and fixed dense, sparse, absent, region, and junction-set answers
are exact. Median v2/v1 ratios are 0.9221 for build wall, 0.9575 for build RSS, 0.9954 for replay
wall, and 1.0129 for replay RSS; dense and junction-set query wall ratios are 1.0769 and sparse is
1.0000, all within the frozen limits. These warm-host measurements support an identity/build-I/O
claim, not a query-speed claim. The root commits to encoded content; it is not a publisher
signature, and ordinary lazy reads verify selected payloads while `--verify-content` verifies all
payloads.

Two earlier runs are retained only as non-confirmatory harness provenance. r1 stopped before its
first collection build because a planned output parent was absent; it has no completion record or
driver artifact manifest. r2 completed the workload but was used diagnostically to harden manifest
path ordering, replay-diagnostic parsing, repeated replay resource parsing, and comparisons of
identity-dependent metadata size. Its legacy manifest binds 813 files and 6,896,522,424 bytes
(`e6d5758a20022c5eb77f6bb13c278fd8c711d84f29d89e5457e6eb5db7cba6f7`). Neither run contributes
reported treatment estimates. Their excluded identities are recorded in
`manifests/excluded-artifacts.tsv`.

After that Gate-A result commit, script `211` read only `chroms` and `index.junctions` from the
eight rooted archives and materialized the Gate-B panel before any route implementation or
measurement. The frozen panel has 96 distinct junctions, exactly 24 from each archive-presence
stratum `{1, 2--3, 4--7, 8}`, and SHA-256
`55980706cd6fbf0c07a1e652b330e5aca6946b1365010c8ef4e8cb6d6c4f1ca1`. Its exact TSV, digest,
metadata, and the frozen exact-span pair-serialization addendum are under `experiments/designs/`
and `experiments/specs/2026-09-01-jshape-routing-serialization-addendum.md`. The eligible stratum
populations are 1,139,256; 114,543; 68,436; and 22,035 junctions, respectively.

At the prospective freeze boundary, each exact-span route was specified to store sorted unique
`(local_shape_id, donor_offset)` pairs and recover `acceptor_offset` by checked addition of the
bucket span. The frozen addendum proves this representation bijective with the original triples
and retains every gate threshold. The 33,897,246 compressed shape bytes and fixed-query source-byte
audits were pre-implementation inputs; the predicted route size and proposed savings remained
projections until the confirmatory run below.

After the route CLI and structured JSON were finalized, scripts `212`--`214` froze the independent
Gate-B campaign without running it. The fixed adapter requires build schema
`gravlax.collection.build.v2`, inspection schema `gravlax.collection.v4`, and point, junction-set,
and region schemas `gravlax.collection.junction.v2`, `gravlax.collection.jset.v2`, and
`gravlax.collection.region.v2`. It compares route-free and `--shape-routes` collections built by
the same binary in fresh, repeated, reversed-input, and A--D/E--H extension forms; executes eight
alternating warm blocks over the fixed controls and 96-row panel; and records source, sidecar,
total-logical I/O, phase timing, wall time, and RSS separately. Positive geometry, representative,
multimapper, and cross-chunk UMI behavior is bound to ten exact named Rust tests. Ten independent
collection mutations exercise codec, source-root, shape-digest, shape-id, reconstruction,
ordinal, span-bucket, truncation, checksum, and trailing-byte failures through the release CLI.

Script `213` is the only component allowed to emit a Gate-B PASS. It independently checks the
complete artifact manifest, source and route bindings, exact normalized scientific JSON, all byte
arithmetic, all 26 fixture predicates, and every unchanged threshold before atomically exposing a
compact result and excluded-artifact inventory. The exact prerequisites, schemas, commands, and
claim boundary are in `experiments/specs/2026-09-02-jshape-routing-gate-b-runbook.md`.

The sole confirmatory campaign is r3. From the project root, its driver and reduction commands
were:

```bash
source gravlax-paper-scripts/environment/setup.sh
unset AIE_CHROM
cd work/gravlax
CARGO_TARGET_DIR="$PWD/../../env/target-jshape-r3" \
  TMPDIR="$PWD/../../env/tmp-root-main" \
  cargo +1.98.0 build --release -p aie
cd ../..

python3 gravlax-paper-scripts/scripts/212_benchmark_jshape_routing.py \
  --project-root "$PWD" \
  --run-dir "$PWD/runs/post-v1/jshape-routing-gate-b-r3" \
  --aie "$PWD/env/target-jshape-r3/release/aie" \
  --interface-status final \
  --fixture-hook "$PWD/gravlax-paper-scripts/scripts/214_jshape_routing_adversarial.py"

python3 gravlax-paper-scripts/scripts/213_validate_jshape_routing.py \
  --run-dir "$PWD/runs/post-v1/jshape-routing-gate-b-r3" \
  --gravlax-commit 7d997c6066b5a3d1efde9c466d5dce66d9fb8743 \
  --paper-scripts-commit 0de74cd901c39ce50898f3a0a23af05ecd23ef4b \
  --out "$PWD/gravlax-paper-scripts/results/post-v1-jshape-routing-gate-b.json" \
  --artifact-inventory-out \
    "$PWD/gravlax-paper-scripts/manifests/jshape-routing-gate-b-artifacts.tsv"
```

The protocol records both repositories clean at those commits. The single 101,954,648-byte release
binary used for every arm has SHA-256
`356d7699c163d31802e89d3ee54d7b6fccd176701a56bb99dbfd7f747088061a`. The protocol and completion
records have SHA-256 `53a9e6fca9914cdeba460a97eece41e8ca0090d3a007272daa09df0c3307af7c` and
`6d5e924b8b9cc44646f8a3b8c9f9a5dac113f3e8ed4ac748eb5fa49be429fb14`. The completion record is
FINAL, COMPLETE, promotable, and contains all 26 required fixture cases. The complete r3 artifact
manifest binds 8,223 files and 647,307,809 bytes with SHA-256
`70b48020f62a99437027332cd24ceb493c9bdd8c4debbbc8819b4693cefef52e`.

The reducer emits **PASS**. The compact result is
`results/post-v1-jshape-routing-gate-b.json`, 268,715 bytes with SHA-256
`e3c7430ebed30d67711a22e70e6bc05a6d22726218cc4da927f83b646d372848`; the corresponding 37-row,
11,456-byte excluded-artifact inventory is `manifests/jshape-routing-gate-b-artifacts.tsv`, SHA-256
`b4e8d11ccbd64442bc2b6f921195998b7db66eb0fec41d64c61c6ed5cd4ffcf4`. A second independent
reduction to temporary under-project outputs returned PASS; those temporary outputs were removed.

All 96 prospectively selected junctions and the five fixed dense, sparse, absent, junction-set, and
region controls agree exactly between routed and fallback execution, fresh root and immutable
chain: 101 exact scientific cases. The fixed answers remain 182 dense UMIs in 180 cells; 4 sparse
UMIs in 4 cells; zero absent UMIs; 182 include-only, 4 exclude-only, and zero both for the
junction-set control; and 74,848 molecules, 74,352 UMIs, and 20,019 sample-cells for the region
control. Full inspection reconstructs 1,335,446 exact-span rows and 8,092,537 pairs from the source
`shapes` dictionaries. Repeated and reversed-input builds are byte-identical, and all 11 accepted
and 15 rejected integrity predicates behave as frozen.

The route-bearing fresh root is 39,286,118 bytes versus 8,931,406 bytes route-free. Its 30,354,712
route premium is 2.289131% of the 1,326,036,691 source bytes; the complete root is 2.962672% of
source. The immutable route chain is 39,823,041 bytes versus 9,470,591 bytes route-free, a
30,352,450-byte or 2.288960% premium; the complete chain is 3.003163% of source. Compressed route
blocks occupy 30,002,228 bytes, 0.885093 times the prospectively frozen 33,897,246 compressed
`shapes` bytes.

Across eight alternating 24-thread warm pairs, route construction takes 1.340 s median versus
0.705 s fallback (1.900709 times) and 607,132 versus 478,126 KiB peak RSS (1.269816 times), passing
the 2.0 and 786,432-KiB limits. Dense, sparse, and junction-set wall ratios are 0.307692, 0.291667,
and 0.307692; their source-byte ratios are 0.105145, 0.129262, and 0.105145 and total-logical-byte
ratios are 0.135978, 0.198297, and 0.142038. Corresponding routed peak RSS values are 76,262,
39,392, and 77,010 KiB. The route-neutral region control has source ratio 1.0,
total-logical-byte ratio 1.010634, wall ratio 1.038462, and RSS ratio 0.996925. Across the 96-row
panel, source-byte ratios have median 0.119379 and p95 0.332421; total-logical-byte ratios have
median 0.174608 and p95 0.528042; aggregate wall time is 0.386838 times fallback. Every frozen
size, construction, I/O, runtime, and memory gate passes.

Gate-B r1 and r2 remain non-confirmatory. r1 is a preflight-only empty skeleton: it has no files,
protocol, command, completion record, or artifact manifest. r2 completed the workload but was used
diagnostically before implementation optimization. It failed the frozen construction wall
(2.528169 versus 2.0), construction RSS (871,494 versus 786,432 KiB), and region wall (1.153846
versus 1.05) gates. Its artifact manifest binds 8,223 files and 647,324,075 bytes with SHA-256
`ce2076e52dacfd085d4fd639cdc83c3b2e877934ff61a475e48eeb81c69f42de`. r2 motivated treatment
optimization and supplies no reported estimate; r3 alone is confirmatory. The three run-tree
identities are recorded separately in `manifests/excluded-artifacts.tsv`.

The PASS licenses exact point-junction and junction-set acceleration over rooted archives. Routes
supply candidates, never molecular counts; source molecule chunks remain authoritative. The result
does not license replay or region acceleration, a cold-filesystem performance claim, a biological
validity claim, or equivalence to fresh annotation-aware alignment.

### 8.28 Preprint follow-up experiments, 2026-09-11 (scripts `215`--`224`; memo `memos/paper-followup-20260911.md`)

Run with the released Gravlax binaries v0.2.2 (`runs/gravlax-release-0.2.2-20260910/aie-0.2.2`)
and, for the discovery rerun, v0.2.3 (`runs/gravlax-release-0.2.3-20260912/aie-0.2.3`). Records
under `runs/paper-followup-20260911/`; compact copies in `results/paper-followup-20260911/`.

- **D2′ GeneFull timing** (`215`): the script-`120` crossover protocol (one advisory-cold arm per
  method, five warm alternating blocks, 24 CPUs) with STARsolo restricted to GeneFull counting
  and no alignment output versus `replay-rows --gene-full`. Warm medians 361.46 s / 32.55 GiB
  (STARsolo) against 9.27 s / 5.79 GiB (replay); all twelve runs exact against their references.
- **D1 geometry fidelity** (`216`): `ingest-archive --geometry-fidelity --zstd-level 19` on the D1
  ingest BAM. Archive 692,731,515 B (14.44 bits/read) versus 544,022,039 B (11.34) compact.
  Gene moved UMI mass 0.253% versus 0.307%; velocity relative L1 spliced 2.354% versus 2.672%,
  unspliced 0.852% versus 0.859%, ambiguous 4.669% versus 6.132%. The compact replays reproduce
  `runs/replay/d1-aie` and `d1-velocity` byte for byte.
- **Annotation stability** (`217`): GENCODE v32 to v49, transcript Jaccard 0.44 by identifier and
  0.42 by intron chain, junction Jaccard 0.61, 379,610/382,880 v32 junctions retained. Observed
  unique junction reads (at least three) on v49-only junctions: D0 0.24%, D1 0.32%, D2′ 0.73%.
- **Junction-seeded ingest, D0** (`218`, `219`): STAR `--sjdbFileChrStartEnd --sjdbOverhang 90`
  on the base index with v32 or v49 junctions, one- and two-pass. Gene v49 deviation 0.244%
  (unseeded two-pass), 0.228% / 0.213% (v32 / v49 two-pass), 0.160% / 0.143% (v32 / v49
  one-pass); velocity ambiguous 3.059% down to 1.506%. Archive bytes within 1.3% of unseeded.
- **Cohort-wide event discovery** (`220`--`224`): `collection find-events --kinds cassette` over
  the eight rooted SEZ archives (1,326,392,531 B) through the route-root and coordinate-only
  collections, novel-versus-v32 or unfiltered. Under v0.2.2, 83.6--88.0 s and 12.3--12.4 GB
  peak RSS; under v0.2.3 (`223`), 6.4--6.5 s total (8.8--9.0 s wall) and 2.96--2.98 GB, with
  output identical after stripping timing fields. 5,413,686 candidate definitions, 2,976,274
  candidate entities, 205,735 retained (novel-v32), 3,442 of them in the missing-junction class,
  including the FNBP1 cassette (44 UMI classes, 44 cells, six donors). The naive baseline
  (`query events --kind cassette` per archive per chromosome, `222`, `224`) takes 23.9 s plus an
  8.9 s merge, retains 1,196 events, and misses FNBP1 because only two archives generate the
  candidate from their own catalogues.

## 9. Current headline results (as claimed in the paper, all reproduced above)

Size 10.5–18.1 bits/read, 9.0–12.7× vs CRAM 3.1 · fidelity 0.217–0.454% (≥5× under signal) ·
replay 2.6 s / 8.1 s vs 16 / 20.4 min fresh STARsolo (equal 24-thread budgets) · open ~10 ms ·
indexed queries 0.05–0.3 s · six-index federation 0.48 s · velocity ≤6.1% (ambiguous,
decomposed) · discovery 100% clean-stratum prospective recall · APA 26/27 replication +
calibrated G-test + IP filter validated on PolyASite · EM pooled 90.9–98.3% masked recovery,
STARsolo-EM replicated at 0.58% · gated extension +20.2% (nuclei) · 43.78-s four-PBMC
whole-genome event federation recovers FYB1/CD47 external positive controls · prospectively
locked FNBP1 candidate replicates at 97.3% inclusion in untouched normal brain · encoded-section
archive roots add 0.026835% and reduce eight-archive collection-build source I/O 129.19-fold while
preserving exact replay and fixed-query answers · local-shape routes add 2.289131% of source bytes,
retain exact root/chain answers across 101 cases, and reduce prospective-panel aggregate junction
wall time to 0.386838 times fallback.

## 10. Regeneration recipes (minimal, D0)

```sh
cd $P && source env/setup.sh && cd src && cargo build --release && AIE=$PWD/target/release/aie
# 1. oracle + ingest alignments (hours):    scripts/10_oracle_starsolo.sh, scripts/40_ingest_align.sh
# 2. build the index once:
$AIE ingest-archive runs/ingest/full/Aligned.sortedByCoord.out.bam \
     --whitelist data/3M-february-2018.txt --genome data/GRCh38.primary_assembly.genome.fa.gz \
     --out /tmp/d0.aie
# 3. the regression that defines correctness (must be byte-identical to runs/replay/d0-v15):
$AIE replay-rows /tmp/d0.aie --gtf annotations/gencode.v49.annotation.gtf \
     --barcodes runs/oracle/full/v49/Solo.out/Gene/raw/barcodes.tsv --out-dir /tmp/replay
cmp /tmp/replay/matrix.mtx runs/replay/d0-v15/matrix.mtx
# 4. capabilities: replay-rows --velocity | query region/junction/apa/apa-test/discover |
#    extend | em [--star|--emit|--plot] | federate | stamp-genome   (all have --help)
```
The full end-to-end battery for a new dataset is `scripts/96_d2_campaign.sh` (copy-adapt);
`scripts/101_replication_pack.sh` re-runs the capability evaluations.

## 11. Known caveats and incidents (full list in the decision log)

- D16: STARsolo counts unique-gene multimappers via secondaries — tag-derived baselines
  undercount by 2.2% (fixed everywhere; recorded).
- D20: UMI classes must be global (cell,value), not per-chromosome (2.9% over-count, fixed).
- D25: layout meta-guards exist because a stream-count mismatch decodes to garbage silently.
- G-C is MARGINAL by its frozen bar (98.7–98.9% vs 99.0) — treated as a measured ceiling and
  reported as such, not hidden.
- Dev-scale (chr19) numbers are for iteration only; every frontier decision was re-measured at
  full scale after three dev-scale extrapolations failed.
- gates.yaml §datasets/annotations predates D1–D4 and the withheld panel; this file supersedes
  it where they disagree.
- The code repo's history was rewritten once (2026-08-30) to purge accidentally committed
  `node_modules` blobs; current main is `d9257ca` with a 630 KiB pack.
