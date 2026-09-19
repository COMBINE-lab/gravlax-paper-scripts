| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |# Annotation-Independent Molecular Evidence for scRNA-seq — Implementation Plan
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**Host:** analysis-host · **Plan date:** 2026-08-28
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**Source doc:** the original project-plan document (archived as `docs/project_plan.md`)
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |---
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |## 1. Context
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |The previous candidate project (certified L∞ likelihood coresets) passed Gate -1 bit-exactly and
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |then died at Gate 0 on a compression ceiling that was a property of the data, not the algorithm.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |The lesson driving this plan's structure: **the cheapest lethal measurement must run first**, and
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |the previous project's fatal number arrived on day ~12 instead of day ~3.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |This project asks whether a compact, molecule-resolved, annotation-free archive built once from
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |raw 10x reads can be re-quantified under future annotations without returning to FASTQ.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Two findings from planning reshape the source doc's gate order, and both are load-bearing:
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**Finding A — gene assignment is provably replayable from coordinates; UMI collapsing provably is not.**
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |I read STAR's source on this host (`${STAR_SOURCE}/`):
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |- `Transcriptome_classifyAlign.cpp::alignToTranscript` and `Transcriptome_geneCountsAddAlign.cpp`
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  classify an alignment using **only** aligned blocks (`exons[][EX_G]`, `EX_L`), splice junctions
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  (`canonSJ`), strand (`a.Str`), and multimapping count. **No sequence, no quality, no alignment
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  score.** So evidence class **E1 is a sufficient statistic for STARsolo's Gene/GeneFull semantics
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  by construction** — Gate -1A is close to a tautology, and the doc's E0→E4 ablation has a
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  largely predictable answer.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |- `SoloFeature_collapseUMIall.cpp` collapses UMIs in a `for (iG=0; iG<nGenes; iG++)` loop —
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  **per (cell, gene)**, i.e. annotation-dependent by construction.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |This is the project's actual crux, and it is a direct tension: **the storage win comes from
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |collapsing reads→molecules (~3–4×), but the fidelity target requires a collapse that needs the
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |annotation.** Everything else is downstream of how well annotation-independent grouping
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |approximates per-gene collapse.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**Finding B — prior art demotes two of the five claimed contributions.** Full scan in
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |`docs/prior-art.md` (to be written from the planning scan). Highlights:
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |- **recount3 / Monorail / Megadepth / Snaptron** (Genome Biology 2021; Bioinformatics 2018, 2021)
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  already own "annotation-free archive, replay under new annotation, no BAM, no realignment" —
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  for **bulk**, sample-level, coverage-domain, no UMI notion. 990 TB reads → >150 TB summaries (6.6×).
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |- **Boiler** (NAR 2016) already did lossy-but-quantification-faithful alignment archiving:
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  39–56× compression, 3–5× better than CRAM. Bulk, coverage-domain, molecule-destroying.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |- **The barcoded 10x/STARsolo BAM is the real competitor**, and it is not a paper so it cannot be
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  fenced off by citation. It already stores intervals, blocks, junctions, MD/NM, and corrected
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  `CB`/`UB`. Re-quantification from it is routine today (STARsolo `--readFilesType SAM SE`,
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  featureCounts `--byReadGroup`, cellCounts). This subsumes contributions (1) and (3) as *ideas*.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |- **Annotation-independent UMI grouping is not novel** (UMI-tools, sc-SPLASH BKC, UMI-nea 2025).
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |- **Malva** (Nature, 26 Aug 2026) is real and fresh: ~74M cells, index <10% of source footprint.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  On *this host* its pbmc_1k_v3 index is **375 MB vs 5.0 GB FASTQ (13×)** — a directly comparable
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  number on our exact pilot dataset.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**Two facts that help:** Cell Ranger v8+ makes `--create-bam` a required argument and the community
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |increasingly sets it false, so the competing BAM often no longer exists; and no one has published
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |cross-**version** GENCODE deltas for scRNA-seq — that measurement is itself a contribution and is
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |our first gate.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |### Reframed claim (adopt this; it survives the prior art)
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |> *Monorail/Boiler at molecule resolution for single cell.* The first molecule-resolved,
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |> transcriptome-free evidence summary that preserves the UMI identity coverage-domain archives
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |> destroy, is an order of magnitude smaller than the barcoded BAM that is increasingly no longer
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |> produced, supports indexed annotation-compatibility queries, and comes with the first
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |> quantitative characterization of annotation-replay fidelity.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Contribution priority, revised: **(2) replay-fidelity taxonomy** [strongest, unoccupied] >
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**(4) compact format + indexed compatibility queries** [novel in single-cell] >
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**(5) storage** [must be benchmarked honestly vs CRAM 3.1] > (1) abstraction [demote to "minimal
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |sufficient statistic"] > (3) UMI grouping [reframe as *error bound*, not capability].
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |### Decisions taken (user-confirmed 2026-08-28)
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Decision | Choice |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately ||---|---|
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Primary oracle | **STARsolo**; alevin-fry/simpleaf as a second oracle family after Gate -1 |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Gate order | **Front-load storage**; cheapest-lethal-first, deviating from the doc's M4 placement |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Language | **Rust from the start** |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |---
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |## 2. Host resources (verified, all paths real)
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Need | Resource |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately ||---|---|
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Genome + annotation | `/shared-resources/references/GRCh38.gencode.v49/` — `GRCh38.primary_assembly.genome.fa.gz` (845 MB), `gencode.v49.annotation.gtf.gz`, `tx2gene.tsv`, `download.sh` with exact EBI URLs |
| D0 dataset | `${MALVA_ROOT}/reads/pbmc_1k_v3_S1_R{1,2}_001.fastq.gz` (5.0 GB; R1 28bp = 16 CB + 12 UMI, R2 91nt — 10x 3' v3). **Copy to the project scratch volume** (the source volume is nearly full) |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Whitelist | `${WHITELIST_SOURCE}/3M-february-2018.txt` (uncompressed) |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Reference result | `/fs/cbcb-lab/rob/pbmc_1k_v3_simpleaf/af_quant/alevin/` — validated simpleaf run. **Measured: 15,519,436 assigned UMIs over 71,123 barcodes** |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Malva baseline | `${MALVA_ROOT}/` — source + `indices/pbmc_1k_v3` (**375 MB**) + `quant/pbmc_1k_v3` (pseudoquant h5ad). Already built |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || STAR source (for rule fidelity) | `${STAR_SOURCE}/` |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Rust | `env/setup.sh` pattern: project-local `CARGO_HOME`/`RUSTUP_HOME`/`TMPDIR`. rustc/cargo 1.98 |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Python | Clone the pattern of `env/py312` (pysam 0.24, polars 1.44, pyarrow 25, numpy 2.5, samtools 1.24, seqkit) — reports/plots only |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Compute | 256 cores, 1.5 TiB RAM (1.1 TiB free), **14 TB free on the project scratch volume** |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Network | Full outbound HTTPS: GENCODE/EBI, 10x, crates.io, PyPI, conda-forge, bioconda, GitHub all verified reachable |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**Must install** (conda 26.3.2 at `${CONDA_ROOT}/bin/conda`, channels conda-forge+bioconda):
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |`STAR 2.7.11b` (host only has 2.7.2b, whose `--soloFeatures` supports just `Gene`/`SJ` — no
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |`GeneFull`/`Velocyto`), `gffread`. **Set `CONDA_PKGS_DIRS` to the project scratch volume** — conda's default cache is
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |on the 88%-full 30 GB home NFS. Same for `UV_CACHE_DIR`.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**Host gotchas** (from memory + survey): `~/.cargo`/`~/.rustup` are dangling symlinks, so
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |`source env/setup.sh` is mandatory; default `TMPDIR` is on a nearly full volume; `cmake` is absent.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |---
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |## 3. Repository
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |New sibling directory `${GRAVLAX_PROJECT_ROOT}/`, mirroring the
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |conventions that worked on the prior project (provenance manifest, dated decision records,
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |gate memos with pre-registered thresholds and a verdict table — see
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |the predecessor likelihood-coresets project's `memos/gate0-decision.md` for the format to copy).
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |```
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |env/setup.sh                  # adapted from the coresets project
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |docs/project_plan.md          # the source doc, verbatim
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |docs/prior-art.md             # the scan, with the reframing decision
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |docs/implementation_status.md
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |docs/decision_log.md
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |docs/decisions/YYYY-MM-DD-*.md
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |experiments/manifests/gates.yaml   # provenance: tool versions, digests, exact commands
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |scripts/                      # numbered driver scripts, chr19 fast path via $AIE_CHROM
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |src/  (cargo workspace)
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |results/ · runs/ · annotations/ · memos/
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |data/                         # symlinks + DIGESTS.sha256
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |```
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |### Cargo workspace (`src/`), binary `aie`
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Crate | Responsibility |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately ||---|---|
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || `evidence-io` | Record types; archive writer/reader; interned dictionaries; zstd framing |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || `ingest` | Annotation-free BAM → molecule evidence; annotation-independent UMI grouping (G2) |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || `anno` | GTF parse; annotation compiler → exon models, junction sets, interval index |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || `replay` | Molecule↔feature compatibility (**port of `alignToTranscript`**); UMI resolution; matrix emit |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || `eval` | Oracle comparison at molecule, cell, and feature level |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Deps from crates.io (verified reachable): `noodles` (bam/sam/gtf), `rayon`, `clap`, `serde`,
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |`zstd`, `anyhow`. Python in the env is for figures and report tables only.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**Evidence record (E4, overcomplete — do not optimize before Gate G-E):**
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |`cb:u32(interned)`, `umi:u32(2-bit×12)`, `chrom:u8`, `strand:u8`, `blocks:[(start,len)]`,
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |`junctions:[(donor,acceptor,motif)]`, `nm`, `mapq`, `AS`, `NH`, softclip lens, all placements
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |above threshold, `n_reads`, residual clipped sequence.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**Evidence ladder** (per doc §21): E0 coordinate-only · E1 +blocks/junctions · E2 +edits ·
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |E3 +residual sequence · E4 near-lossless. Per Finding A, expect E1 to suffice for STARsolo
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |semantics; the ablation's real job is to **confirm E0 fails on junction concordance** and to
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |find where E2/E3 become necessary (they should: alevin-fry's selective alignment uses transcript
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |sequence, so the second oracle family is where sufficiency gets genuinely hard).
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |### Fast iteration loop
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |The expensive step (alignment) runs **once** per (dataset × annotation-insertion) config; everything
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |downstream is minutes. The dev loop is therefore *align once, replay forever* — which is also
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |literally the project's thesis. All scripts take `$AIE_CHROM` (default `chr19`, gene-dense,
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |~1400 genes) for a whole-pipeline turnaround in minutes; unset for whole-genome.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |---
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |## 4. Annotation panel
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || ID | Annotation | Purpose |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately ||---|---|---|
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || **A1** | GENCODE v49 (on disk) | Baseline. `v50` also available at EBI |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || **A2-far** | GENCODE v32 (= CellRanger 2020-A era) | Realistic "archive built then, re-quantified now". Headline case |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || **A2-near** | GENCODE v48 | Hard case: adjacent releases |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || **A3** | Perturbed v49, programmatic + recorded ground truth | Unseen features |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || **A4** | Gene / GeneFull / Velocyto spliced-unspliced | Alternative semantics — nearly free from STARsolo |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |A3 changes (≥5 classes, sited where D0 actually has coverage; harvest real cases from the v32↔v49
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |diff so they are biologically genuine): added cassette exon · novel splice junction (seed from
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |observed `SJ.out.tab` novel junctions) · 3′ UTR extension/truncation (**highest signal for 10x 3′
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |data**) · gene split · gene merge · new isoform · antisense gene addition (strand semantics).
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |---
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |## 5. Gate ladder (pre-registered thresholds — freeze before running, per prior-project practice)
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Each gate ends in `memos/gate-<id>.md` with a verdict table and an explicit **GO** / **STOP**.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |### G-A · Novelty and reframing — *substantially complete*
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Record the scan in `docs/prior-art.md` and the reframing in `docs/decisions/`.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**Residual task (must do before writing):** obtain Malva's full text + Extended Data to confirm it
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |does **not** retain per-UMI coordinate/junction evidence. *If it does, contribution (1) collapses.*
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**Verdict: GO with the §1 reframing.**
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |### G-B · Is the capability worth anything? (Day 1–2)
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Run the STARsolo oracle on D0 under A1, A2-far, A2-near. Diff the matrices.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Measure: fraction of UMI mass changing gene assignment; genes with >10% count change; cells with
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |>1% total change; the taxonomy of *why* each changed.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**PASS** ≥3% of UMI mass moves v32→v49 **and** ≥5% of expressed genes change >10%.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**STOP** <1% of UMI mass moves and <2% of genes change >10% — annotation churn does not move the
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |matrix, so there is nothing worth replaying.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |*This is also the missing measurement no one has published; keep it publishable-grade.*
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |### G-C · Alignment invariance (Day 2, parallel with G-B) — **highest-risk cheap gate**
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Align D0 three ways off one base index: (i) no GTF + `--twopassMode Basic`, (ii) sjdb from A1,
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |(iii) sjdb from A2-far. Compare per-read evidence tuples (blocks, junctions, strand, primary placement).
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**PASS** ≥99.0% of gene-assigned reads have identical tuples (i) vs (ii); ≥98.5% (i) vs (iii).
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**STOP** <97% — annotation-aware alignment differs enough that an annotation-free archive has an
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |irreducible fidelity ceiling no amount of stored evidence can fix.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |*Also answers a genuinely interesting question in its own right: does 2-pass annotation-free
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |alignment recover what sjdb insertion provides?*
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |### G-D · Annotation-independent UMI grouping (Day 3)
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |On the *same* A1 alignments, compare G1 (STARsolo per-(cell,gene) collapse, from `CB`/`UB`/`GX`)
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |against G2 (collapse on cell + UMI + genomic locus/strand, per doc §20). Report the fraction of
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |(CB,UMI) keys spanning >1 gene or >1 locus.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**PASS** median per-cell L1 relative error <1% **and — the criterion that matters — the G1↔G2
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |discrepancy is strictly smaller than the A1↔A2 annotation delta measured in G-B.** Our error must
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |sit below the signal we claim to capture.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**STOP** G2 error ≥ annotation delta. Prior art requires this be framed as an *error bound*, not
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |a capability claim.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |### G-E · Storage projection (Day 4–6) — **the gate that killed the last project; do it now**
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Build the real E1 archive for D0 with a real encoder (interned dictionaries, delta coding, zstd).
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Benchmark against: compressed FASTQ (5.0 GB, known) · **CRAM 3.1 of the annotation-free BAM with
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |`--capture-all-tags` so CB/UB are retained** · Genozip · Malva's 375 MB on this same dataset.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**PASS** ≥10× smaller than CRAM 3.1 **and** ≤375 MB.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**MARGINAL** 5–10×: continue only with an explicit written capability argument.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**STOP** <5× vs CRAM 3.1.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |*Sanity check from measured numbers: 15.5M assigned UMIs, but the archive must also retain
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |intronic/intergenic/antisense molecules that gene-level quantification discards — precisely because
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |a future annotation may claim them. Realistic record count is well above 15.5M; at 12 B/molecule
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |×~30M ≈ 360 MB, which is roughly Malva's 375 MB. **The bar is tight, and that is the honest read.***
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |### Gate -1A/B/C · Replay fidelity (Day 7–12) — the source doc's gates, unchanged
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |A1 self-replay → A2 cross-replay → A3 unseen features. Classify **every** disagreement before
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |advancing; the taxonomy is the paper's spine (contribution 2). Report all metrics separately on the
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |annotation-sensitive subset S₁₂; correlation alone is not acceptable.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Targets per doc: >99.5% agreement on high-confidence molecules; small cell×gene L1;
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |high precision/recall on A3 features.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |### Gate 0 · Full compression and query audit (Day 13–16)
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Only if G-E was PASS, not MARGINAL. Repetition structure (coordinate/cell-local/cross-cell),
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |alternative-placement multiplicity, residual burden, replay wall-clock vs fresh STARsolo.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Negative controls per doc §16: matrix-only baseline (must fail), Malva pseudocount baseline
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |(already built on this host), full BAM/CRAM replay.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |---
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |## 6. Datasets
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || ID | Dataset | Status |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately ||---|---|---|
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || **D0-dev** | pbmc_1k_v3, chr19 only | Minutes per iteration |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || **D0** | pbmc_1k_v3 whole genome, ~1k cells, ~66M reads | **On disk** |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || **D1** | pbmc_5k_v3 | Download from 10x (verified reachable) |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || **D2** | Annotation-sensitive tissue (brain/testis — high lncRNA churn) | Only after G-E |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || **D3** | Single-nucleus or 5′ | Only after Gate -1; high intronic content makes A4 velocity semantics a strong test |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |---
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |## 7. Execution order — Day 0
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |1. Create the tree; copy `env/setup.sh` and adapt `PROJ_ROOT`. Copy the source doc to
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |   `docs/project_plan.md` and the prior-art scan to `docs/prior-art.md`.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |2. `conda create` a project-local env on the project scratch volume (**`CONDA_PKGS_DIRS` first**) with
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |   `star=2.7.11b gffread samtools`; verify `STAR --version` and that `--soloFeatures GeneFull Velocyto` parses.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |3. Stage D0 FASTQs and the whitelist onto the project scratch volume; write `DIGESTS.sha256`.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |4. Fetch GENCODE v32/v48/v49 GTFs (v49 already local); record URLs and digests in
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |   `experiments/manifests/gates.yaml`.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |5. Build the base STAR index (no GTF) from `GRCh38.primary_assembly.genome.fa`. ~1 hr at 32 threads,
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |   ~30 GB. Per-annotation sjdb is inserted **at mapping time** (`--sjdbGTFfile`) so every run shares
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |   one base index — which is also the correct experimental control.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |6. Scaffold the cargo workspace; `cargo build` to confirm the crates.io path works end to end.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |7. **Freeze G-B/G-C/G-D/G-E thresholds into `docs/decisions/2026-08-28-frozen-gates.md` before any
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |   measurement is run.**
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |Then G-B and G-C in parallel (independent STAR runs, 256 cores available), G-D, G-E.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**Thread-count discipline:** shared interactive box, no scheduler. Cap at 32–64 threads
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |(`CARGO_BUILD_JOBS=32`, `STAR --runThreadN 32`) as the prior project did.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |---
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |## 8. Verification
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |- **Oracle sanity:** STARsolo A1 gene counts vs the existing validated simpleaf run
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  (`/fs/cbcb-lab/rob/pbmc_1k_v3_simpleaf/`, 15,519,436 UMIs). Different mappers so exact agreement
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  is not expected; large disagreement means the oracle is misconfigured. Gate this before G-B.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |- **Replay unit tests:** hand-built synthetic transcript models exercising each
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  `alignToTranscript` return value (Concordant / ExonIntron / Intron / ExonIntronSpan / incompatible),
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  asserting the Rust port matches the C++ semantics case for case.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |- **Round-trip:** archive write→read is bit-identical; E4→E1 projection is deterministic.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |- **End-to-end:** on D0-dev (chr19), `aie ingest` then `aie replay --annotation A1` must reproduce
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  the STARsolo A1 chr19 submatrix within the Gate -1A threshold.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |- **Provenance:** every gate memo cites `experiments/manifests/gates.yaml` (tool versions, digests,
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |  exact commands), as the prior project did.
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |---
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |## 9. Principal risks
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Risk | Where it surfaces | Mitigation |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately ||---|---|---|
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Annotation-free alignment can't support annotation-specific decisions | **G-C** | Highest-risk cheap gate; 2-pass mode is the mitigation, and if it fails the project stops on day 2 |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Storage advantage vs CRAM 3.1 is <5× | **G-E** | Front-loaded specifically because this killed the previous project late |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || "Just keep the BAM/CRAM" | Framing throughout | Answer with size + query speed + `--create-bam false` adoption, not with novelty of the idea |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || Malva retains per-UMI coordinate evidence | **G-A residual** | Read the Nature Extended Data before writing |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately || UMI grouping error exceeds the annotation signal | **G-D** | Pre-registered as the pass criterion, not an afterthought |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |**Overall:** fidelity is very likely to pass (Finding A makes Gate -1A near-tautological for
| The nearly full source volume holding D0 and STAR source is cleaned | Day 0 | Copy both to the project scratch volume immediately |STARsolo). The project lives or dies on **G-C, G-D, and G-E**, all reachable within one week.
