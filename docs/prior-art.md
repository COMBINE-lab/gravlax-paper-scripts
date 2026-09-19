# Prior-art scan and novelty reframing

**Date:** 2026-08-28 · **Gate G-A** · Scan performed during planning; verdict **GO with reframing**.

## Verdict

Not already done, but the source doc's framing is over-claimed: **two of its five stated
contributions are effectively occupied.** The project proceeds under the reframing in §4 below.

## 1. Malva — the stated prior art, confirmed

León-Periñán D, Karaiskos N, Rajewsky N. *Ultrafast and reference-free sequence discovery in
single-cell data.* **Nature, 26 Aug 2026.** https://www.nature.com/articles/s41586-026-10975-w

A k-mer/sequence-space index over raw single-cell + spatial reads; ~74M cells; ~80% of the Human
Cell Atlas. Reference-free, species-agnostic. Index is **<10% of source storage footprint**.
Its CLI has a `quant` command that *pseudo*quantifies against a reference transcriptome and emits
`.h5ad`, and it accepts a barcode list — but it is a sequence index, not a genome-coordinate or
molecule-interval store. No evidence it retains per-molecule intervals, blocks, junctions, or
per-UMI edit information; no evidence it reproduces STARsolo-equivalent UMI counts.

**On this host** (`${MALVA_ROOT}/`, already built by us on pbmc_1k_v3):
index = **375 MB** vs 5.0 GB FASTQ (**13×**); `quant/pbmc_1k_v3/pseudoquant.h5ad` = 26 MB.
This is a directly comparable baseline on our exact pilot dataset.

**Open risk (must close before writing):** obtain the full text + Extended Data and confirm Malva
does not retain per-UMI coordinate/junction evidence. *If it does, contribution (1) collapses.*

Related and equally recent: **sc-SPLASH / SPLASH2**, *Nat Biotech 2026*
(https://www.nature.com/articles/s41587-026-03084-6). Reference-free discovery on barcoded
scRNA-seq; its **BKC** submodule does barcode/UMI-aware preprocessing ~50× faster than UMI-tools.

## 2. The closest conceptual precedent — bulk annotation-free archives

| Work | Year | What it does |
|---|---|---|
| intropolis | 2016 | Junctions across 21,504 SRA samples, annotation-free |
| Snaptron | 2018, *Bioinformatics* | Indexed query engine over junctions/coverage, tens of thousands of samples |
| **recount3 / Monorail** | 2021, *Genome Biology* | >750k samples; junction + bigWig coverage summaries "compiled **without use of an annotation**" |
| Megadepth | 2021, *Bioinformatics* | The re-quantification engine: summarize coverage against any new GTF, repeatedly |

recount3 states outright that Monorail **"can rapidly re-quantify datasets using new or updated gene
annotations starting from the coverage files, bypassing the need to re-run the aligner or to store
large BAM files."** That is the source doc's contributions (1), (2)-adjacent and (5), in bulk form.

**Number to beat:** recount3 turned 990 TB of compressed reads into >150 TB of summaries (**6.6×**).

**Why the gap is still real:** recount3 is (a) bulk, (b) sample-level not molecule-level, (c) has no
UMI/barcode notion, (d) coverage-domain — so per-molecule assignment is impossible by construction.
No "recount/Snaptron for single cell" appears to exist. But we are now explicitly *"Monorail at
molecule resolution"*, not a new idea.

## 3. Existing single-cell intermediates — all annotation-committed

| Format | Committed at build time? | Re-quantifiable under a new annotation? |
|---|---|---|
| 10x `molecule_info.h5` | **Yes, totally** — rows are (barcode, umi, **feature_idx**, count) | No. Genes already assigned; no genomic coordinate at all |
| kallisto\|bustools BUS | **Yes** — records carry **equivalence class** IDs into the index | Only via a new t2g over the *same* transcript set/index order |
| alevin-fry RAD | **Yes** — header carries `refCount`/`refNames` from the mapping index | t2g regrouping and USA-mode splicing status only; bounded by the splici index. Not genome-coordinate |
| STARsolo `--soloFeatures SJ` | **Partly not** — counts annotated *and* novel junctions per cell | Junction counts alone cannot rebuild a gene matrix: no body coverage, no molecule identity |
| UMI-tools `group` | No — position+UMI | A step, not an archive |
| "sc-RAD" | — | **Does not exist**; phantom reference |

**So the specific object — a molecule-indexed, transcriptome-free, genome-coordinate single-cell
evidence record — is genuinely unoccupied among named formats.**

## 4. THE central threat: the barcoded BAM already is the archive

This is not a paper, which is exactly why it is dangerous — a reviewer raises it in one sentence
and there is no citation to argue with.

A Cell Ranger / STARsolo BAM already stores, per alignment: genomic interval, strand, aligned blocks
and junctions (CIGAR N), edit info (MD/NM), unmapped reads, **and corrected `CB`/`UB` tags**. That is
item-for-item the source doc's list of "annotation-independent evidence, molecule-indexed". And
re-quantification from it without touching FASTQ is **routine practice today**: STARsolo
`--readFilesType SAM SE --soloInputSAMattrBarcodeSeq CR UR`; featureCounts `--byReadGroup`;
cellCounts (*Bioinformatics* 2023); dropEst; velocyto; SCExecute.

**This subsumes contributions (1) and (3) as concepts.** What it does *not* do: be small (tens of GB
per run), support fast annotation→molecule compatibility queries, or come with any characterization
of replay fidelity. **The paper must be framed against the BAM, not against `molecule_info.h5`.**

Two facts that help us:
- Cell Ranger v8+ makes `--create-bam <true|false>` a **required** argument; the community
  increasingly sets it false, so the competing BAM frequently no longer exists.
- *Persistent hindrances to data re-use in single-cell genomics*, bioRxiv 2025
  (10.1101/2025.10.02.680150): only ~40% of GEO scRNA-seq studies provide usable processed counts.
  This is the motivation section, with a citation.

## 5. Compression prior art — the storage claim must be tightened

| Tool | Result |
|---|---|
| **Boiler** (NAR 2016) | Lossy RNA-seq alignment archive preserving what quantification needs: **39–56×**, 3–5× better than CRAM. Bulk, coverage-domain, **molecule-destroying** |
| CRAM 3.1 (2022) | ~40–50% over BAM. Lossless; keeps CB/UB with `--capture-all-tags`. **Our primary baseline** |
| Genozip 14/15 | ~16% better than samtools on BAM; ~6:1 on fastq.gz |
| recount3 | 6.6× |
| Malva | ~10× (13× measured here) |

The source doc's "5–100×" is a wide, weak band: **Boiler already achieved 39–56× vs BAM in 2016.**
The claim needs an explicit, fair baseline (CRAM 3.1 with tags retained) or it reads as cherry-picked.

## 6. Does annotation version actually matter?

Strong evidence exists but is mostly **bulk and cross-*source*** (Ensembl vs RefSeq), not
cross-*version* — a soft spot, since our value hinges on version-to-version deltas.

- Zhao & Zhang, *BMC Genomics* 2015: of 21,598 shared genes, only **16.3% quantified identically**;
  **9.3% differ by ≥50%** (RefGene vs Ensembl).
- Chisanga, Liao & Shi, *BMC Bioinformatics* 2022: effective gene-length differences up to **64-fold**;
  RefSeq transcripts grew **3.6M (2015) → 7.8M (2020)**.
- Singletrome, *Sci Rep* 2025: single-cell-specific — adding 91,215 lncRNA genes changes read
  assignment for lncRNAs *and* protein-coding genes via sense/antisense overlap.

**No published study measures cross-GENCODE-version deltas for scRNA-seq.** That measurement is
Gate G-B, and it is itself a contribution that closes the weakest hole in the pitch.

## 7. Annotation-independent UMI grouping is NOT novel

Position-based collapsing is the default mental model: UMI-tools (*Genome Research* 2017), zUMIs,
umis, gencore, fgbio, Picard. Recent: **UMI-nea** (*Bioinformatics* 2025) — "fast, robust,
**reference-free** UMI deduplication"; **sc-SPLASH BKC** — annotation-free, ~50× faster than UMI-tools.

Counterpoint that defines our job: alevin/alevin-fry deliberately resolve UMIs **within gene
context** because gene-agnostic collapsing merges molecules from distinct genes sharing a UMI+locus
and splits molecules whose reads land at different positions.

**Reframe contribution (3) from "we do it" to "we bound the error it introduces."** That is exactly
what Gate G-D is built to measure.

## 8. Ranked threats

1. **The barcoded BAM + STARsolo-from-BAM / featureCounts `--byReadGroup` / cellCounts.** Not a
   paper, so unanswerable by citation. Subsumes (1) and (3) as ideas; leaves (2), (4), (5) open.
2. **recount3 / Monorail / Megadepth / Snaptron.** Owns the framing, in bulk, at 6.6×.
3. **Boiler.** Pre-empts "compact lossy quantification-faithful alignment archive" at 39–56×.
4. **Malva.** Same domain, two days old, high-profile, <10% footprint, has a `quant` command.
5. **sc-SPLASH/BKC and UMI-nea.** Kill contribution (3) in its original phrasing.

Cite-and-distinguish (non-threatening): alevin-fry/RAD, BUS, `molecule_info.h5`, SICILIAN,
Singletrome, scBaseCount (the brute-force reprocessing alternative we argue against).

## 9. Adopted reframing

> *Monorail/Boiler at molecule resolution for single cell.* The first molecule-resolved,
> transcriptome-free evidence summary that (a) preserves the UMI identity coverage-domain archives
> destroy, (b) is an order of magnitude smaller than the barcoded BAM that is increasingly no longer
> produced, (c) supports indexed annotation-compatibility queries, and (d) comes with the first
> quantitative characterization of annotation-replay fidelity — including the missing measurement of
> how much scRNA-seq count matrices actually change between GENCODE releases.

**Revised contribution priority:**

| Rank | Contribution | Status |
|---|---|---|
| 1 | Replay-fidelity taxonomy across annotation-change classes | **Unoccupied. The paper's spine.** |
| 2 | Compact format + indexed molecule-compatibility queries | Novel in single-cell (Snaptron did junctions, not molecules) |
| 3 | Storage advantage | Defensible only against CRAM 3.1 with tags; benchmark honestly |
| 4 | The abstraction itself | **Demoted** to "the minimal sufficient statistic for annotation replay" |
| 5 | Annotation-independent UMI grouping | **Demoted** to an error bound (Gate G-D) |
