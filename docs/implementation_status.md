# Implementation status

**Project:** Annotation-independent molecular evidence for scRNA-seq
**Started:** 2026-08-28 · **Host:** analysis-host

## Current position

**The source doc's hard continuation gates A–F are all demonstrated** (A–D in
`memos/gate-1BC-decision.md`; E in `memos/gate-E-decision.md`; F provisionally, ~12× from a BAM
stand-in). The principal unbuilt artifact is the block-structured archive file format — every
design constraint for it is now measured: alternative placements (2.2% of matrix mass), the UMI
adjacency graph (1.2 bits/molecule), the frozen query list, and Gate F's real measurement.

**All five week-1 gates are decided. None fired a STOP; the project continues.** Full D0 (66.6M
reads) is processed under all three annotations, the annotation-free ingest alignment is built, and
the archive itself exists: **19,918,212 molecules in 119.5 MB, 11.76× smaller than CRAM 3.1.**

The storage question that ended the predecessor project is answered affirmatively, and in the
conservative variant that keeps UMIs. Remaining week-1 item: G-E (storage), now unblocked.

**The capability is now demonstrated end-to-end** (`memos/gate-1A-replay.md`): annotation-free
ingest + GTF at query time reproduces fresh STARsolo within **0.552% of UMI mass** at a **4.4×
margin** under the annotation signal, with replay in 55 s vs ~16 min fresh. Getting there surfaced
that STARsolo counts unique-GENE multimappers via secondary alignments (2.2% of matrix mass, D16) —
alternative placements are mandatory archive evidence, and earlier tag-derived baselines need
re-basing. Error decomposition now: barcode correction 0.68pp > alignment 0.15pp > assignment 0.03pp.

**Post-week-1 result (`memos/umi-adjacency-graph.md`): storing the UMI *adjacency graph* instead of
UMI *values* makes the archive 2.05× smaller and 4.05× more accurate simultaneously** — 58.4 MB at
24.06× vs CRAM 3.1, with grouping error 0.232%. Not a trade-off: the graph defers merge decisions to
replay rather than freezing one at ingest. The graph costs 1.20 bits/molecule against 27 for raw
values (19.4×), and is sparse for birthday-problem reasons.

**The headline risk has moved.** Fidelity was expected to be the easy part (decision D4: E1 is a
sufficient statistic for STARsolo's gene assignment). It is not — because the *alignment* is
annotation-sensitive even though the assignment rule is not. The error budget is now:

| Source | Magnitude | Gate |
|---|---:|---|
| Alignment ceiling | 1.13% of gene-assigned reads (v32, operating regime) | G-C |
| Grouping error | **0.23%** of per-cell count mass (graph replay) | G-D |
| **Signal to be captured** | **2.43% of UMI mass** | G-B (far pair) |

The cross-sample junction catalogue was attempted and is a recorded **negative result**
(`memos/junction-catalogue-negative.md`): a D0+D1 catalogue of 185k data-supported junctions LOWERS
evidence-tuple agreement to 98.21% (baseline 98.70%) — more discovered junctions means more
divergence from an annotation-aware oracle. It also turned out not to be needed: the ceiling's
matrix-level footprint is ~0.15pp, and end-to-end fidelity landed at 0.217% without it.

## Gate ledger

| Gate | Question | Status |
|---|---|---|
| **G-A** | Is this novel after Malva? | **GO with reframing** — `docs/prior-art.md` |
| **G-B** | Does annotation version move the count matrix? | **GO** — far pair 2.43% mass moved, **21.8% of genes change >10%** · `memos/gate-B-decision.md` |
| **G-C** | Does annotation-free alignment support the same decisions? | **MARGINAL, no STOP** — 98.87% vs v32 (operating regime), 98.70% vs v49; bar 99.0% · `memos/gate-C-decision.md` |
| **G-D** | How much error does annotation-free UMI grouping add? | **PASS** — **0.232%** via the adjacency graph (0.939% with the frozen grouping; residual was UMI *policy* mismatch, not missing information) · `memos/gate-D-decision.md`, `memos/umi-adjacency-graph.md` |
| **G-E** | Is the archive ≥10× smaller than CRAM 3.1? | **PASS** — **58.4 MB, 24.06×** vs CRAM 3.1 with the UMI-adjacency graph (119.5 MB / 11.76× with raw UMIs) · `memos/gate-E-decision.md`, `memos/umi-adjacency-graph.md` |
| Gate -1A/B/C | Replay fidelity | **ALL PASS** — one archive, 4 annotations, no rebuild: 0.207/0.212/0.217/0.246% (v32/v48/v49/w1), 0 cells >1% everywhere; change capture 99.94% sign, 1.17% delta error; withheld recovery 99.88% vs 91% vanishing from the matrix · `memos/gate-1BC-decision.md` |
| Gate 0 | Full compression + query audit | not started |

Thresholds for all of the above were frozen **before** any measurement, in
`docs/decisions/2026-08-28-frozen-gates.md`.

## Done

- Repo at `${GRAVLAX_PROJECT_ROOT}`, conventions copied from the
  predecessor project (provenance manifest, dated decisions, gate memos with verdict tables).
- `env/setup.sh`: project-local `CARGO_HOME`/`RUSTUP_HOME`/`TMPDIR`/`CONDA_PKGS_DIRS`, thread caps,
  and `AIE_CHROM` for the fast path. Rust homes symlink to the predecessor's (rustc/cargo 1.98).
- Conda env `env/sc`: **STAR 2.7.11b**, samtools 1.23.1, gffread 0.12.9. Verified that
  `--soloFeatures GeneFull Velocyto` parses — the host-wide STAR 2.7.2b supports only `Gene`/`SJ`
  and would have silently limited the annotation panel.
- Data staged to the project scratch volume with `data/DIGESTS.sha256`: pbmc_1k_v3 FASTQs (5.0 GB), 3M whitelist,
  GRCh38 primary assembly. A 4M-read `dev` subsample for minute-scale iteration.
- Annotation panel fetched and validated. **v49 is byte-identical to EBI release_49**
  (md5 `0ef4a024…`), confirming the `/shared-resources` copy is official.
- **Annotation-free base STAR index** built (28 GB, no `--sjdbGTFfile`) — the archive's ingest path.
- Cargo workspace: `evidence-io`, `ingest`, `anno`, `replay`, `eval`, `aie`. **22 tests passing.**
  `aie gate-c` implemented.

## Measurements already in hand (before any pipeline run)

| Quantity | Value | Source |
|---|---|---|
| D0 compressed FASTQ | 5.0 GB | on disk |
| D0 assigned UMIs | 15,519,436 | summed from the validated simpleaf run |
| D0 barcodes / nnz | 71,123 / 5,709,557 | same |
| Malva index, same dataset | **375 MB (13.4×)** | `${MALVA_ROOT}/indices/pbmc_1k_v3` |
| GENCODE v32 → v49 | 60,609 → 78,691 genes (+29.8%); 227,462 → 507,365 transcripts (+123%) | GTF parse |
| GENCODE v48 → v49 | +5 genes; 385,669 → 507,365 transcripts (+31%) | GTF parse |

## Key technical findings

1. **STAR's gene assignment reads no sequence.** `Transcriptome_classifyAlign.cpp::alignToTranscript`
   and `Transcriptome_geneCountsAddAlign.cpp` use only aligned blocks, splice junctions, strand and
   `NH`. So **E1 is a sufficient statistic for STARsolo `Gene`/`GeneFull` by construction** — Gate
   -1A is near-tautological, and the project's risk lies elsewhere.
2. **UMI collapsing is annotation-dependent.** `SoloFeature_collapseUMIall.cpp` collapses in a
   `for (iG=0; iG<nGenes; iG++)` loop, i.e. per (cell, gene). This is the project's crux: the
   storage win comes from collapsing reads→molecules, but the fidelity target needs a collapse
   that consults the annotation. Gate G-D measures the resulting error.
3. Block semantics must split at CIGAR `N` **only**; STAR walks across indel separators when
   matching exons (`alignToTranscriptMinOverlap`), so a `D` must extend a block, not split it.
   Encoded and unit-tested in `ingest::cigar`.

## Deviations from the source plan

- **Gate order.** Storage (G-E) is front-loaded ahead of the full fidelity work, rather than sitting
  at M4. Fidelity is very likely to pass (finding 1); storage is what killed the predecessor
  project, and it killed it late.
- **Contributions (1) and (3) demoted** on prior-art grounds — see `docs/prior-art.md` §9.
- **Per-annotation STAR indices pre-built** (`scripts/05`) instead of inserting junctions on every
  oracle run.

## Log

- **2026-08-28** Project created. G-A closed. Day-0 infrastructure complete; base index built;
  workspace green at 22 tests. Dev-scale oracle validated the STARsolo parameter set — it caught
  that STAR refuses to emit `CB`/`UB` into an unsorted BAM, a failure that would have cost hours at
  full scale. Annotated indices for v49/v32/v48 building in parallel.
