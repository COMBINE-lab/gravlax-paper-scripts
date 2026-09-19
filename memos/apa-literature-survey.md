# APA from 10x 3′ scRNA-seq — literature survey and upgrade plan for the Gravex APA capability

Date: 2026-08-30. Scope: methods for poly(A)-site calling, quantification, and differential
usage from 3′-tagged droplet scRNA-seq, assessed against our index-resident algorithm
(`aie query apa` + `scripts/100_apa_replicate.py` / `scripts/105_ptprc.py`; paper §"Differential
3′-end usage between cell populations"; gates in D24/D32).

## 0. Our algorithm, as implemented (baseline for the comparison)

- **Site calling:** single-linkage clustering of per-molecule 3′-most span-extreme coordinates,
  strand-aware, gap ≤ 24 bp (frozen). No peak-shape model, no cluster splitting, no
  internal-priming filter, no reference sequence consulted.
- **Quantification:** exact per-site, per-cell/per-group UMI counts straight from the index
  (this part is *stronger* than the field: no re-counting, no pseudo-alignment, molecule-exact).
- **Differential test:** per gene, |Δ weighted-mean 3′ position| between two groups
  (MIN_UMIS = 30/group, top-250 expressed genes); significance = exceeding the 95th percentile of
  a *single* random within-population half-split, pooled across genes. No per-gene p-value, no
  dispersion model, no FDR.
- **Validation:** 26/27 cross-dataset replication; 30/35 genes with both dominant sites ≤200 bp
  of GENCODE transcript ends.
- **Index-resident evidence available** (format-spec v1.5): per-molecule span extremes
  (most-contained + most-extended read per junction chain), junction chains, per-molecule cell id,
  read counts, paralog patterns. **Not stored:** soft-clipped residual sequence (E3 "residual
  sequence" is a *planned optional section*, never implemented) — so no direct poly(A)-tail
  evidence; no base-level sequence of any kind.

## 1. Comparative table of published methods

| Tool (year, venue) | Site calling model | Priming-artifact / IP filtering | Soft-clip poly(A) use | Junction aware | Per-cell quantification | Differential test |
|---|---|---|---|---|---|---|
| **Sierra** (2020, Genome Biol) | Gaussian NLS fit to coverage per gene (~600 bp peaks); coverage split into within-/across-junction subsets before fitting | *Annotates* (not auto-filters) A-rich context: 13 consecutive A (≤1 mismatch); AAUAAA motif; optional pre-test filter | No | **Yes** (regtools junctions) | UMI-per-peak count matrix | **DEXSeq** on n=6 random pseudo-replicates of cells; padj<0.01, |LFC|>0.5 |
| **scAPA** (2019, NAR) | HOMER `findPeaks` on cluster-pooled BAMs; wide peaks split | Downstream A-stretch peak filter | No | No | UMI counts per peak per cluster (pooled) | **Chi-square goodness-of-fit** across clusters + FDR; proximal Peak Usage Index (PUI) |
| **scAPAtrap** (2021, Brief Bioinform) | Annotation-free: coverage `findPeaks` **+ `findTails`** — reads with terminal soft-clipped A/T stretches give single-base PA positions, rescuing low-coverage sites | Tail evidence itself discriminates; A-rich filtering downstream | **Yes (core feature)** | No | umi_tools per-cell counts per site | None built in (delegates, e.g. movAPA/DEXSeq) |
| **polyApipe** (Monash, unpub. pipeline) | Sites *only* from reads with soft-clipped poly(A); ends grouped into peaks | Misprime flag from genomic A-rich context at clip site | **Yes (only signal)** | No | UMI counts per peak per cell | **weitrix**: per-gene "shift" scores, precision-weighted linear models, calibrated dispersion |
| **MAAPER** (2021, Genome Biol) | No de-novo calling: probabilistic assignment of "nearSite" reads to **PolyA_DB v3** PAS; Gaussian-kernel read→PAS distance learned on single-PAS genes | Inherited from 3′READS-based DB (internal-priming-free by construction) | No | No | Pseudobulk per group (not per-cell) | **Likelihood-ratio test** (gene-level PAS proportions), Fisher's exact (top-2 PAS), paired F-test; REDu/RLDu, REDi/RLDi indices |
| **SCAPTURE** (2021, Genome Biol) | HOMER-style peak calling, then **DeepPASS** neural net classifies each peak as true PAS from sequence context (prob>0.5) | Learned (replaces hard A-rich rules) | No | Partially (transcript-guided) | PAS-transcript per-cell counts | Seurat-style tests on PAS usage between populations |
| **scUTRquant / scUTRboot** (2024, Nat Commun) | No calling: kallisto-bustools pseudoalignment to truncated **cleavage-site atlas** (Microwell-seq, 200+ cell types, +40% vs GENCODE); sites <200 bp merged (declared unresolvable) | Atlas built with **cleanUpdTSeq** naive-Bayes IP classifier; 30-nt merge | In atlas construction (tail-clipped reads) | Via transcriptome | Per-cell CS-isoform counts | **Cell-label bootstrap** of WUI (position-weighted UTR index) / usage index / Wasserstein distance; BH-FDR, effect floor ΔWUI>0.1 |
| **scDaPars** (2021, Genome Res) | DaPars2 linear-regression coverage breakpoint (2-site model) → PDUI per cell; **NNLS imputation** over kNN cell graph | Weak (coverage model) | No | No | PDUI (distal fraction) per cell | Downstream on imputed PDUI (clustering, Wilcoxon-style) |
| **Infernape** (2023, Genome Res) | Candidate peaks from coverage; **Gaussian mixture** deconvolution of overlapping peaks; peak→PA distance calibrated on single-PA-single-peak genes | Exclude peaks with **≥8 consecutive A in +10..+140 nt** downstream of 3′ edge | Yes (apex refinement) | No | Per-cell counts per PA | Cluster-level tests on usage fractions |
| **SCAPE** (2022, NAR) | Bayesian mixture: each PA site a Gaussian over *insert-size-implied* 3′ ends (μ, σ, π) + uniform noise component; EM | Noise component absorbs artifacts; downstream A filter | No (uses insert length instead) | No | Per-cell posterior site assignment | Downstream (e.g. chi-square on weights) |
| **SCINPAS** (2023, NAR GaB) | Dedup modified to keep distal reads; PAS from reads with **non-templated poly(A) tails**; validated vs PAPERCLIP motif profiles | Tail evidence + motif checks | **Yes (core)** | No | Per-cell site counts | n/a (site catalog focus) |
| **scTail** (2025, Genome Biol) | Uses **read-1-preserved** chemistry: R1 carries the tail junction → single-base PAS; clusters R1 ends; **CNN classifier filters internal-priming peaks** | Learned classifier (essential step) | Yes (R1 tails) | No | R2 reads quantify each PAS per cell | BB/GLM-style differential usage (statBiomed framework) |
| **PASTA / CPA-Perturb-seq** (2024, Cell) | Sites from external catalog/polyApipe-style counting | Upstream | No | No | Per-cell site counts | **Dirichlet-multinomial** background → per-site per-cell Pearson "polyA residuals"; linear models on residuals for supervised DE |
| Newer 2025–26: **spvAPA** (supervised multimodal APA), **metaAPA** (consensus of tools' PAS sets), **scPAISO** (2025 preprint), **scraps** (R1-based high-resolution) | — | — | — | — | — | — |

**Benchmark literature** (Guidelines bioRxiv 2024.11.29.626111; NAR GaB 2025 lqaf056): tool
outputs overlap poorly (no gene consistently called by all); Sierra = fewer sites but high
precision *because* junction-aware; scAPAtrap/soft-clip methods = single-base resolution but
false positives grow with read length; SCAPE best at low coverage; recommendation = combine
evidence types and filter low-confidence sites. Sites <200 bp apart are the common failure mode
(scAPA resolves ~5% of them; scUTRquant refuses to try).

**Validation databases:** **PolyA_DB v3** (NAR 2018; 3′READS-based, internal-priming-free,
~85K human PAS, PSE/RPM usage metrics, cross-species conservation flags), **PolyASite 2.0**
(NAR 2020; Zavolan atlas of clustered 3′-end-seq sites), **PolyA_DB v4** (NAR 2026; adds
long-read support; ~130K PAS shared with v3.2/PolyASite — sites in multiple DBs have the most
canonical motif profiles), Mayr-lab cleavage-site atlas (via scUTRquant). We currently validate
only against GENCODE transcript ends — the weakest of the available references.

## 2. Honest assessment: where we are naive

1. **No internal-priming control at all (biggest weakness).** Every serious tool filters or
   models oligo(dT) mispriming on genomic A-stretches; it is the dominant false-positive source
   and specifically *creates* spurious intronic/exonic "sites". Our 1.78M-site table and the
   D24 discover results (57% of "novel" mass ≤5 kb downstream of 3′ ends) are unfiltered. The
   26/27 replication survives because mispriming replicates too — it is sequence-driven, so
   cross-dataset replication does **not** rebut it; a reviewer who knows this field will ask
   within minutes. The 30/35 annotation-concordance is our only current defense, and it covers
   dominant sites of 35 genes.
2. **Statistical test is a single-permutation heuristic.** One random half-split, one seed, P95
   pooled across genes: no per-gene p-value, no per-gene depth calibration (high-UMI genes have
   tighter nulls — pooling P95 across genes conflates them), no FDR, mean-shift-only. The field
   uses DEXSeq (NB GLM), chi-square, LRT, Dirichlet-multinomial residuals, or cell-label
   bootstrap. Also self-inconsistent: the paper's own headline is that regulation is
   *proportional redistribution among shared termini* (27/30), yet |Δ weighted mean| is blind to
   any redistribution symmetric about the mean and underpowered for the rest.
3. **Site boundaries: single-linkage with a 24 bp gap.** Chains arbitrarily (dense genes can
   percolate into one mega-site), cannot split adjacent PAS, and has no notion of a peak apex or
   cleavage position. Field baseline is at minimum mixture/valley splitting (Infernape, SCAPE)
   or a declared 200 bp resolution floor (scUTRquant). We also don't use our own junction chains
   here — ends on different terminal exons can merge if within 24 bp genomically.
4. **Our "3′ end" is an alignment end, not a cleavage site.** 10x R2 reads stop upstream of the
   poly(A) junction by a library-dependent offset; MAAPER/SCAPE/Infernape explicitly model the
   read-end→PAS distance. We report raw span-extreme coordinates and then compare to annotation
   with a 200 bp tolerance that silently absorbs this offset.
5. **No poly(A)-tail evidence.** The one signal that *proves* a 3′ end (terminal soft-clipped
   A's; scAPAtrap/polyApipe/SCINPAS/scTail) is discarded at ingest — the .aie stores placements,
   not clip residues. This is a real format gap, not just an algorithm gap.
6. Smaller: analysis scoped to top-250 expressed genes with GTF gene spans (fine for a demo,
   but the capability claim is "annotation-free" — gene spans could come from our own discover
   loci); MIN_UMIS=30 is ad hoc; no effect-size floor.

What we are *not* naive about: quantification. Molecule-exact, per-cell, per-site counts with
sub-second access, federated across datasets, is genuinely ahead of every pipeline above (they
all re-scan BAMs or re-pseudoalign; none can do this from a 100 MB artifact). The pitch should
be "the index makes the *best* statistics cheap", and then we should actually use them.

## 3. Prioritized improvements (index-implementable)

### P1 — Real per-gene test: multinomial site-usage + cell-label permutation FDR
- **Sketch:** for each gene, S×2 table of per-site UMI counts by group → G-test/chi-square
  statistic (captures redistribution, not just mean shift). Null: B permutations of group labels
  *over cells* (not one half-split), preserving per-cell UMI totals → per-gene empirical p,
  BH-FDR. Optionally Dirichlet-multinomial (PASTA-style) to absorb cell-level overdispersion —
  fit α by moment matching per gene; or scUTRboot-style bootstrap of a weighted-index statistic
  as the secondary "shift" effect size (keep the weighted-mean as the reported effect, drop it
  as the test).
- **Evidence needed:** per-site per-*cell* counts. **Already in the index** (cell ids are stored
  per molecule; `apa` currently aggregates to groups — add `--per-cell` or permute inside the
  query).
- **Validation:** calibration — many within-population splits must give uniform p (replaces the
  single-seed null); re-run D0/D1 replication at q<0.05 and report it with real FDR.
- **Cost:** small; pure analysis-layer. This is the minimum bar for a methods venue.

### P2 — Internal-priming filter from genomic A-content (no format change)
- **Sketch:** strand-aware, for each called site take the reference window +1..+20 nt (and the
  Infernape-style +10..+140 nt scan) downstream of the site's 3′ edge; flag if ≥6 consecutive A
  or ≥12/20 A (classic Beaudoing/Tian rule), or ≥8 consecutive A in the wide window. Report all
  sites with an `ip_flag`; exclude flagged sites from differential testing by default. Optionally
  annotate canonical PAS hexamer (AATAAA/ATTAAA in −50..−10 of the inferred cleavage zone) as a
  positive-confidence tier (Sierra's annotate-don't-delete philosophy fits our
  "evidence-preserving" story).
- **Evidence needed:** site coordinates (index) + reference FASTA (already on disk for other
  gates). The archive stays sequence-free; the *analysis* consults the genome — same relationship
  the discover→GTF→replay loop already has.
- **Validation:** flagged-site fraction; distance-to-PolyASite-2.0/PolyA_DB-v3 distributions for
  flagged vs unflagged (flagged should be depleted near DB sites and of AATAAA); effect on the
  1.78M-site table and on D24's readthrough mass.

### P3 — Poly(A)-tail evidence bit (format v1.6 candidate; the E3 gap, priced honestly)
- **Sketch:** at ingest, when a read's 3′-terminal soft clip is ≥6 nt and ≥90% A (strand-
  corrected), set a per-molecule `tail` flag (optionally 4-bit clip length, capped). Store as one
  new chunk stream; rANS will crush it (low entropy, mostly 0 — expect ≪0.2 bits/mol, i.e. ~a few
  hundred KB on D0). Sites then classify as **tail-anchored** (≥k flagged molecules; cleavage
  position = flagged molecules' end coordinate, single-base) vs **peak-only**. This is exactly the
  scAPAtrap findPeaks+findTails combination, but index-resident and per-molecule exact — and it
  makes P2 partially learnable from our own data (tail-anchored sites = positive training set for
  an A-content threshold, SCAPTURE-style without the CNN).
- **Evidence needed:** BAM soft-clip at ingest. **Not currently stored** — requires a format bump
  (guarded like every other stream change; byte-identical regression on all four datasets minus
  the new stream).
- **Validation:** distance of tail-anchored cleavage positions to PolyASite 2.0 (expect a sharp
  mode at 0, the SCINPAS/scTail figure); fraction of our current sites that are tail-anchored;
  re-run APA gates with tail-anchored-only sites.

### P4 — Junction-aware site splitting + offset calibration
- **Sketch:** (a) cluster 3′ ends *within* junction-chain compatibility groups (two molecules
  whose chains end on different terminal exons never merge, regardless of genomic gap) — uses
  evidence no other 10x tool has per molecule; (b) inside clusters, split at density valleys
  (count minimum < f·min(flanking peaks), f≈0.1, or a 2-component fit) instead of raw 24 bp
  linkage; (c) calibrate the alignment-end→cleavage offset distribution on single-site genes
  (MAAPER's trick) — against tail-anchored positions once P3 lands, against PolyA_DB before —
  and report site positions with that offset, tightening the 200 bp concordance window.
- **Evidence needed:** junction chains + span extremes + per-site counts. **All already stored.**
- **Validation:** resolution benchmark — fraction of annotated PAS pairs 50–200 bp apart that we
  separate (field baseline: scAPA ~5%, scUTRquant punts); concordance at 100 bp and 50 bp
  windows, not just 200 bp.

### P5 — Validate against real PAS databases, not GENCODE ends
- **Sketch:** download PolyASite 2.0 (human) + PolyA_DB v3.2; recompute the concordance analysis
  as precision/recall vs DB sites (with the P4 offset), stratified by tail-anchored/peak-only and
  ip-flagged; add the PSE-weighted version (sites used in many samples should be hit first).
  One afternoon of work, converts an internal sanity check into the field's standard evidence.
- **Evidence needed:** nothing new.

### P6 (opportunistic) — Federated APA as the differentiator
- With P1's test in place, run the same contrast across all six indexes via `federate` and
  meta-analyze (Stouffer over per-dataset p) — no published 10x APA tool can even phrase a
  six-dataset per-molecule query. Cheap capability-section win; only do it after P1/P2 so the
  statistics survive review.

Suggested order: P1+P2 (analysis-only, unblock the paper's claims) → P5 (evidence upgrade) →
P4 → P3 (format bump, schedule with the next version window) → P6.

## 4. Citations to add (refs.bib-ready)

```bibtex
@article{patrick2020sierra,
  author={Patrick, Ralph and Humphreys, David T. and Janbandhu, Vaibhao and others},
  title={Sierra: discovery of differential transcript usage from polyA-captured single-cell RNA-seq data},
  journal={Genome Biology}, volume={21}, pages={167}, year={2020}}
@article{shulman2019scapa,
  author={Shulman, Eldad David and Elkon, Ran},
  title={Cell-type-specific analysis of alternative polyadenylation using single-cell transcriptomics data},
  journal={Nucleic Acids Research}, volume={47}, number={19}, pages={10027--10039}, year={2019}}
@article{wu2021scapatrap,
  author={Wu, Xiaohui and Liu, Tao and Ye, Congting and Ye, Wenbin and Ji, Guoli},
  title={scAPAtrap: identification and quantification of alternative polyadenylation sites from single-cell RNA-seq data},
  journal={Briefings in Bioinformatics}, volume={22}, number={4}, pages={bbaa273}, year={2021}}
@article{li2021maaper,
  author={Li, Wei Vivian and Zheng, Dinghai and Wang, Ruijia and Tian, Bin},
  title={MAAPER: model-based analysis of alternative polyadenylation using 3' end-linked reads},
  journal={Genome Biology}, volume={22}, pages={222}, year={2021}}
@article{li2021scapture,
  author={Li, Guo-Wei and Nan, Fang and Yuan, Guo-Hua and others},
  title={SCAPTURE: a deep learning-embedded pipeline that captures polyadenylation information from 3' tag-based RNA-seq of single cells},
  journal={Genome Biology}, volume={22}, pages={221}, year={2021}}
@article{fansler2024scutrquant,
  author={Fansler, Mervin M. and Mitschka, Sibylle and Mayr, Christine},
  title={Quantifying 3'UTR length from scRNA-seq data reveals changes independent of gene expression},
  journal={Nature Communications}, volume={15}, pages={4050}, year={2024}}
@article{gao2021scdapars,
  author={Gao, Yipeng and Li, Lei and Amos, Christopher I. and Li, Wei},
  title={Analysis of alternative polyadenylation from single-cell RNA-seq using scDaPars reveals cell subpopulations invisible to gene expression},
  journal={Genome Research}, volume={31}, number={10}, pages={1856--1866}, year={2021}}
@article{kang2023infernape,
  author={Kang, Bowei and Yang, Yalan and Hu, Kaining and others},
  title={Infernape uncovers cell type-specific and spatially resolved alternative polyadenylation in the brain},
  journal={Genome Research}, volume={33}, number={10}, pages={1774--1787}, year={2023}}
@article{zhou2022scape,
  author={Zhou, Ran and Xiao, Xin and He, Ping and others},
  title={SCAPE: a mixture model revealing single-cell polyadenylation diversity and cellular dynamics during cell differentiation and reprogramming},
  journal={Nucleic Acids Research}, volume={50}, number={11}, pages={e66}, year={2022}}
@article{moon2023scinpas,
  author={Moon, Youngbin and Burri, Dominik and Zavolan, Mihaela},
  title={Identification of experimentally-supported poly(A) sites in single-cell RNA-seq data with SCINPAS},
  journal={NAR Genomics and Bioinformatics}, volume={5}, number={3}, pages={lqad079}, year={2023}}
@article{hou2025sctail,
  author={Hou, Ruiyan and Huang, Yuanhua},
  title={scTail: precise polyadenylation site detection and its alternative usage analysis from reads 1 preserved 3' scRNA-seq data},
  journal={Genome Biology}, volume={26}, pages={236}, year={2025}}
@article{kowalski2024pasta,
  author={Kowalski, Madeline H. and Wessels, Hans-Hermann and Linder, Johannes and others},
  title={Multiplexed single-cell characterization of alternative polyadenylation regulators},
  journal={Cell}, volume={187}, number={16}, pages={4408--4425}, year={2024}}
@article{wang2018polyadb3,
  author={Wang, Ruijia and Nambiar, Ram and Zheng, Dinghai and Tian, Bin},
  title={PolyA\_DB 3 catalogs cleavage and polyadenylation sites identified by deep sequencing in multiple genomes},
  journal={Nucleic Acids Research}, volume={46}, number={D1}, pages={D315--D319}, year={2018}}
@article{herrmann2020polyasite,
  author={Herrmann, Christina J. and Schmidt, Ralf and Kanitz, Alexander and Artimo, Panu and Gruber, Andreas J. and Zavolan, Mihaela},
  title={PolyASite 2.0: a consolidated atlas of polyadenylation sites from 3' end sequencing},
  journal={Nucleic Acids Research}, volume={48}, number={D1}, pages={D174--D179}, year={2020}}
```
(Author lists marked "and others" abbreviated; SCAPTURE/SCAPE/Infernape middle-author lists and
the SCAPTURE first-author spelling should be double-checked against PubMed when pasted into
`refs.bib`. polyApipe has no journal paper — cite the GitHub/Monash docs in a footnote if used.
Optional extras: the 2024 benchmark preprint bioRxiv 2024.11.29.626111 ("Guidelines for
alternative polyadenylation identification tools...") and PolyA_DB v4, NAR 2026.)

Where to cite in `main.typ`: the APA section currently cites nothing from this field. Minimum:
one sentence acknowledging the tool landscape (Sierra/scAPA/scAPAtrap/MAAPER/scUTRquant) and one
acknowledging internal priming with our mitigation (after P2 lands), plus PolyASite/PolyA_DB
when P5 replaces the GENCODE-ends concordance.

## 5. Sources

- Sierra: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02071-7 (methods via PMC7341584); repo https://github.com/VCCRI/Sierra
- scAPA: https://academic.oup.com/nar/article/47/19/10027/5566587 ; pipeline https://github.com/ElkonLab/scAPA
- scAPAtrap: https://academic.oup.com/bib/article/22/4/bbaa273/5952304 ; https://github.com/BMILAB/scAPAtrap
- polyApipe/weitrix: https://github.com/MonashBioinformaticsPlatform/polyApipe ; docs https://monashbioinformaticsplatform.github.io/polyApipe/polyApipe.html
- MAAPER: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02429-5 (PMC8356463)
- SCAPTURE: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02437-5 (PMC8353616)
- scUTRquant/scUTRboot: https://pmc.ncbi.nlm.nih.gov/articles/PMC11094166/ ; https://github.com/Mayrlab/scUTRquant ; https://github.com/mfansler/scutrboot
- scDaPars: https://genome.cshlp.org/content/31/10/1856.full ; https://github.com/YiPeng-Gao/scDaPars
- Infernape: https://genome.cshlp.org/content/33/10/1774.full
- SCAPE: https://academic.oup.com/nar/article/50/11/e66/6548409 ; https://github.com/LuChenLab/SCAPE
- SCINPAS: https://academic.oup.com/nargab/article/5/3/lqad079/7269179
- scTail: https://link.springer.com/article/10.1186/s13059-025-03710-7 ; https://github.com/StatBiomed/scTail
- PASTA / CPA-Perturb-seq: https://www.cell.com/cell/fulltext/S0092-8674(24)00645-7 ; https://github.com/satijalab/PASTA ; vignette https://satijalab.org/seurat/articles/pasta_vignette.html
- Benchmarks: Guidelines preprint https://www.biorxiv.org/content/10.1101/2024.11.29.626111v1.full ; NAR GaB 2025 https://academic.oup.com/nargab/article/7/2/lqaf056/8131109 (PMC12076406)
- Databases: PolyA_DB 3 https://dx.doi.org/10.1093/nar/gkx1000 ; PolyASite 2.0 https://academic.oup.com/nar/article/48/D1/D174/5588346 ; PolyA_DB v4 https://academic.oup.com/nar/article/54/D1/D247/8356008
- Newer/adjacent: spvAPA https://pmc.ncbi.nlm.nih.gov/articles/PMC11724721/ ; metaAPA https://academic.oup.com/bioinformaticsadvances/article/6/1/vbag147/8694770 ; scPAISO https://www.biorxiv.org/content/10.1101/2025.08.20.669565.full.pdf ; scraps https://www.biorxiv.org/content/10.1101/2022.08.22.504859.full.pdf
