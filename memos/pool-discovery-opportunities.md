# Memo — Pool et al. 2023 deep read: what their reference optimization actually does, and what our discover→replay should steal from it

**Date:** 2026-08-30 · **Sources:** Pool, Põldsam, Chen, Thomson, Oka, *Nature Methods* 20:1506–1515
(2023), **doi 10.1038/s41592-023-02003-w** (note: PMID 37697162; the doi `…-01965-1` circulating in
our notes is wrong — the bib entry in `paper/refs.bib` carries no doi and is otherwise correct);
full Methods read from the CC-BY preprint (bioRxiv 2022.04.26.489449, "Enhanced recovery of
single-cell RNA-sequencing reads for missing gene expression data"); ReferenceEnhancer R package
(github.com/PoolLab/ReferenceEnhancer) and Generecovery repo; follow-on tools surveyed below.
**Our side:** `scripts/102_pool_extension.py` (blanket 3 kb replay), `src/crates/aie/src/querycmd.rs`
(`Discover`, `Apa`), decision log D24/D32/D40/D41, `memos/gate-apa-discover.md`.

---

## 1. What Pool et al. actually do (the details that matter)

Three measured sources of read loss (their numbers, uniquely-mapped reads, Cell Ranger 6.1.2 /
Ensembl v98): mouse MnPO brain 71.8% exonic / 19.5% intronic / 8.7% intergenic of 272M reads;
human PBMC 69.9 / 23.5 / 6.7%. Adding intronic, 3′-intergenic-within-10kb, or both raises detected
genes by +13.6 / +25.8 / +33.6% (mouse) and +19.9 / +23.2 / +39.2% (human). ~25% of *all*
intergenic reads sit within 10 kb of an annotated 3′ end (~2× enrichment over genomic expectation)
— that is their core empirical argument that "intergenic" is mostly unannotated 3′ UTR.

Their optimization is a **three-step annotation edit**, then one ordinary Cell Ranger run:

**(a) Hybrid pre-mRNA reference (intronic reads).** They benchmarked three intron strategies and
found each defective alone: (i) the "traditional pre-mRNA reference" (rewrite every `transcript`
row as a single `exon`) captures intron mass but *destroys exon boundaries* — spliced reads
spanning exon–exon junctions no longer match, and nested genes inside a big gene's span become
multigene and are wiped; (ii) `--include-introns` (≈ STARsolo GeneFull) avoids the overlap blowup
but misses some intronic reads the traditional reference catches. Their hybrid: keep the normal
exonic gene structures **and add one full-transcript-span exon per gene** (skipping genes on the
overlap list), then run with `--include-introns`. Both read classes then match. Note: this is a
*mapping* fix — in our world GeneFull replay already gives the intron mass without any of this.

**(b) Gene-overlap resolution (readthrough / premature-start).** Same-strand gene overlaps in
Ensembl v98: 2,035 mouse genes (6.3%), 5,195 human (14.2%) — multigene reads at these loci are
discarded. Two stereotyped pathologies: **readthrough transcripts** (an upstream gene's transcript
swallows a downstream gene's exons → downstream gene blinded) and **premature-start transcripts**
(a downstream gene's transcript starts 5′ of the upstream gene's *terminal exon* → since 3′-tag
reads concentrate on terminal exons, the *upstream* gene is blinded — the worst case). Dozens of
gene pairs *share* a terminal exon and are mutually invisible. Their resolution is a decision tree
(automated triage → manual IGV curation): exons don't overlap → keep both; one of the pair is a
weak model (mouse `Gm*`/`*Rik`, human `AC*`/`AL*`) with exon overlap → delete the weak model;
both well-supported with exon overlap → count exon-overlap incidences per exon: downstream gene's
exon is the offender → delete its premature-start transcript(s), upstream gene's exon is the
offender → delete its readthrough transcript(s); fully shared terminal exon → delete one gene and
rename the survivor. Deletions are confirmed against RefSeq (delete only if unsupported in both
Ensembl and RefSeq) and against observed read mapping (zero reads at the locus).

**(c) 3′ boundary redefinition (false intergenic reads).** Extract unique, barcode-valid intergenic
reads from the BAM (with a documented Cell Ranger ≤6 quirk: reads antisense to exons carry
`RE="E"` + `AN` tag and are really intergenic — they re-classify these); `bedtools closest -s -D a
-fu` assigns each to the nearest same-strand upstream gene end; **rank genes by intergenic read
count within 10 kb of their 3′ end**; then *manually* accept an extension per gene if at least one
of: (i) reads **splice into known exons** of the gene (their strongest criterion), (ii) continuous
coverage from the annotated end, (iii) **RefSeq already ends further 3′** (then adopt the RefSeq
end wholesale, done genome-wide), (iv) orthogonal expression truth (Allen ISH atlas, Human Protein
Atlas). The new end is eyeballed in IGV, recorded per gene. Scale: they "convincingly" extended
**under 1,000 genes per species** — supervised curation is their bottleneck, stated as such.

**Assembly & results.** An R script deletes listed transcripts/genes, rewrites 3′ ends, appends
the hybrid pre-mRNA exons, and emits one GTF (`cellranger mkref` + `count --include-introns`).
Gains: mouse >3,000 newly detected genes, +14.8% registered reads, ≈+600 genes/cell median (>20%),
13→16 neuron types (recovering Ptger3 → identifies the warmth-activated type); human >4,500 genes,
+21% reads, ≈+400 genes/cell, up to 6 additional T-cell subtypes with canonical markers (Tregs,
MAIT, Th2/Th17…). Tools: Cell Ranger 6.1.2, STAR 2.7.9a (GeneFull with CR-matched flags),
bedtools 2.30.0, GenomicAlignments, IGV. The published version packages steps (b)/(c) triage as
**ReferenceEnhancer** (`IdentifyOverlappers` → `OverlapResolutions` → `IsolateIntergenicReads` →
`GenerateExtensionCandidates` (10 kb rank order) → `OptimizedAnnotationAssembler`).

**Follow-on work (post-2023).**
- **GeneExt** (Zolotarov, Grau-Bové, Sebé-Pedrós; Bioinformatics btaf094-line, preprint bioRxiv
  2023.12.05.570120): the automated version of step (c) for any species. MACS2 stranded peaks from
  the BAM; *genic* peak coverage distribution sets a percentile threshold (`--peak_perc`) that
  filters *intergenic* peaks; each gene extends to its **most downstream surviving peak** within
  `-m` (recommended ≈1–2× median gene length); **orphan peaks** (no gene within reach) are kept and
  clustered into putative unannotated genes (inter-peak gap ≤ 75th-pct intron size, cluster span
  ≤ median gene length); optional 5′ clipping of the downstream gene at overlaps.
- **peaks2utr** (Haese-Hill et al., Bioinformatics 39:btad112, 2023): MACS2 peaks → formal 3′ UTR
  features added to a canonical annotation; uses soft-clipped polyA evidence in reads to place ends.
- **10x's own reference update** (GRCh38/GRCm39 "2024-A" bundles) adopted the same lesson: filters
  readthrough transcripts and problematic overlapping models at the vendor level; Cell Ranger ≥7
  made `--include-introns` the default. Pool's step (b) is being absorbed upstream — worth saying
  in the paper.
- **GENCODE 2025** (Mudge et al., NAR 53:D949) cites Pool for the 3′-end deficit; their fix is
  long-read (CLS/polyA-seq) re-annotation — annotation authorities are the slow path, per-dataset
  evidence the fast one. **scCensus** (He, Mount, Patro, bioRxiv 2024) systematizes what off-target
  (intergenic/antisense/intronic) reads contain; **Forseti** (He et al., Bioinformatics 40, ISMB
  2024) models splicing status of individual reads; **ELATUS** (Goñi et al., Nat Commun 2024)
  shows functional lncRNAs hiding in exactly the read classes standard pipelines drop.

---

## 2. Gap analysis

**What Pool exploits that our discovery does not (yet):**

1. **Per-gene evidence gating of 3′ extension.** Our `102_pool_extension.py` extends *every*
   terminal exon by a flat 3 kb clipped at the next gene span (any strand). Pool ranks genes by
   observed downstream mass in a 10 kb window and accepts each extension on evidence; GeneExt picks
   the most downstream *supported peak* as the new cleavage site. Consequences of our bluntness,
   already visible in our own results: (i) the paper's admitted nuclei caveat — extending every
   isoform's terminal exon sweeps in intronic mass near internal isoform ends (the +126% on nuclei
   blends true 3′ recovery with pre-mRNA capture); (ii) 3 kb under-reaches for exactly the
   NRXN1/RBFOX1 class Pool targets (their window is 10 kb; some validated ends are further).
2. **Splice-linkage as attachment evidence.** Pool's strongest manual criterion — downstream reads
   splicing into known exons — is precisely what our junction catalogue holds queryably, yet
   `Discover` never consults junctions: candidates are span-clustered molecules, emitted as
   single-exon models, never *attached to* the upstream gene they extend. Our own DISC-3 measured
   the cost: 57% of "novel" no-overlap mass sits ≤5 kb downstream of an annotated 3′ end — Pool
   would call nearly all of it unannotated UTR of the neighbor, not novel loci.
3. **Overlap-pathology resolution.** We replay whatever annotation we're handed; readthrough and
   premature-start transcripts blind genes in our exact-replay matrices just as they do in Cell
   Ranger (STARsolo Gene semantics discards multigene reads — we reproduce the discard perfectly).
   Our D33 audit already priced this: 3.82% of read mass multigene-discarded, 691k gene-informative
   classes fully dropped. And `Discover`'s claim-generously rule (any-strand span overlap) means
   overlap pathologies also *suppress discovery* (the 82.7% stratum in DISC-1).
4. **Orthogonal-truth validation per edit** (RefSeq cross-reference, ISH/Protein Atlas). We
   validate in aggregate (withheld panel, annotation history); we don't yet score individual
   discovered ends against external 3′-end truth (polyA atlases, RefSeq ends).
5. **Internal-priming awareness.** Pool cites the finding that much intronic signal is aberrant
   priming at intronic polyA tracts; GeneExt filters low-coverage intergenic peaks partly for this
   reason. We currently apply no A-content filter to `apa` sites or discover candidates.

**What we exploit that Pool (and GeneExt/peaks2utr) cannot:**

- **No re-alignment, ever.** Their loop costs a Cell Ranger run per annotation iteration; ours is a
  7–35 s replay — so evidence thresholds can be *swept and validated*, which is why our discovered
  loci get quantified better than the clustering that found them (7.6% vs 26.3%).
- **Genome-wide unbiased discovery**, including truly intergenic loci and the cell dimension
  (per-candidate UMI *and cell* counts); federation gives replication as a confidence signal
  (79.4% of D0 intergenic loci recur in D1) — no analogue in any of these tools.
- **Junctions and 3′ sites are already indexed** (jpost genome-wide support totals; `apa` 24 bp
  site clustering) — Pool needed IGV sessions per gene for the same evidence.
- **Intron handling is free**: GeneFull replay = their step (a) without touching the reference.

---

## 3. Concrete improvements (prioritized)

### I1 — Evidence-gated 3′ extension: `aie query extend` (replaces the blanket 3 kb) — HIGH, ~1–2 days
**Compute from the index:** one genome scan (reuse the `Discover` chunk walk). For each annotated
gene g on strand s: collect unclaimed molecules in (end(g), end(g)+W], W=10 kb, same strand;
compute (i) total UMIs/cells, (ii) the contiguity profile — largest zero-evidence gap between the
annotated end and each candidate stopping point, (iii) 3′-site clusters in the window (`apa`
machinery, 24 bp gap), (iv) junction linkage: any junction with donor inside g's terminal exon and
acceptor in the window (jpost lookup, zero-decode). **Rule:** new end = downstream edge of the most
downstream 3′-site cluster with ≥U UMIs and ≥C cells (U=5, C=3 to start) reachable from end(g)
without an evidence gap >G (G=2 kb), clipped at the next *same-strand* gene start minus 1 (keep the
any-strand clip as an optional stricter mode); extend **only the 3′-most isoform's terminal exon**
(not every isoform — this alone removes the nuclei intronic-blending caveat we currently footnote).
Emit the optimized GTF plus a per-gene evidence table (Pool's 'gene extension candidates' file,
computed instead of curated). **Benefit:** per-gene validated extension with longer reach than 3 kb
where evidence supports it and *no* extension where it doesn't; turns the paper's "we note this
rather than tune it away" caveat into a designed result; directly comparable to Pool/GeneExt.
**Validate:** (a) distance from chosen ends to PolyASite 2.0 / PolyA_DB v3 clusters vs the blanket
3 kb ends and vs unextended v49 ends; (b) agreement with RefSeq ends where RefSeq is 3′ of GENCODE
(Pool's own criterion (iii)); (c) replay: retain ≥90% of the blanket run's recovered mass on the
Pool marker genes (NRXN1, RBFOX1, KCNIP4…) while adding measurably less intronic/internal mass on
nuclei (score: recovered mass whose molecules have junction/contiguity linkage to the terminal
exon); (d) D0↔D1 replication of chosen end positions (±200 bp, the D41 convention).

### I2 — Junction-stitched candidate models and splice-attachment in `discover` — HIGH, ~2–4 days
**Compute:** for each candidate cluster, gather member molecules' junction chains. (i) *Stitch*:
merge clusters connected by a junction (donor in one, acceptor in the other) regardless of the
1,000 bp merge gap — data-driven gene span instead of fixed-gap single linkage; (ii) *structure*:
emit multi-exon models where junctions partition the span (exons = covered blocks between
junctions), instead of the current single-exon `AIENOVEL` span; (iii) *attach*: if a candidate's
junctions splice into an annotated gene's exons, classify it as `extension_of=<gene>` (feeding I1)
rather than novel. **Benefit:** automates Pool's strongest manual evidence; splits DISC-3's
57%-downstream bucket into "UTR extension of neighbor" vs "novel" *by mechanism, not distance*;
should cut the 26.3% boundary-driven quantification error and the 58% median error on v32→v49
novel lncRNAs (their structures are multi-exon; we currently hand replay a rectangle); emitted GTFs
become usable for velocity replay. **Validate:** withheld panel — exon-level Jaccard of emitted vs
true models, plus DISC-2 rerun (expect median err < 26.3%); v32→v49 prospective discovery quant
error (expect < 58%); fraction of downstream-of-3′ candidates that junction-attach (report as the
honest precision split).

### I3 — Annotation-pathology audit + transcript-blacklist replay (Pool step (b) as a query) — MED-HIGH, ~2 days
**Compute:** for each transcript in the supplied GTF: (i) its junction set's index support from
jpost (a readthrough transcript's *bridging* junction — spanning the gap between two gene loci —
with 0 observed support is the smoking gun); (ii) presence of an `apa` 3′-site cluster at its
annotated end; (iii) structural flags computable from the GTF alone: transcript overlapping ≥2
gene loci (readthrough), transcript starting 5′ of an upstream gene's terminal exon
(premature-start), fully shared terminal exons. Emit a ranked blacklist (structural flag AND zero
distinctive index evidence); replay with the blacklist removed; report per-gene recovered UMIs.
**Benefit:** recovers part of the 3.82% multigene-discarded mass with *evidence per deletion*
(Pool needed IGV; 10x's 2024-A does it by fiat); complements the EM layer (EM shares ambiguous
mass; this removes the ambiguity at its annotation source). **Validate:** overlap of our
structural readthrough calls with GENCODE `readthrough_transcript` tags and 10x 2024-A exclusion
lists (external truth for the classifier); replay recovered genes vs Pool's recovered PBMC marker
set on D0; regression: blacklist-replay deviation on non-overlapping genes must be 0.

### I4 — Internal-priming flags on 3′ sites and candidates — MED, <1 day
**Compute:** for every `apa` site cluster and every discover candidate 3′ end, A-fraction of the
genomic +1..+20 window (strand-aware, from `ref/`); flag ≥12/20 A or ≥6 consecutive A (the
scAPA/polyApipe convention). Report flags in `apa` and `discover` output; never silently drop.
**Benefit:** quantifies how much of the 3,151 intergenic loci and of APA-diff signal could be
priming artifact — the reviewer question we will get; sharpens the D41 3′-concordance claim.
**Validate:** flagged fraction at PolyASite-matched sites vs unmatched vs random intergenic
positions (expect low/mid/high); D0↔D1 replication rate of flagged vs unflagged loci (real sites
should replicate more — if flagged loci replicate equally, that is itself a finding to report).

### I5 — Strand-resolved claiming: an antisense candidate class — MED, ~1 day
**Compute:** split `Discover`'s claim into same-strand (claimed) vs antisense-only overlap; report
antisense-overlap molecules as a separate `antisense` candidate class (with the sense gene named)
instead of discarding them into "claimed". Keep the conservative rule for the headline numbers.
**Benefit:** opens the lncRNA/antisense discovery surface (ELATUS: functional lncRNAs live exactly
here; Pool documented Cell Ranger's antisense-to-exon misclassification), and partially addresses
the DISC-1 82.7% overlap stratum without weakening the frozen claim rule. **Validate:** v32→v49
antisense/lncRNA additions overlapping v32 genes on the opposite strand (currently absorbed by the
conservative rule — this is the stratum where our prospective-discovery recall is untested);
D0↔D1 recurrence of antisense candidates.

### I6 — (byproduct) Ship the ranked extension-candidate table as a paper artifact — LOW, free with I1
Pool's `gene_extension_candidates.csv` per dataset, from one index scan, no BAM: gene, downstream
UMIs/cells within 10 kb, contiguity, best 3′ site, junction linkage, priming flag. It is the
"reference optimization as a replay" section's figure-ready evidence and a deliverable no
BAM-discarding pipeline can produce retroactively.

---

## 4. Citations to add to `refs.bib`

- **Keep** `pool2023recovery` (entry correct; if we add dois anywhere, theirs is
  10.1038/s41592-023-02003-w — *not* `…-01965-1`).
- Zolotarov G, Grau-Bové X, Sebé-Pedrós A. *GeneExt: a gene model extension tool for enhanced
  single-cell RNA-seq analysis.* Bioinformatics (2026; preprint bioRxiv 2023.12.05.570120). — the
  automated peak-based competitor for I1/I2; we should compare against it.
- Haese-Hill W, Crouch K, Otto TD. *peaks2utr: a robust Python tool for the annotation of 3′
  UTRs.* Bioinformatics 39(3):btad112 (2023).
- Herrmann CJ et al. *PolyASite 2.0.* Nucleic Acids Res 48:D174 (2020) — external 3′-end truth for
  I1/I4 (and/or Wang R et al., *PolyA_DB 3*, NAR 46:D315, 2018).
- Mudge JM et al. *GENCODE 2025.* Nucleic Acids Res 53:D966 (2025) — annotation authorities citing
  Pool's deficit; frames "annotations churn, indexes shouldn't".
- Optional, discussion-grade: He D, Mount SM, Patro R, *scCensus* (bioRxiv 2024) — taxonomy of
  off-target read information; Goñi E et al., *ELATUS*, Nat Commun (2024) — functional lncRNAs in
  discarded read classes (motivates I5); 10x Genomics 2024-A reference note (readthrough filtering
  adopted by the vendor — motivates I3).
