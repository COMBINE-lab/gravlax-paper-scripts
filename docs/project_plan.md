# Project Plan: Annotation-Independent Molecular Evidence for Single-Cell RNA-seq

**Date:** 2026-08-28  
**Status:** New primary feasibility candidate  
**Primary goal:** Develop and evaluate a compact intermediate representation between raw single-cell sequencing reads and already-quantified count matrices that preserves enough molecule-level evidence to support accurate re-quantification against future or alternative annotations without returning to FASTQ.

---

# 1. Executive summary

Standard single-cell RNA-seq workflows collapse information aggressively:

\[
\text{FASTQ}
\rightarrow
\text{mapping}
\rightarrow
\text{molecular resolution}
\rightarrow
\text{cell}\times\text{feature count matrix}.
\]

The final count matrix is compact and convenient, but it is tightly coupled to the annotation and feature definitions used during the original analysis. Once the raw sequencing evidence has been reduced to counts, changing the annotation, transcript definitions, spliced/unspliced interpretation, or feature catalogue can require remapping or reprocessing the original reads.

Recent work such as **Malva** demonstrates that retaining sequence-level information between FASTQ and the count matrix is scientifically valuable. Malva enables sequence, mutation, junction, pathogen, and other arbitrary queries over extremely large single-cell collections. However, its quantitative outputs are sequence-derived pseudocounts rather than UMI-resolved molecule quantification, and it is fundamentally a sequence-search index rather than a representation designed to replay annotation-specific molecular inference.

This project asks whether there is a useful missing representation:

\[
\boxed{
\text{raw reads}
\rightarrow
\text{annotation-independent molecular evidence}
\rightarrow
\text{arbitrary future quantification}
}
\]

The central target is a compact archive that:

1. is built once from raw single-cell sequencing data;
2. preserves cell barcode and molecule/UMI identity;
3. avoids committing to a specific transcript or gene annotation at ingest;
4. preserves enough sequence and/or genomic event evidence to later determine molecule–feature compatibility;
5. supports re-quantification against alternative or updated annotations;
6. reproduces fresh annotation-specific processing as closely as possible;
7. is substantially smaller and faster to reuse than raw FASTQ/BAM;
8. can eventually scale to very large cell atlases.

The initial project must **not** start by designing a sophisticated compressed index.

Instead, begin with an intentionally overcomplete and inefficient molecular evidence representation on modest datasets. The first question is information sufficiency:

> Can an annotation-independent molecule record reproduce the molecular assignments and counts obtained by fresh processing against multiple different annotations?

Only if that succeeds should the project move to data-structural compression and large-scale indexing.

The project is attractive as a methods contribution if it can ultimately provide either:

- a fundamentally new capability: near-faithful annotation replay without raw reads; or
- a dramatic efficiency advantage: orders-of-magnitude faster and/or substantially smaller reanalysis than remapping/reprocessing.

---

# 2. Novelty boundary after Malva

The recently published Malva work changes the novelty claims that are available.

We should **not** claim:

- the first annotation-independent single-cell sequence index;
- the first ability to query arbitrary sequences across large cell collections;
- the first sequence-level representation between reads and count matrices;
- the first reference-free single-cell search system.

Malva already establishes those capabilities at very large scale.

The proposed project must instead focus on the distinction between:

\[
\text{sequence search}
\]

and

\[
\text{molecular re-quantification}.
\]

The proposed representation should treat the **captured molecule / UMI** as the fundamental unit, not merely the occurrence of a sequence k-mer in a cell.

A successful system should answer questions such as:

> Given a new annotation \(A'\), what molecules in each cell are compatible with each feature in \(A'\), and what count matrix would a fresh quantification pipeline have produced?

This is a stronger target than approximate sequence-derived pseudocounts.

---

# 3. Core research question

Let:

- \(D\) be the original sequencing experiment;
- \(A\) be an annotation or feature catalogue;
- \(Q(D,A)\) be the result of a trusted fresh annotation-specific quantification pipeline.

We seek an annotation-independent representation:

\[
E(D)
\]

such that for a useful family of annotations \(\mathcal A\),

\[
\hat Q(E(D),A)
\approx
Q(D,A)
\qquad
\forall A\in\mathcal A.
\]

The desired representation should be much smaller and cheaper to query than retaining and replaying the complete raw sequencing experiment.

The first project phase should answer:

1. What information must \(E(D)\) contain?
2. How faithfully can it reproduce fresh quantification?
3. Which annotation changes can be replayed exactly or nearly exactly?
4. Which changes fundamentally require information not preserved in \(E(D)\)?
5. How much redundancy exists between read-level evidence and molecule-level evidence?
6. Is there enough redundancy to justify a specialized compact representation?

---

# 4. Fundamental design principle

The project should distinguish:

## Reference independent

The representation requires no genome/reference at ingest.

This maximizes flexibility but likely requires storing substantial sequence evidence.

## Annotation independent

The representation may use a genome/reference at ingest but must not commit to gene or transcript annotations.

This is the preferred initial target.

A genome is relatively stable compared with annotations, and genomic coordinates, orientation, splice junctions, and local sequence differences are highly compressible.

Therefore the initial project should prioritize:

> **genome-aware but annotation-independent molecular evidence.**

Reference-free residual sequence evidence can be added later for cases not faithfully representable by genomic events.

---

# 5. Candidate evidence representation

For each inferred molecule \(m\), define a record:

\[
M_m =
(c_m,
u_m,
R_m,
E_m)
\]

where:

- \(c_m\): corrected cell barcode;
- \(u_m\): corrected or grouped UMI identity;
- \(R_m\): read/molecule metadata;
- \(E_m\): annotation-independent sequence/genomic evidence.

A deliberately overcomplete initial \(E_m\) may contain:

- genomic mapping interval(s);
- orientation;
- exon-like aligned blocks;
- splice junction coordinates;
- alignment score / edit summary;
- start/end positions;
- local indel or mismatch information;
- sequence witnesses/minimizers for portions not explained by the genome;
- all alternative genomic placements above a permissive threshold;
- per-placement confidence or score;
- number of raw reads supporting the molecule.

Do not optimize this representation during Gate -1.

The objective is first to determine what is sufficient.

---

# 6. Two-layer long-term representation

If feasibility is established, the likely final design is a two-layer archive.

## Layer 1: genomic/event evidence

For evidence well explained by the reference genome, store:

- genomic intervals;
- strand;
- splice junctions;
- coordinate transforms;
- local variant/edit information;
- molecule identity.

This should be compact and efficient.

## Layer 2: residual sequence evidence

Store sequence-derived witnesses only for evidence not adequately explained by Layer 1, including:

- unmapped sequence;
- insertions;
- fusion breakpoints;
- pathogen sequence;
- novel sequence;
- poorly mapped fragments;
- potentially annotation-relevant sequence variants.

This layer enables later discovery and increases robustness to future annotations.

The initial feasibility work should measure how much evidence falls into each layer.

---

# 7. Quantification replay model

Given an annotation \(A\), define feature set:

\[
F_A.
\]

Each feature \(f\in F_A\) consists of some combination of:

- genomic intervals;
- splice junctions;
- strand constraints;
- transcript path constraints;
- sequence-specific constraints.

For each stored molecule \(m\), evaluate:

\[
C_A(m)
=
\{f\in F_A:
E_m\text{ is compatible with }f\}.
\]

The replay system should construct molecule–feature compatibility sets without remapping raw reads.

Depending on the target quantification model, it can then perform:

- winner-take-all assignment;
- equivalence-class counting;
- UMI resolution;
- EM / probabilistic assignment;
- gene-level counting;
- transcript-level counting;
- spliced/unspliced counting.

The first experiments should use the simplest target necessary to test replay fidelity.

---

# 8. Scope of the first feasibility phase

Do not start with massive atlases.

Use modest datasets that permit repeated complete reprocessing.

Recommended scale:

- 5k–20k cells;
- standard 10x Chromium 3' datasets;
- ideally 50M–250M raw reads or less for the earliest pilots;
- one small enough subset for extremely rapid iteration.

Initial development sequence:

1. 500-cell or 1k-cell subset for code/debugging;
2. 5k-cell dataset for Gate -1;
3. 10k–20k-cell confirmation dataset only after the design passes.

Do not perform atlas-scale engineering until the replay capability is proven.

---

# 9. Annotation panel

Feasibility requires genuinely different annotations.

At minimum construct:

## A1 — baseline annotation

For example, GENCODE v49 or v50 using the standard lab preprocessing conventions.

## A2 — nearby real annotation release

Use a newer/older GENCODE or Ensembl release with meaningful differences in:

- transcript models;
- exon boundaries;
- gene definitions;
- 3' UTRs;
- overlapping genes.

## A3 — deliberately perturbed annotation

Construct a controlled annotation containing synthetic or curated changes:

- added exon;
- removed exon;
- alternate 3' end;
- new splice junction;
- split gene;
- merged gene;
- new transcript isoform.

The key property is that these features were **not known to the molecular archive at ingest**.

## A4 — optional alternative feature semantics

Examples:

- gene-only quantification;
- spliced/unspliced gene-state features;
- transcript-level features;
- custom targeted feature definitions.

---

# 10. Trusted fresh-processing oracle

Every replay result must be compared against a fresh processing oracle.

For each annotation \(A_i\):

\[
Q_i^{\mathrm{fresh}}
=
Q(D,A_i).
\]

Use an established trusted pipeline, preferably one already maintained by the lab, and freeze:

- mapper;
- barcode correction;
- UMI rules;
- score thresholds;
- transcript-to-gene rules;
- quantification policy.

The archive is evaluated by:

\[
Q_i^{\mathrm{archive}}
=
\hat Q(E(D),A_i).
\]

The central comparison is:

\[
Q_i^{\mathrm{archive}}
\quad\text{vs}\quad
Q_i^{\mathrm{fresh}}.
\]

Do not use correlation alone.

---

# 11. Gate -1: does the required information exist?

This is the first and most important feasibility gate.

Build an intentionally overcomplete archive.

No compression.
No clever index.
No final file format.

Then perform annotation replay.

## Gate -1A: baseline self-replay

Build \(E(D)\) without using annotation \(A_1\).

Replay \(A_1\).

Measure agreement against fresh \(A_1\).

Required initial target:

- >99.5% agreement for high-confidence uniquely assignable molecules;
- very high agreement in gene-level UMI counts;
- discrepancies fully classifiable.

If even the baseline annotation cannot be replayed accurately, stop and diagnose the missing evidence.

---

## Gate -1B: cross-annotation replay

Replay \(A_2\) from the same archive.

Compare with fresh \(A_2\).

The archive must not be rebuilt.

Measure:

- molecule compatibility agreement;
- cell × gene count agreement;
- total molecules assigned;
- molecules whose assignment changes between \(A_1\) and \(A_2\);
- whether those changes are reproduced correctly.

The important subset is:

\[
\{m:Q(D,A_1)\neq Q(D,A_2)\}.
\]

Overall agreement can be misleading because most molecules may be unaffected by annotation changes.

Require very strong accuracy specifically on annotation-sensitive molecules.

---

## Gate -1C: unseen-feature replay

Replay controlled annotation \(A_3\).

Ask whether molecules are correctly associated with features that did not exist at ingest.

This is the strongest evidence of true annotation independence.

For every synthetic/curated change measure:

- recall of supporting molecules;
- false support;
- count error;
- affected-cell agreement.

If the representation cannot recover genuinely new features without returning to raw reads, characterize exactly why.

---

# 12. Gate -1 stop conditions

Stop the project immediately if:

1. near-faithful replay requires storing essentially the complete read sequence/alignment record;
2. meaningful annotation changes cannot be reproduced from annotation-independent evidence;
3. UMI/molecule grouping itself fundamentally depends on the final annotation in a way that cannot be replayed;
4. archive-based replay differs materially from fresh processing on annotation-sensitive molecules;
5. the amount of evidence that must be retained is nearly as large as BAM/FASTQ with no compensating capability advantage.

A partial capability may still be useful, but it should not remain the primary target without a stronger reformulation.

---

# 13. Gate 0: measure the information/compression opportunity

Only after Gate -1 passes.

Measure the redundancy of the deliberately overcomplete archive.

For dataset \(D\), report:

- raw FASTQ size;
- compressed FASTQ size;
- BAM/CRAM size;
- RAD-like intermediate size;
- number of reads;
- number of aligned reads;
- number of molecules;
- reads per molecule;
- genomic intervals per molecule;
- junctions per molecule;
- residual sequence bytes per molecule.

Important ratios:

\[
R/M
\]

where \(R\) is usable reads and \(M\) is molecules,

and:

\[
\frac{\text{read-level evidence bytes}}
{\text{molecule-level evidence bytes}}.
\]

Do not infer a storage win merely from reads-per-UMI.

Measure actual entropy/repetition.

---

# 14. Gate 0 structural questions

## 14.1 Coordinate repetition

How often do molecules share:

- exon blocks;
- splice junctions;
- alignment path shapes;
- genomic intervals?

## 14.2 Cell-local repetition

Within each cell, how many molecules reuse the same genomic event structures?

## 14.3 Cross-cell repetition

Across cells, how many event structures recur?

## 14.4 Alternative-mapping multiplicity

How many molecules require multiple genomic placements to preserve replay fidelity?

## 14.5 Residual sequence burden

What fraction of molecules require storing any sequence witnesses at all?

If >90% of ordinary quantification can be replayed from coordinate/event evidence, that is encouraging.

---

# 15. Gate 0 success criteria

Proceed to a real compressed representation only if:

### Fidelity

Replay remains near-oracle on:

- baseline annotation;
- alternative annotation;
- controlled novel features.

### Storage potential

A plausible compressed molecular archive is projected to be substantially smaller than raw sequence data.

Initial target:

\[
\ge 5\times
\]

smaller than compressed FASTQ/BAM for the quantification-relevant representation.

A stronger result:

\[
10\times-100\times
\]

would make the project especially compelling.

### Query potential

Replaying annotation \(A\) from the archive should avoid expensive read mapping.

Even the unoptimized prototype should show that replay work scales with molecules/evidence rather than raw reads.

### Capability

At least one important class of annotation update must be accurately supported that cannot be recovered from an ordinary count matrix.

---

# 16. Gate 0 negative controls

## Matrix-only baseline

Attempt the same annotation changes from the original count matrix.

This should fail for annotation-dependent changes.

This demonstrates what information the archive preserves.

## Malva-like sequence pseudocount baseline

Where practical, compare sequence-derived feature counts from k-mer search with molecule-aware replay.

The objective is not to beat Malva at arbitrary sequence search.

The question is whether molecule-aware evidence gives substantially more faithful quantification.

## Full BAM replay

Use BAM/CRAM as a high-information baseline.

If the proposed archive is no more efficient than BAM/CRAM and offers no stronger semantics, novelty is weak.

---

# 17. Accuracy metrics

Do not rely on Pearson/Spearman correlation alone.

## Molecule-level metrics

- exact candidate-set agreement;
- best-feature agreement;
- molecule recall;
- molecule false-positive rate;
- equivalence-class Jaccard;
- disagreement stratified by ambiguity class.

## Cell × gene count metrics

- total count L1 error;
- maximum absolute count error;
- per-cell count error;
- per-gene count error;
- fraction of identical entries;
- false-positive count mass;
- false-negative count mass.

## Annotation-sensitive subset

Define:

\[
S_{12}
=
\{m:
Q(D,A_1)\neq Q(D,A_2)\}.
\]

Report all key metrics on \(S_{12}\) separately.

This subset matters more than overall correlation.

## Novel-feature metrics

For each feature added only in \(A_3\):

- supporting-molecule precision;
- supporting-molecule recall;
- count error;
- number of recovered cells.

---

# 18. First reference implementation

The first reference implementation should prioritize correctness.

Suggested components:

## Evidence extraction

Rust preferred if integrating with existing lab tools.

Input:

- FASTQ;
- cell-barcode whitelist;
- reference genome.

Output:

- molecule records in a simple Parquet/Arrow or binary table.

## Replay engine

Can initially be Python or Rust.

Input:

- evidence archive;
- annotation.

Output:

- molecule × feature compatibility;
- cell × feature count matrix.

## Evaluation harness

Python is sufficient.

Must compare fresh vs replay at:

- molecule level;
- cell level;
- feature level.

---

# 19. UMI grouping strategy

This is a critical design issue.

Annotation independence is weakened if UMI grouping depends strongly on gene assignment.

Therefore investigate at least two grouping models.

## G1 — conventional annotation-dependent grouping oracle

Use the trusted pipeline's existing grouping.

This provides the target result but cannot itself define the archive.

## G2 — annotation-independent grouping

Group reads by:

- cell barcode;
- UMI;
- genomic locality / alignment compatibility;
- strand;
- possibly splice path.

The goal is to infer molecule identity without assigning the molecule to a gene.

Gate -1 should explicitly compare G2 with G1.

This may become one of the main algorithmic contributions.

If annotation-independent molecular grouping is not sufficiently accurate, that is a serious project risk.

---

# 20. Possible molecule-key formulation

A candidate annotation-independent molecule key might combine:

\[
(c,u,L,o)
\]

where:

- \(c\): cell barcode;
- \(u\): UMI;
- \(L\): genomic locus or compatible genomic interval set;
- \(o\): orientation.

Reads with identical \(c,u\) but incompatible genomic placements remain separate molecule hypotheses.

The project should preserve uncertainty rather than forcing an early merge.

---

# 21. Evidence sufficiency hierarchy

Define increasingly rich archive classes.

## E0 — coordinate-only

Store:

- cell;
- UMI;
- genomic start/end;
- strand.

## E1 — splice-aware

Add:

- aligned blocks;
- splice junctions.

## E2 — edit-aware

Add:

- mismatch/indel summaries;
- variant witnesses.

## E3 — residual-sequence-aware

Add:

- sequence sketches/witnesses for unexplained sequence.

## E4 — near-lossless replay

Store any additional evidence required to reproduce fresh mapping decisions.

Gate -1 should determine the smallest level that achieves required fidelity.

This gives a clean experimental curve:

\[
\text{archive richness}
\quad\text{vs}\quad
\text{replay accuracy}.
\]

---

# 22. Potential theory problem

If the empirical concept passes, formulate a family of annotation queries.

Let \(\mathcal A_k\) be annotations expressible using features defined by:

- genomic intervals;
- at most \(k\) splice junctions;
- orientation;
- optional finite sequence witnesses.

Ask:

> What molecular evidence is sufficient to reproduce compatibility with every \(A\in\mathcal A_k\)?

Possible results may characterize:

- sufficient event summaries;
- lower bounds on information required;
- impossibility results for coordinate-only representations;
- exactness conditions.

Do not force a theorem before the empirical abstraction is stable.

---

# 23. Potential data-structure problem

Once the molecule evidence semantics are fixed, the main indexing problem becomes:

> Store a large set of molecule records and support rapid annotation-to-molecule compatibility queries.

Possible primitive indexes:

- interval index;
- splice-junction index;
- path/event dictionary;
- sequence witness index;
- molecule posting lists;
- cell posting lists.

Repeated genomic event patterns may allow strong interning/factorization.

Important hierarchy:

\[
\text{event}
\rightarrow
\text{molecules}
\rightarrow
\text{cells}.
\]

Do not collapse directly to cells if doing so destroys UMI-level count semantics.

---

# 24. Large-scale architecture — deferred

Only after Gates -1 and 0.

Potential long-term architecture:

## Reference/event dictionary

Global interned dictionary of:

- genomic blocks;
- splice junctions;
- variants;
- residual sequence witnesses.

## Molecule table

Each molecule references a compact event signature.

## Cell table

Each cell references molecule IDs/ranges.

## Annotation compiler

Transforms a new annotation into constraints over event signatures.

## Query engine

Maps annotation features to compatible molecule cohorts.

## Quantifier

Runs desired molecule-resolution/counting policy.

This architecture may ultimately support both:

- rapid re-quantification;
- arbitrary sequence/event queries.

---

# 25. Malva comparison plan

Malva is relevant evidence and prior art, not merely a competitor.

The final paper should explicitly compare capabilities.

## Tasks Malva is expected to excel at

- arbitrary k-mer/sequence queries;
- pathogen search;
- mutation search;
- massive-scale indexing.

We should not attempt to beat Malva on these unless the architecture unexpectedly provides a large advantage.

## Tasks the proposed archive should excel at

- UMI/molecule counting;
- re-quantification under changed annotation;
- ambiguity-preserving molecule compatibility;
- annotation-sensitive count replay;
- re-use of one molecule across multiple feature definitions.

A successful paper should demonstrate qualitative differences, not only runtime differences.

---

# 26. Initial datasets

Use small/moderate datasets first.

Recommended:

## Dataset D0 — development subset

- 500–1,000 cells;
- standard 10x 3' PBMC or similarly simple dataset;
- small enough for complete repeated processing in minutes.

Purpose:

- implementation debugging;
- annotation replay;
- UMI-grouping tests.

## Dataset D1 — first real Gate -1

- ~5k cells;
- standard 10x 3' dataset;
- good existing annotations and known pipeline behavior.

## Dataset D2 — annotation-sensitive confirmation

- 5k–10k cells;
- choose a dataset where annotation differences meaningfully affect genes/transcripts.

## Dataset D3 — single-nucleus or 5' confirmation

Only after success.

This tests whether the representation generalizes beyond 3' Chromium.

Do not start with 100k+ cells.

---

# 27. Required experiment matrix

For D0 and D1:

| Archive | A1 replay | A2 replay | A3 replay |
|---|---:|---:|---:|
| E0 coordinate-only | yes | yes | yes |
| E1 splice-aware | yes | yes | yes |
| E2 edit-aware | yes | yes | yes |
| E3 residual sequence | yes | yes | yes |
| fresh raw-read oracle | reference | reference | reference |

For each, report:

- molecular fidelity;
- count fidelity;
- storage;
- build time;
- replay time.

This immediately reveals how much evidence is actually needed.

---

# 28. Fast milestone schedule

The agent should work in short, gated phases.

## Milestone M0 — 1–2 days: oracle and annotation panel

Deliver:

- frozen fresh-processing pipeline;
- D0 dataset;
- A1/A2/A3 annotations;
- molecule-level comparison tooling.

No new archive yet.

## Milestone M1 — 2–4 days: deliberately overcomplete archive

Deliver:

- annotation-independent molecule grouping prototype;
- near-lossless evidence record;
- A1 replay.

Decision:

Can A1 be reproduced?

If no, stop and diagnose.

## Milestone M2 — 2–4 days: cross-annotation replay

Deliver:

- A2 and A3 replay;
- annotation-sensitive subset metrics.

Decision:

Does annotation independence actually work?

If no, stop.

## Milestone M3 — 2–3 days: evidence ablation

Test E0/E1/E2/E3.

Decision:

What is the minimal useful evidence class?

## Milestone M4 — 2–3 days: compression opportunity audit

Measure:

- evidence volume;
- repetition;
- theoretical storage opportunity;
- replay runtime.

Decision:

Is there enough headroom to justify a new data structure?

Only after M4 should optimized indexing begin.

---

# 29. Hard continuation gates

Proceed to full development only if all of the following hold.

## Gate A — capability

One archive supports at least A1/A2/A3 without rebuilding.

## Gate B — fidelity

For ordinary gene-level quantification, replay is nearly indistinguishable from fresh processing.

Target:

- >99% agreement on high-confidence molecule assignments;
- very small cell × gene count L1 error.

## Gate C — annotation-sensitive fidelity

On molecules whose assignment changes across A1/A2, archive replay captures the fresh-processing change accurately.

## Gate D — genuinely new features

A3 features introduced after archive construction can recover supporting molecules with high precision/recall.

## Gate E — storage headroom

The minimal sufficient archive has substantial projected storage advantage over raw read/alignment representations.

## Gate F — replay speed headroom

Replay avoids read mapping and is projected to be dramatically faster than fresh processing.

If Gates A–D pass but E/F do not, the project may still be scientifically interesting but needs a different novelty argument.

---

# 30. Stop conditions

Stop as a primary target if:

- faithful replay requires near-complete read retention;
- annotation-independent UMI grouping is fundamentally unreliable;
- annotation-sensitive molecules cannot be reproduced accurately;
- most annotation changes require sequence information absent from the compact representation;
- replay is only modestly faster than reading/reprocessing BAM;
- Malva-style sequence search plus conventional quantification already provides essentially the same capability.

---

# 31. Potential framing if successful

Possible title:

**An annotation-independent molecular evidence representation for single-cell RNA-seq**

Core claim:

> Single-cell count matrices discard sequence and molecule information needed to reinterpret an experiment under future annotations, whereas raw reads retain substantial redundant information and are expensive to revisit. We introduce an annotation-independent molecule-level representation that preserves the evidence necessary for future feature assignment while avoiding commitment to a particular transcript annotation. A single archive can be re-quantified under multiple gene/transcript definitions without remapping the original reads.

Potential contributions:

1. **New abstraction:** annotation-independent molecule evidence.
2. **Correctness/fidelity characterization:** which future annotations can be replayed.
3. **Algorithm:** annotation-independent UMI grouping and evidence compilation.
4. **Data structure:** compact event/molecule representation and rapid query.
5. **Application:** rapid re-quantification under updated annotations.
6. **Scale:** eventual application to very large single-cell collections.

---

# 32. Relation to Malva

Malva provides strong evidence that preserving sequence-level information beyond the count matrix is valuable.

The proposed project is complementary but must remain clearly distinct.

Malva's fundamental abstraction is sequence-indexed cell evidence.

The proposed abstraction is:

\[
\boxed{\text{molecule-indexed biological evidence}}
\]

designed specifically to preserve UMI-aware quantification semantics.

The final novelty claim should be:

> Not merely “which cells contain this sequence?” but “which physical captured molecules support this feature under the annotation I am asking about now?”

That is the intended gap.

---

# 33. Concrete first tasks for the implementation agent

Perform these in order.

## Task 1 — create repository structure

```text
docs/
experiments/
src/
scripts/
tests/
results/
annotations/
```

Add this plan as:

```text
docs/project_plan.md
```

Create:

```text
docs/implementation_status.md
docs/decision_log.md
```

## Task 2 — freeze D0

Select a small 500–1,000-cell Chromium dataset.

Record:

- FASTQ paths;
- genome;
- chemistry;
- whitelist;
- expected cell set;
- existing fresh-processing commands.

## Task 3 — define A1, A2, A3

Create exact annotation manifests.

A3 must include at least five controlled changes representing different classes:

- new exon;
- new splice junction;
- altered transcript end;
- split/merged gene;
- new transcript.

## Task 4 — build fresh oracle outputs

Run the trusted pipeline independently under A1/A2/A3.

Export molecule-level assignments where possible, not only count matrices.

## Task 5 — implement annotation-independent molecule grouping

Do not use transcript/gene IDs.

Use only:

- cell barcode;
- UMI;
- genomic evidence;
- strand;
- local compatibility.

Produce diagnostics against the fresh oracle.

## Task 6 — implement overcomplete E4 record

Retain everything plausibly needed.

Use a simple format.

Do not compress.

## Task 7 — implement A1 replay

Compare molecule by molecule.

Classify every disagreement.

Do not move to A2 until the disagreement taxonomy is understood.

## Task 8 — implement A2 replay

Focus evaluation on annotation-sensitive molecules.

## Task 9 — implement A3 replay

Test genuinely unseen features.

## Task 10 — evidence ablation

Generate E0/E1/E2/E3 from E4.

Measure replay degradation.

## Task 11 — Gate -1 report

Create:

```text
experiments/reports/gate_minus1.md
```

Answer:

1. Can annotation-independent evidence reproduce fresh quantification?
2. Which evidence fields are essential?
3. Can genuinely new features be replayed?
4. Does annotation-independent UMI grouping work?
5. What failures remain?

End with:

```text
GO
```

or

```text
STOP
```

## Task 12 — only after GO: compression audit

Measure storage and repetition structure.

No optimized index before this point.

---

# 34. Final philosophy

The last several feasibility attempts failed because they sought additional gains from already highly optimized methods.

This project should be treated differently.

The initial question is not:

> Can we make quantification 20% faster?

It is:

> Can we preserve a compact intermediate biological representation that lets us ask future annotation-dependent questions that a count matrix cannot answer, without paying the cost of returning to raw reads?

That is the capability to prove first.

If the capability exists, the data-structure and scaling work can follow.

If it does not, the project should be stopped quickly before substantial engineering effort is invested.
