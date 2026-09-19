# Gravlax post-v1 red-team review and eight-week roadmap

Date: 2026-08-30  
Frozen base: `gravlax-v1-8-30-2026`  
Review code: `b5a791e681e7f1f93d46459745eedba59ac841d6`  
Review paper: `1a1ef0bf3501ee57c009e37661c773eff88ff046`

## Executive verdict

The project has a publishable central contribution: a compact, queryable molecular evidence
abstraction that makes compatible-annotation recomputation interactive. No finding invalidates
that contribution. The strongest paper is not “a universally minimal sufficient statistic” and
not “exact STARsolo without reads.” It is:

> Conditional on fixed placements and barcode correction, Gravlax retains the evidence consumed
> by named annotation-assignment policies, applies an explicitly measured molecule reduction,
> and turns repeated annotation analysis into a small indexed replay problem.

The v1 artifact was technically strong but its broadest words outran its evidence. The post-v1
branch has already corrected the exactness taxonomy, removed the unproved minimality claim,
hardened malformed-input handling, replaced a stale runtime denominator, and established a much
fairer D0 storage baseline. The next paper-changing priority is peak memory, not another entropy
codec: a 111 MB D0 index currently reaches about 6.95 GB RSS during replay.

## Findings by severity

### Claim-blocking, now corrected

1. **Three different comparisons were called exact.** Archive- and BAM-sourced Gravlax replay is
   byte-identical. The two-representative reduction has incremental error. Fresh STARsolo differs
   by 0.22–0.45% of UMI mass. These are now separate claims.
2. **“Minimal” had no lower-bound proof.** Necessary-field counterexamples support conditional
   sufficiency for a named consumer, not global minimality. The headline is corrected.
3. **The 150–370× speed headline was not traceable to the current D0 oracle command.** The frozen
   oracle log itself spans about 192 s, while the manuscript used a roughly 16-minute D0
   denominator. Five current Gene-only no-BAM runs give a defensible D0 median of 83.68 s versus
   2.42 s, or 34.6×. The stale headline is removed.

### Submission blockers

1. **Storage capability matching is incomplete.** The D0 reduced CRAM reproduces the archive and
   matrices byte-for-byte and is still 6.83× larger, but it retains raw correction evidence,
   names, and every alignment instance. Call it ingest-equivalent, not a lower bound.
2. **Runtime replication is incomplete.** Five warm D0 Gene-only runs exist, but arms were
   sequential, a controlled cold-cache observation is missing, and D1/D2' are unmeasured.
3. **Generalization is all 3′ human.** One public mouse 5′ experiment is the minimum credible
   chemistry/species extension.
4. **Peak RSS weakens the “interactive” story.** Replay is fast but eager archive materialization
   expands D0 to roughly 6.95 GB and was previously measured at about 18.5 GB on D1.

### Important engineering risks

- The seekable reader previously trusted future versions, section-directory extents, rANS
  residual lengths, and several integer conversions. The post-v1 code fails closed on these and
  adds an end-to-end tagged-BAM CLI golden test; fuzzing and a permanent corruption corpus remain.
- `replay` and `eval` are reserved empty crates. Public documentation now says so; do not present
  them as implemented libraries.
- The CLI and format documentation had stale examples and described nonexistent fields. Corrected
  post-v1, but documentation should become a release gate.
- GTF parsing now rejects truncated, reversed, zero, and invalid-strand exons. Coordinate-system,
  contig-alias, and liftover tests remain sparse.

### Speculative, not blockers

- Cross-sample dictionary sharing may help atlas storage, but the measured federated open/query
  path is already fast and no redundancy audit yet justifies a container redesign.
- Sparse sequence evidence could enable allele/fusion/edit queries, but it should remain an
  opt-in E2/E3 layer tied to a demonstrated biological question.
- A formal information-theoretic minimum over all future annotations is unlikely to be useful;
  a consumer-conditional theorem plus counterexamples is sufficient and more honest.

## Prior-art positioning

- CRAM is the general, reference-based alignment compressor; Gravlax must beat a reduced CRAM at
  a declared capability vector, not just an all-tag CRAM. CRAM 3.1 reports substantial gains over
  BAM but preserves an alignment-file abstraction ([Bonfield et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8896640/)).
- Boiler is the closest compression analogy: it intentionally discards per-read information in
  favor of coverage vectors and empirical distributions, trading exact read identity for common
  bulk-RNA analyses ([Pritt and Langmead](https://pmc.ncbi.nlm.nih.gov/articles/PMC5027496/)).
  Gravlax's differentiator is barcode/UMI molecule identity and annotation-policy replay.
- Monorail/recount3 demonstrates atlas-scale annotation-agnostic coverage and junction summaries,
  but exposes study-level bulk summaries rather than per-cell molecular evidence
  ([Wilks et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8628444/)).
- BUS and RAD are reduced pseudoalignment/intermediate representations. They validate the general
  “store reusable evidence” direction, but their target-set semantics and transcriptome-index
  dependency differ from genome-coordinate annotation replay
  ([BUS](https://academic.oup.com/bioinformatics/article/35/21/4472/5487510),
  [alevin-fry/RAD](https://pmc.ncbi.nlm.nih.gov/articles/PMC8933848/)).
- sc-SPLASH and Malva target reference-free sequence discovery, a richer and different query
  family. They are complements or capability-rich baselines, not substitutes for compatible-GTF
  replay ([sc-SPLASH](https://www.nature.com/articles/s41587-026-03084-6),
  [Malva](https://www.nature.com/articles/s41586-026-10975-w)).

## Ranked size opportunities

1. **Called-cell and ambient tiers.** Separate called-cell primary evidence from background
   barcodes while preserving a compact ambient rescue tier. Measure before designing the format;
   this has the largest plausible dataset-dependent gain.
2. **Post-correction matched representation.** Encode the exact corrected cell/UMI and placement
   relations needed for current replay in BAM/CRAM or a neutral columnar table. This makes the
   storage claim reviewer-proof even if the ratio narrows.
3. **Adaptive evidence richness.** Keep E1 geometry for ordinary molecules and promote only
   edit/soft-clip/ambiguous records to E2/E3. Require a query that consumes the promoted data.
4. **Depth curves and saturation-aware coding.** Quantify bytes/read and bytes/molecule at fixed
   FASTQ prefixes before tuning class, cell, or pattern codes.
5. **Cross-sample immutable dictionaries.** Price shared genome/shape/pattern dictionaries in an
   atlas without making individual archives dependent on an unavailable global object.
6. **Index ablation.** Remove each index in turn and report bytes saved and queries lost. Existing
   D0 index sections are small enough that this is likely a clarity result, not the main win.

## Ranked speed and memory opportunities

1. **Streaming replay.** Replace eager `read_archive` materialization with chunk iteration plus
   bounded reducers. Gate: byte-identical outputs, D0 RSS below 2 GB, D1 below 6 GB, and wall time
   no worse than 1.25× current.
2. **Out-of-core cell sharding.** Spill `(cell,gene,UMI)` reducers by a stable cell partition and
   merge deterministic shards. This is the fallback if streaming alone does not bound aggregation.
3. **Compiled annotations.** Serialize contig-normalized exon intervals and feature metadata.
   Measure cold parse/compile separately from replay; key by GTF and genome digest.
4. **Batched query execution.** Share chunk decoding across region/junction/APA/set requests and
   expose a query-plan summary so atlas workloads do not repeat decompression.
5. **EM allocation/profile pass.** Reuse target buffers, compress class-to-cell traversal, and
   evaluate hierarchical sharing only after a flamegraph identifies the actual hot path.
6. **Direct CRAM input.** Either support CRAM in the BAM-sourced replay arm or report conversion
   time explicitly; do not silently benchmark a preconverted file.

## Ranked capability opportunities

1. Stable molecule export with archive/annotation/provenance identifiers.
2. Cell and named-group filters on every replay and query path.
3. Junction enumeration plus PSI/set queries with explicit denominators.
4. Hierarchical EM with shrinkage across related cell groups and held-out likelihood reporting.
5. Mouse 5′ ingestion/replay, including 10-base UMI geometry and strand validation.
6. Coordinate liftover with rejected/ambiguous interval accounting.
7. Sparse E2/E3 allele, editing, fusion, or soft-clip queries tied to orthogonal truth.

## Experiments, in execution order

### 1. Fair storage and runtime baselines

- **Hypothesis:** `.aie` remains at least 2× smaller and 10× faster than the smallest comparator
  retaining the named replay capability.
- **Implementation:** use the accepted CR/CY/UR/NH reduced transform, add a post-correction neutral
  representation, direct CRAM or explicit conversion, Gene-only no-BAM STARsolo, and randomized
  warm/cold timing.
- **Datasets/controls:** D0, D1, D2'; identical FASTQ, genome, GTF, whitelist, cell list, threads;
  byte-identical within-pipeline regressions and nonzero fresh-STARsolo deviation reported apart.
- **Metrics/gates:** bytes, bytes/molecule, wall/CPU/RSS/I/O; five warm and one targeted cold run;
  replace all storage/runtime headlines with these values.
- **Failure meaning:** near-sized CRAM shifts the paper from compression to queryability; a small
  speedup shifts it to avoided repeated alignment and capability.
- **Compute/value:** 2–4 host-days; critical.

### 2. Controlled sequencing-depth scaling

- **Hypothesis:** molecule collapse makes `.aie` grow sublinearly after library saturation while
  BAM/CRAM remains closer to read-linear.
- **Implementation:** deterministic read-name-hash prefixes at 10/25/50/75/100% of D1; align each
  prefix once, build all artifact tiers, fit size and wall models against reads and molecules.
- **Controls/metrics:** identical cell whitelist/policy; three independent hash seeds; confidence
  intervals for slope, bits/read, bytes/molecule, duplication, and called-cell saturation.
- **Failure meaning:** linear scaling means current compression is structural but not increasingly
  advantageous with depth.
- **Compute/value:** 1–2 host-days and about one D1 BAM of temporary space; high.

### 3. Called-cell versus full-background profiles

- **Hypothesis:** most D2' archive overhead comes from ambient barcodes that can be tiered without
  damaging called-cell replay or ambient rescue.
- **Implementation:** build called-only, full, and called-plus-ambient-sketch archives from the
  same extracted rows; freeze the cell caller; support later promotion from the ambient tier.
- **Controls/metrics:** called-cell matrix identity, recovered cells after alternate filtering,
  bytes by section, molecules retained, false promotion rate.
- **Failure meaning:** a small size delta rules out a format split and favors simple filtering.
- **Compute/value:** 1 host-day after a 2–3 day prototype; high if D2' dominates future atlases.

### 4. Streaming replay

- **Hypothesis:** most RSS is eager materialization rather than irreducible UMI state.
- **Implementation:** expose chunk iterators, stream annotation assignment, bound per-cell shards,
  and perform a deterministic external merge only when necessary.
- **Controls/metrics:** byte-identical outputs for D0/D1/D2'; current versus streaming wall, RSS,
  decoded bytes, spill bytes; gates are <2/<6 GB D0/D1 and ≤1.25× wall.
- **Failure meaning:** if aggregation state dominates, advance cell sharding; if decode dominates,
  retain eager mode for small archives.
- **Compute/value:** 3–5 engineering days plus 1 benchmark day; critical.

### 5. Mouse 5′ generalization

- **Hypothesis:** the evidence abstraction is chemistry/species independent once barcode/UMI
  geometry and strandedness are parameterized.
- **Dataset:** 10x “5k Mouse PBMCs,” Universal 5′ GEM-X, GRCm39 with current GENCODE/Ensembl mouse
  annotation; retain the older C57BL/6 5′ dataset as a compatibility fallback
  ([dataset](https://www.10xgenomics.com/jp/datasets/5k_Mouse_PBMCs_5p_gem-x)).
- **Implementation/controls:** chemistry manifest drives CB/UMI positions; fresh STARsolo/Cell
  Ranger-compatible oracle; no hard-coded 12-base UMI or human contigs; run storage/replay gates.
- **Metrics/gates:** ≤1% Gene UMI-mass deviation, archive/BAM serialization identity, successful
  velocity strata where supported, section size and runtime.
- **Failure meaning:** chemistry-specific assumptions become explicit format/tool limitations.
- **Compute/value:** 2–3 host-days plus download/index time; critical for generality.

### 6. Downstream clustering, marker, and velocity fidelity

- **Hypothesis:** the reported UMI-mass deviations do not alter biological conclusions.
- **Implementation:** fixed Scanpy pipeline over oracle/replay matrices; transfer a single PCA/NN
  parameter set; evaluate Velocyto components separately.
- **Datasets/controls:** D0/D1/D2'; same called cells, normalization, seeds, and feature filters;
  include a deliberately perturbed positive control.
- **Metrics/gates:** ARI/NMI, neighbor Jaccard, marker top-k overlap/effect correlation, pseudobulk
  correlation, velocity cosine/transition agreement; pre-register acceptable deltas.
- **Failure meaning:** small mass error may be structured; report affected cell types/features and
  adjust the scientific claim rather than averaging it away.
- **Compute/value:** 1–2 days; high paper value.

### 7. Complete-denominator discovery with orthogonal truth

- **Hypothesis:** discovery gains persist when every eligible event is counted, not only events
  emitted by one method.
- **Implementation:** union candidate universe from Gravlax, STAR junctions, PolyASite, and a
  sequence-aware caller; score emitted and non-emitted events under one denominator.
- **Controls/metrics:** synthetic spike-ins plus independent junction/3′ catalogues; precision,
  recall, PR-AUC, calibration, support-stratified error, and missing-denominator accounting.
- **Failure meaning:** high precision/low recall narrows the claim to targeted rediscovery.
- **Compute/value:** 3–5 days; high, especially for the discovery story.

### 8. Less biased EM and hierarchical sharing

- **Hypothesis:** cross-cell sharing improves sparse cells without leakage from the evaluation
  truth construction.
- **Implementation:** simulate from held-out genes/classes and mask unique evidence by a mechanism
  independent of the fitted model; compare per-cell, pooled, group-hierarchical, and oracle priors.
- **Controls/metrics:** held-out log likelihood, top-1 recovery, calibration, abundance error, rare
  cell-type stratification, negative-control unrelated grouping.
- **Failure meaning:** pooled improvement caused by truth leakage or dominant-cell bias removes the
  biological sharing claim but not EM functionality.
- **Compute/value:** 3–5 days; high.

### 9. Sparse E2/E3 biological prototype

- **Hypothesis:** promoting only anomalous reads enables one sequence-dependent question at a
  modest storage premium.
- **Implementation:** retain edit positions/bases and residual soft clips only for a predeclared
  target such as expressed SNVs or RNA editing; build an indexed allele query.
- **Controls/metrics:** matched BAM plus orthogonal variant/edit catalogue; precision/recall,
  allele-count concordance, bytes promoted, query latency; storage premium gate ≤25%.
- **Failure meaning:** poor truth or >25% overhead keeps sequence evidence out of v2.
- **Compute/value:** 1 week; medium, do only after the base paper gates.

### 10. Formal conditional sufficiency and counterexamples

- **Hypothesis:** equality can be proved for named consumers over fixed placements, while dropping
  each retained field yields a concrete annotation/policy counterexample.
- **Implementation:** formalize consumer input projections; add property tests over synthetic
  placements/GTFs and minimal counterexamples for strand, blocks, junctions, NH, cell, and UMI.
- **Controls/metrics:** proof statement matches code; every necessary-field counterexample fails
  after field removal and passes with the field restored.
- **Failure meaning:** narrow the consumer family until code and statement agree; never restore the
  global “any annotation/minimal” wording.
- **Compute/value:** 2–3 days; high conceptual value, low compute.

## Do not pursue now

- Another generic entropy coder: fitted rANS moved total size by at most about 1%.
- A v2 container rewrite without a measured ≥2× size win or a memory/query capability win.
- GPU acceleration before ingest/replay/EM profiles identify a suitable dominant kernel.
- Sequence-rich evidence by default; it undermines the compact abstraction and duplicates FASTQ.
- A universal minimality theorem or “any annotation” claim.
- Atlas-wide shared dictionaries before measuring redundancy and failure/recovery semantics.
- New public CLI surface in the empty `replay`/`eval` crates before real consumers exist.

## Decision-complete eight-week roadmap

| week | work | exit gate |
|---|---|---|
| 1 | snapshot, claim audit, corruption guards, CLI golden, D0 reduced BAM/CRAM and counts-only pilots | complete: tags verified; 86 tests; D0 6.83× reduced-CRAM ratio; 34.6× warm Gene-only ratio |
| 2 | randomized/cold D0 benchmark; post-correction comparator; D1/D2' reduced storage runs | exact capability table; five warm + one targeted cold; no stale denominator |
| 3 | controlled D1 depth series and called/background section accounting | section-level curves with three seeds; decide whether tiering earns a prototype |
| 4 | streaming reader and deterministic streaming replay | byte identity; D0 <2 GB, D1 <6 GB RSS; ≤1.25× wall or documented fallback |
| 5 | mouse 5′ ingest/replay plus chemistry manifest | full oracle/archive gates; no human/12-base assumptions |
| 6 | clustering/marker/velocity fidelity and complete-denominator discovery | biological metrics and failure strata frozen; weak stories removed |
| 7 | unbiased/hierarchical EM evaluation; conditional-sufficiency property tests | held-out likelihood/recovery report; theorem language matches code |
| 8 | only if earned: sparse E2/E3 prototype; otherwise figures, manuscript, release hardening | v2 change requires ≥2× or unique capability; full tests/docs/reproduction pack clean |

At the end of each week, update `experiments/manifests/post-v1-gates.yaml`. A failed gate changes
the paper claim or stops the branch; it does not silently change the metric.
