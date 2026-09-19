# Gravlax paper reproducibility record

This repository archives the scripts, decisions, compact results, and provenance behind the
first complete Gravlax manuscript and `.aie` implementation. It intentionally does **not** copy
the implementation, manuscript source, raw reads, alignments, reference indexes, or molecular
archives.

The first complete snapshot (2026-08-30), which the documents below call `gravlax-v1-8-30-2026`,
is identified by its source revisions:

- Gravlax code: [`d9257ca2c81569aa14dfeb2899beb4b90ebb5a9d`](https://github.com/COMBINE-lab/gravlax/commit/d9257ca2c81569aa14dfeb2899beb4b90ebb5a9d)
- Paper: [`67da17775aed087a16ab9dbd1aea91e3eca7160a`](https://github.com/COMBINE-lab/moleceular-evidence-store-paper/commit/67da17775aed087a16ab9dbd1aea91e3eca7160a)

The bioRxiv preprint (submitted 2026-09-18) reports those results together with the post-v1 work
below and the 2026-09-11 follow-up experiments. This repository is published as a single snapshot
at the preprint state; its authoritative source revisions are:

- Gravlax code: release [`v0.2.3`](https://github.com/COMBINE-lab/gravlax/releases/tag/v0.2.3)
  (`75b8d6c01064ba92af295543d50230429774e170`); the follow-up experiments were first run with
  release `v0.2.2` and the cohort discovery was rerun with `v0.2.3`.
- Paper: [`b41a35462a80f0f92617b426a8611dbb81ab639b`](https://github.com/COMBINE-lab/moleceular-evidence-store-paper/commit/b41a35462a80f0f92617b426a8611dbb81ab639b)
- Reproducibility scripts: this repository at `gravlax-biorxiv-2026-09`

## What is here

- `PROVENANCE.md` — the self-contained experiment and artifact catalogue.
- `scripts/` — the exact experiment drivers used on the analysis host.
- `docs/decisions/` — thresholds frozen before the corresponding measurements.
- `docs/decision_log.md` and `memos/` — chronological decisions, failures, negative results,
  and gate verdicts.
- `results/` — compact machine-readable outputs used to support the reported conclusions.
- `logs/` — command transcripts and timing output; these are small and retained as evidence.
- `experiments/manifests/gates.yaml` — the original machine, dataset, and gate manifest.
- `environment/setup.sh` — the exact host-local environment setup used for the snapshot.
- `manifests/` — immutable input identities and inventories of excluded large artifacts.

The host-local setup script preserves its original toolchain layout below the project root. It is a
record of the executed environment, not a portable installer. See `manifests/software.tsv` and
`PROVENANCE.md` when recreating the environment elsewhere.

## Project root variable

Every script locates the analysis project through one environment variable,
`GRAVLAX_PROJECT_ROOT`, the directory that holds the excluded `data/`, `annotations/`, `ref/`,
and `runs/` trees described in the manifests. Shell drivers fail immediately if it is unset, and
Python drivers read it from `os.environ`. Two optional variables select tools: `GRAVLAX_PYTHON`
(the interpreter used for the analysis environment, default `python3`) and `GRAVLAX_ZSTD` (the
`zstd` binary used by script `211`, default `zstd` on `PATH`).

```bash
export GRAVLAX_PROJECT_ROOT=/path/to/analysis-project
export GRAVLAX_PYTHON=/path/to/env/bin/python   # optional
```

Records under `results/`, `logs/`, `PROVENANCE.md`, and `memos/` write the same root as the
literal `${GRAVLAX_PROJECT_ROOT}` wherever a command transcript or artifact identity contains an
absolute path. Relative paths in those records are relative to that root, which is what `$P`
denotes in `PROVENANCE.md`. A few records name locations outside the project root with the same
placeholder style: `${GRAVLAX_SRC}` (the Gravlax source clone), `${GRAVLAX_PAPER_SRC}` (the
manuscript clone), `${GRAVLAX_TOOLS}` (the Node toolchain used for the docs site),
`${CONDA_ROOT}`, `${MALVA_ROOT}` (the Malva 1.0.0 checkout used as a storage baseline),
`${STAR_SOURCE}` (the STAR source tree read for rule fidelity), and `${WHITELIST_SOURCE}` (the
directory holding the 10x `3M-february-2018.txt` whitelist). The analysis machine is written as
`analysis-host`; its hardware is described in `PROVENANCE.md`.

## Result-to-command map

Run commands from a project root containing the excluded `data/`, `annotations/`, `ref/`, and
`runs/` trees described in the manifests, with `GRAVLAX_PROJECT_ROOT` exported and after adapting
or sourcing `environment/setup.sh`.

| Result family | Primary drivers |
|---|---|
| STAR reference indexes and fresh oracles | `scripts/05_build_anno_indices.sh`, `scripts/10_oracle_starsolo.sh` |
| Annotation-signal and alignment gates | `scripts/20_gatec_align_configs.sh`, `scripts/30_gateb_compare.py`, `scripts/31_delta_capture.py` |
| Annotation-free ingest and storage baselines | `scripts/40_ingest_align.sh`, `scripts/50_gatee_baselines.sh` |
| Withheld-annotation and discovery tests | `scripts/70_withhold_annotation.py`, `scripts/80_apa_gate.py`, `scripts/81_disc_apa_v2.py`, `scripts/82_disc3_classify.py`, `scripts/92_discr_gate.py`, `scripts/98_prospective.py` |
| Independent D1/D2/D2' campaigns | `scripts/90_d1_reverify.sh`, `scripts/96_d2_campaign.sh`, `scripts/97_d2p_campaign.sh` |
| Velocity replay | `scripts/91_velo_gate.py`, `scripts/94_velo2_strata.py` |
| APA replication and figures | `scripts/93_apa_diff.py`, `scripts/100_apa_replicate.py`, `scripts/105_ptprc.py` |
| EM recovery and external comparison | `scripts/95_em0_gate.py`, `scripts/99_em_af_concordance.py`, `scripts/101_replication_pack.sh`, `scripts/104_em_de.py` |
| Reference extension and federation | `scripts/102_pool_extension.py`, `scripts/103_atlas_build.sh` |
| Function-matched molecule BAM/CRAM | `scripts/117_stream_molecule_cram_replay.sh`, `scripts/118_post_correction_storage.sh`, `scripts/119_native_aie_replay.sh` |
| Bounded versus eager replay | `scripts/121_streaming_replay_gate.sh`, `scripts/122_summarize_streaming_replay.py` |
| Packed, cell-sharded recovery EM | `scripts/133_compare_em_layers.py`, `scripts/134_packed_sharded_em_gate.sh`, `scripts/135_summarize_packed_sharded_em.py` |
| Mouse 5′ v2 generalization | `scripts/123_fetch_mouse_5p.sh`, `scripts/124_mouse_5p_gate.sh`, `scripts/125_summarize_mouse_5p.py` |
| Downstream clustering, markers, and velocity | frozen stress test: `scripts/126_downstream_fidelity.py`, `scripts/127_downstream_fidelity_gate.sh`, `scripts/128_downstream_cluster_diagnostics.py`; stable macro-partitions: `scripts/129_multiresolution_select.py`--`scripts/132_multiresolution_evaluate.sh` |
| Federated splice-event biology and FNBP1 confirmation | screen: scripts `179`--`181`; untouched confirmation: scripts `182`--`184`; paper figure: script `185` |
| D5 broad-cell FNBP1 localization | frozen label/query workflow: scripts `186`--`189`; compact result and memo: `results/post-v1-fnbp1-broad-cell-localization.*`, `memos/fnbp1-broad-cell-localization.md` |
| Multi-donor FNBP1 validation | prospective GSE234790 gate and protocol: `experiments/manifests/fnbp1-multidonor-sez-validation.json`, `experiments/specs/fnbp1-multidonor-sez-validation.md`; event-blind feasibility, verified acquisition, streamed runner, and registered summarizer: scripts `190`--`193`; compact STOP result and audit: `results/post-v1-fnbp1-multidonor-sez-validation.*`, `memos/fnbp1-multidonor-sez-validation.md` |
| SEZ transcript-end atlas | prospective protocol-native APA/ALE gate: `experiments/manifests/sez-transcript-end-atlas-v1.json` and `experiments/specs/sez-transcript-end-atlas-v1.md`; published-cluster feasibility and frozen donor/group design: script `196`, `experiments/designs/sez-transcript-end-atlas-v1.tsv`, and `results/post-v1-sez-transcript-end-feasibility/` |
| SEZ PolyASite mixture | method-amended protocol-aware 3′ analysis: scripts `197`, `199`--`203`; cross-fitted controls and NTRK2 terminal-architecture result: `results/post-v1-sez-polyasite-mixture*`, with interpretation in `memos/sez-polyasite-mixture.md` |
| Federated query workflow optimization | concurrent driver and streaming summarizer: scripts `180`--`181`; exact engine/full-workflow gates: `results/post-v1-federated-query-workflow-optimization.json`, `memos/federated-query-workflow-optimization.md` |
| Sparse federated cohort output | fail-closed dense-JSON reconstruction: script `205`; compact equivalence/size/timing record: `results/post-v1-federated-sparse-output-validation.json`, with interpretation in `memos/post-v1-federated-sparse-output.md` |
| Optional federated collection index | final benchmark driver and fail-closed reducer: scripts `206`--`207`; deterministic root, immutable extension, exact routed queries, and matched naive timings: `results/post-v1-collection-index-v2.json`, `memos/post-v1-collection-index-v2.md`; exact commands: `PROVENANCE.md` §8.25 |
| Adaptive evidence-frontier negative pilot | dedicated frozen gate: `experiments/specs/adaptive-evidence-frontier.md`; compact STOP result and formal audit: `results/post-v1-adaptive-evidence-frontier.json`, `memos/adaptive-evidence-frontier.md`; excluded artifact identities: `manifests/adaptive-evidence-frontier-artifacts.tsv`; exact commands: `PROVENANCE.md` §8.23 |
| Bounded-memory PolyASite cohort analysis | fail-closed summary: script `204`; exact-table/RSS/runtime gate: `results/post-v1-polyasite-bounded-memory.json`; implementation and limits: `memos/post-v1-polyasite-bounded-memory.md`; exact commands: `PROVENANCE.md` §8.24 |
| Sparse terminal-tail evidence | prospective pilot/promotion gates: `experiments/specs/2026-09-01-sparse-terminal-tail*.md`; compact three-way verdict: `results/post-v1-sparse-terminal-tail.json`; biological and representation audit: `memos/sparse-terminal-tail-evidence.md`; exact excluded identities: `manifests/sparse-terminal-tail-artifacts.tsv`; commands: `PROVENANCE.md` §8.26 |
| Archive root and local-shape routing | Gate A PASS: exact encoded-section commitments add 0.026835% and reduce collection-build source I/O 129.19-fold. Gate B confirmatory r3 PASS: all 101 root/chain routed-versus-fallback cases and 26 integrity predicates pass; the route premium is 2.289131% of source bytes, build wall is 1.900709× fallback at 607,132 KiB, and prospective-panel aggregate wall is 0.386838× fallback. Drivers `209`--`214`; compact results and memos `post-v1-{archive-root-gate-a,jshape-routing-gate-b}`; exact commands in `PROVENANCE.md` §8.27. |

Archive-format v1.3-v1.5 command transcripts are retained in `logs/v13-pipeline.sh`,
`logs/v14-pipeline.sh`, and `logs/v15-pipeline.sh`. The tagged Gravlax repository is the only
authoritative implementation source.

## Post-v1 work

Development after the frozen tag is tracked on the `post-v1-week2` branch. The current
specifications and results include:

- `docs/formal-conditional-sufficiency.md` — exactness vocabulary, proof scope, and
  necessary-field counterexamples;
- `docs/post-v1-week1-status.md` — local review commits, claim audit, and validation record;
- `docs/post-v1-red-team-review.md` — severity-ranked independent review, experiment slate,
  “do not pursue” list, and decision-complete eight-week roadmap;
- `experiments/specs/fair-baselines.md` — decision-complete matched storage/runtime experiment;
- `experiments/manifests/post-v1-gates.yaml` — machine-readable acceptance gates.
- `scripts/110_make_replay_minimal_bam.sh` and `scripts/111_d0_fair_baseline_pilot.sh` — the
  accepted D0 ingest-equivalent transform and its fail-closed reproduction driver;
- `results/post-v1-d0-fair-baseline-pilot.json` and
  `memos/post-v1-d0-fair-baseline-pilot.md` — exact artifact identities, equivalence checks,
  preliminary timings, and interpretation.
- `scripts/112_d0_starsolo_counts_only.sh`, `scripts/113_d0_warm_counts_comparison.sh`, and
  `results/post-v1-d0-counts-only-pilot.json` — the no-BAM STARsolo arms and the five-run D0
  correction to the speed headline.
- `experiments/specs/streaming-replay.md`, `results/post-v1-streaming-replay-{d0,d1}.json`, and
  `memos/post-v1-streaming-replay.md` — the exact bounded-memory replay gate and verdict.
- `experiments/specs/mouse-5p-generalization.md`, `experiments/manifests/mouse-5p.yaml`, and
  `results/post-v1-mouse-5p-pilot.json` — the public mouse 5′ v2 Gene pass, marginal velocity
  boundary, matched-CRAM comparator, and complete rerun protocol.
- `experiments/specs/downstream-fidelity.md`, `results/post-v1-downstream-fidelity.json`, and
  `results/post-v1-downstream-cluster-diagnostics.json` — the preregistered downstream gate,
  including its negative independent-partition result, stable marker/velocity consequences, and
  separately labeled post-hoc failure localization.
- `experiments/specs/candidate-normalized-convex-em.md` and
  `results/post-v1-candidate-normalized-convex-em-d0.json` — the D0 development grid that repaired
  whole-transcriptome pseudo-count dilution but found only a small real-over-shuffled group gain.
- `experiments/specs/dirichlet-proxy-screen.md` and
  `results/post-v1-dirichlet-proxy-screen-d0.json` — the locked posterior-mean adaptive-pooling
  screen, including matched fixed-convex controls, evidence-depth strata, and its marginal STOP
  decision before a full hierarchical Dirichlet implementation.
- `experiments/specs/depth-gated-hybrid-em.md` and
  `results/post-v1-depth-hybrid-em-d0.json` — the one-parameter monotone transition between fixed
  convex and proxy predictions, its successful depth repair, and the single predeclared
  real/shuffled interaction criterion that retained a marginal verdict.
- `experiments/specs/depth-hybrid-em-confirmation.md` and
  `results/post-v1-em-confirmation-d4/` — the remotely locked untouched D4 factorial confirmation,
  replay-only candidate-excluded group selection, raw seed-level metrics, and the sole Brier
  guardrail failure that closes the branch before D3 or a full hierarchical Dirichlet model.
- `experiments/manifests/cell-group-junction-set-query.yaml`, scripts `175`--`176`, and
  `results/post-v1-cell-group-junction-set-query.json` — the strict shared cell/group scope,
  conservative class-level junction-set semantics, D0/D1 exactness checks, preserved unscoped
  interfaces, and the passing union-query runtime/memory gate.
- `experiments/manifests/federated-biological-screen.yaml`, scripts `179`--`181`, and
  `results/post-v1-federated-biological-*` — the exploratory four-PBMC and six-shard
  whole-genome event screens, exact STAR-SJ reducer audits, literature-positive controls, ranked
  compact tables, and the prospectively locked untouched FNBP1 normal-brain confirmation.
- `results/post-v1-fnbp1-normal-brain-confirmation.json` and scripts `182`--`184` — the passing
  untouched test: 74 informative UMIs, 97.3% archive inclusion, identical targeted
  annotation-free STAR archive totals, and a 96.8% STAR-SJ inclusion proxy. The recorded limits
  distinguish prospective biological replication from a whole-transcriptome mapping test.
- `experiments/manifests/fnbp1-broad-cell-localization.json`, scripts `186`--`189`, and
  `results/post-v1-fnbp1-broad-cell-localization.{json,tsv}` — frozen marker-rule localization of
  the fixed event in D5. Rates use all 10,261 H5 nuclei; the two exclusion molecules fall in the
  four-molecule microglia/immune stratum, a descriptive lead rather than a differential claim.
- `results/post-v1-federated-query-workflow-optimization.json`, scripts `180`--`181`, and the
  associated memo — exact row-denominator pushdown, packed group reduction, bounded chromosome
  concurrency, and one-event-at-a-time summary parsing, with byte-identical locked candidate TSVs.
- `scripts/205_validate_federated_sparse_output.py`,
  `results/post-v1-federated-sparse-output-validation.json`, and the associated memo — an
  end-to-end reconstruction of the matched six-archive cohort JSON from the versioned sparse
  tables and explicit dimensions. All 543 events and 3,258 event–sample observations match; the
  final `annotation_genes_json` bundle is 76,970 bytes, 45.27× smaller than the 3,484,594-byte
  dense JSON, at effectively neutral single-run wall time and RSS.
- `scripts/206_benchmark_collection_index.sh`, `scripts/207_validate_collection_index.py`,
  `results/post-v1-collection-index-v2.json`, and the associated memo — an optional exact sidecar
  over eight 1,326,036,691-byte SEZ archives. The deterministic fresh root is 8,931,229 bytes
  (0.674%); the immutable four-plus-four chain is 9,469,771 bytes (0.714%). Forty ordinary
  per-archive comparisons and five fresh-root/chain comparisons are exact. Sparse and absent
  routes prune six and eight archives, and a support bound prunes all eight without opening one.
- `experiments/manifests/fnbp1-multidonor-sez-validation.json`, its protocol specification,
  scripts `190`--`193`, and `results/post-v1-fnbp1-multidonor-sez-*` — a remotely registered,
  donor-replicated fixed-event gate for eight adult human SEZ donors. Event-blind feasibility and
  checksum-verified acquisition precede all event queries. The result is an underpowered STOP:
  one donor reaches the ten-molecule primary minimum and none reaches the microglial minimum.
  Descriptively, seven donors contribute 39 inclusion and zero exclusion molecules; the independent
  STAR-SJ reduction agrees exactly per donor. See `memos/fnbp1-multidonor-sez-validation.md`.
- `experiments/manifests/molecular-splice-graph-v1.json`, script `194`, and
  `results/post-v1-molecular-splice-graph-v1.json` — the passing exact per-archive graph primitive:
  strand-aware junction edges, class-deduplicated molecular path fragments, strict scopes, hard
  overflow, and explicit lower-bound/non-transcript semantics. The real donor-G FNBP1 check finds
  three two-edge patterns supported by six UMIs in a 4.8-ms median query; this is a capability
  observation, not a population contrast. See `memos/molecular-splice-graph-v1.md`.
- `experiments/manifests/replicate-aware-splice-graph-v1.json`, its protocol specification, the
  strict eight-donor design, script `195`, and
  `results/post-v1-replicate-aware-splice-graph-v1.json` — the passing population-layer
  engineering gate. Exact zero-explicit sample × path/edge matrices feed a biological-sample
  beta-binomial contrast; duplicate archives, low depth, malformed designs, reference mismatch,
  and overflow fail closed. The real FNBP1 panel remains counts-only and preserves the registered
  donor-level STOP. See `memos/replicate-aware-splice-graph-v1.md`.
- `experiments/manifests/sez-polyasite-mixture-posthoc-v1.json`, its method-audit specification,
  scripts `197`, `199`--`203`, and `results/post-v1-sez-polyasite-mixture*` — a protocol-aware
  replacement for invalid cleavage-coordinate interpretation of fragmented 10x 3′ boundaries.
  PolyASite identities, genomic internal-priming rejection, leave-one-donor-out empirical kernels,
  paired inference, 12/24/48-bp controls, a shuffled-label null, and simulation all pass. The
  optimized primary run is byte-identical and takes 22.02 s, 5.18× faster than eight legacy scans.
  The exploratory NTRK2 result shifts from 20.0% to 78.1% MANE/full-length-like terminal-group
  usage across all eight donor pairs; the preregistered global-lengthening threshold remains failed.
  See `memos/sez-polyasite-mixture.md`.
- `experiments/specs/adaptive-evidence-frontier.md`,
  `results/post-v1-adaptive-evidence-frontier.json`, and
  `memos/adaptive-evidence-frontier.md` — a rejected annotation-independent geometry pilot. The
  all-chain frontier improves D0 Gene movement by 0.02833 percentage points and mean Velocity
  error by 0.41326 points, but costs 19.272% more bytes and 1.252× replay time. The sole
  pre-declared junction-only rescue clears the resource limits (+2.810%, 1.025×) but loses the
  accuracy gain, so the frozen verdict is **STOP**. No implementation was merged into Gravlax and
  the default two-extreme archive bytes remain unchanged.
- `scripts/204_summarize_bounded_polyasite.py`,
  `results/post-v1-polyasite-bounded-memory.json`, and
  `memos/post-v1-polyasite-bounded-memory.md` — the exact eight-donor PolyASite analysis with a
  chromosome/chunk/sample working-set frontier. The three scientific tables are byte-identical;
  peak RSS falls from 3,660,468 to 2,262,980 KiB (38.2%) while wall time remains inside the frozen
  30-s gate at 28.21 s. The bound is measured and structural, but not asymptotically constant in
  sequencing depth because one chromosome's selected molecule classes remain resident.
- `experiments/specs/2026-09-01-sparse-terminal-tail*.md`,
  `results/post-v1-sparse-terminal-tail.json`, and
  `memos/sparse-terminal-tail-evidence.md` — a sequence-derived optional-evidence pilot whose
  retained sidecars store no nucleotide sequence or quality, with three intentionally separate
  outcomes. The eight-donor general signal passes: 53,346 external-site
  witnesses, 10.299× external/control odds, +0.127980 AUC over clip length, internal-priming
  depletion, and 0.126073% sidecar/archive bytes, with every donor clearing the frozen gates. The
  NTRK2 check fails with only three witnesses in one donor. The proposed one-u16-per-molecule
  archive layout also hard-stops because current molecule reduction can discard the positive
  read's endpoint and can map multiple anchors to one `MolRec`; the exact redesign is a sparse
  per-selected-ordinal event list carrying signed anchor deltas and one u16 per event. No pilot
  code was merged into Gravlax and default archive bytes remain unchanged.
- `experiments/specs/2026-09-01-archive-root-and-jshape-routing.md` and
  `experiments/manifests/archive-root-jshape-routing.yaml` — two ordered systems gates. Gate A
  passes: its exact encoded-section commitment adds 355,840 bytes (0.026835%) across the eight SEZ
  archives, reduces collection-build source I/O from 1,336,022,077 to 10,341,234 bytes (129.19×),
  and preserves exact Gene/Velocity replay and all fixed query answers. The compact result,
  engineering interpretation, and 59-row excluded-artifact inventory are
  `results/post-v1-archive-root-gate-a.json`, `memos/post-v1-archive-root-gate-a.md`, and
  `manifests/archive-root-gate-a-artifacts.tsv`. After the Gate-A result commit, script `211`
  materialized the prospectively selected 96-junction Gate-B panel: 24 rows in each of the one,
  two-to-three, four-to-seven, and eight-archive presence strata, with SHA-256
  `55980706cd6fbf0c07a1e652b330e5aca6946b1365010c8ef4e8cb6d6c4f1ca1`. The
  `experiments/designs/archive-root-jshape-routing-panel.*` trio and
  `experiments/specs/2026-09-01-jshape-routing-serialization-addendum.md` are frozen together. No
  Gate-B route code or measurement existed at this freeze boundary. After the route interface was
  finalized, scripts `212`--`214` locked the same-binary fallback/route campaign, strict reducer,
  named Rust predicate evidence, and independent malformed-collection fixtures. Exact commands
  and schemas are in `experiments/specs/2026-09-02-jshape-routing-gate-b-runbook.md`. The sole
  confirmatory run, r3, passes: the optional route adds 30,354,712 bytes (2.289131% of source
  bytes), returns exact root/chain answers for the five controls and all 96 prospective panel
  junctions, reduces aggregate panel wall time to 0.386838× fallback, and remains inside the
  construction limit at 1.900709× wall and 607,132 KiB. r1 was an empty preflight and r2 a
  diagnostic performance failure used before optimization; neither supplies a reported estimate.
  The compact result, interpretation, and exact excluded inventory are
  `results/post-v1-jshape-routing-gate-b.json`, `memos/post-v1-jshape-routing-gate-b.md`, and
  `manifests/jshape-routing-gate-b-artifacts.tsv`.

## Preprint follow-up experiments (2026-09-11)

Five experiments were added for the bioRxiv preprint after the post-v1 work above. All are
summarized in `memos/paper-followup-20260911.md`, and their compact records (summaries,
`/usr/bin/time -v` output, ingest reports, run logs) are in `results/paper-followup-20260911/`.
The large run trees (STAR alignments, archives, per-chromosome discovery JSON) are excluded under
the large-artifact policy.

| Result family | Primary drivers | Where it appears |
|---|---|---|
| D2′ GeneFull timing (STARsolo GeneFull-only versus native GeneFull replay, protocol of script `120`) | `scripts/215_d2p_genefull_timing.sh` | Supplementary runtime protocol |
| D1 `--geometry-fidelity` archive: size, Gene replay, velocity components | `scripts/216_d1_geometry_fidelity.sh` | Supplementary Note S1 (read reduction) |
| GENCODE v32/v48/v49 transcript and junction stability, observed support on release-absent junctions | `scripts/217_annotation_junction_stability.py` | Supplementary annotation-stability note |
| D0 junction-seeded ingest alignment (v32/v49 seed, one- and two-pass) | `scripts/218_junction_seeded_ingest.sh`, `scripts/219_junction_seeded_v49_onepass_arm.sh` | Supplementary junction-seeding table |
| Cohort-wide cassette discovery over the eight-donor SEZ collection versus a naive per-archive discovery-and-merge baseline; rerun with release v0.2.3 | `scripts/220_collection_discovery.sh`, `scripts/221_collection_discovery_rerun.sh`, `scripts/222_naive_discovery_merge.py`, `scripts/223_collection_discovery_v023.sh`, `scripts/224_naive_discovery_v023.sh` | Section 3.4 and Supplementary Note S4 |

## Large-artifact policy

No large artifact is tracked by Git. `manifests/publication-artifacts.tsv` records the current
publication archives and CRAM baselines with SHA-256 checksums. Regenerable BAMs and public FASTQ
bundles are recorded by dataset identity and byte size in `manifests/excluded-artifacts.tsv`.
`manifests/input-digests.sha256` contains the pre-existing checksums of the D0 reads, whitelist,
genome, and compressed GENCODE inputs. The rejected adaptive-frontier run tree is excluded; its
logical paths, sizes, and SHA-256 identities are retained in
`manifests/adaptive-evidence-frontier-artifacts.tsv`. The final excluded sparse-output and
collection-index run trees, source archive panel, sidecars, and executable are identified in
`manifests/excluded-artifacts.tsv` and their compact result records.

The final archive-root Gate-A binary, sealed archives, collection sidecars, replay matrices, and
confirmatory r3 run tree are pinned in `manifests/archive-root-gate-a-artifacts.tsv`. The incomplete
r1 and diagnostic-only r2 run trees are recorded in `manifests/excluded-artifacts.tsv`; neither is
used as confirmatory evidence.

The Gate-B binary, frozen panel, rooted input archives, confirmatory r3 run tree, and routed and
fallback collections are pinned in `manifests/jshape-routing-gate-b-artifacts.tsv`. Gate-B r1
(preflight only), r2 (diagnostic only), and r3 (the sole confirmatory PASS) are separately recorded
in `manifests/excluded-artifacts.tsv`.

The terminal-tail BAMs, archives, experimental sidecars, audits, final pilot executable, reference,
catalogue, and group maps are likewise excluded and pinned in
`manifests/sparse-terminal-tail-artifacts.tsv`.
