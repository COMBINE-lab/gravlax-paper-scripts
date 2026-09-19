# Gravlax paper scripts

This repository archives the experiment drivers, compact result records, and provenance behind
the Gravlax manuscript, *Gravlax: an annotation-independent molecular evidence archive for
single-cell RNA-seq*. It accompanies the software at
[COMBINE-lab/gravlax](https://github.com/COMBINE-lab/gravlax); the results in the paper were
produced with release [v0.2.3](https://github.com/COMBINE-lab/gravlax/releases/tag/v0.2.3)
(the follow-up experiments in `scripts/215`--`224` were first run with v0.2.2 and the cohort
discovery was rerun with v0.2.3; the two produce identical output).

The repository deliberately does **not** contain the Gravlax implementation, the manuscript
source, raw reads, alignments, reference indexes, or molecular archives. It contains what is
needed to see exactly how each reported number was obtained and to rerun the analysis given the
public inputs listed in `PROVENANCE.md`.

## Layout

| Path | Contents |
|---|---|
| `scripts/` | The experiment drivers, exactly as executed. Numbered in the order they were written. |
| `results/` | Compact machine-readable outputs (JSON, TSV, small summaries) that the reported numbers are read from. |
| `logs/` | Command transcripts and `/usr/bin/time -v` records for the timing and storage measurements. |
| `PROVENANCE.md` | The experiment catalogue: datasets, references, how the STARsolo oracle was run, and one entry per experiment with its inputs, command lines, and results. |
| `experiments/` | Pre-registered specifications (`specs/`), machine-readable acceptance criteria (`manifests/`), and frozen query panels (`designs/`) written before the corresponding measurements. |
| `docs/decisions/` | Thresholds fixed before each measurement. |
| `docs/`, `memos/` | Chronological decision log, method notes, and negative results. |
| `manifests/` | SHA-256 identities of inputs and of the large artifacts that are excluded from Git. |
| `environment/setup.sh` | The environment setup used on the analysis machine, kept as a record rather than as an installer. |
| `tests/` | Small fixture tests for the fail-closed result reducers. |

## Running the scripts

Every driver locates the analysis project through one environment variable,
`GRAVLAX_PROJECT_ROOT`, the directory holding the `data/`, `annotations/`, `ref/`, and `runs/`
trees described in `PROVENANCE.md`. Shell drivers fail immediately if it is unset. Two optional
variables select tools: `GRAVLAX_PYTHON` (the interpreter for the analysis environment, default
`python3`) and `GRAVLAX_ZSTD` (the `zstd` binary used by script `211`, default `zstd` on `PATH`).

```bash
export GRAVLAX_PROJECT_ROOT=/path/to/analysis-project
export GRAVLAX_PYTHON=/path/to/env/bin/python   # optional
```

Records under `results/`, `logs/`, `PROVENANCE.md`, and `memos/` write the same root as the
literal `${GRAVLAX_PROJECT_ROOT}` wherever a transcript or artifact identity contains an
absolute path; `$P` in `PROVENANCE.md` means the same thing. A few records name locations
outside the project root in the same style: `${GRAVLAX_SRC}` (the Gravlax source clone),
`${GRAVLAX_PAPER_SRC}` (the manuscript clone), `${GRAVLAX_TOOLS}` (the Node toolchain for the
docs site), `${CONDA_ROOT}`, `${MALVA_ROOT}` (the Malva 1.0.0 checkout used as a storage
baseline), `${STAR_SOURCE}` (the STAR source tree consulted for counting-rule fidelity), and
`${WHITELIST_SOURCE}` (the directory holding the 10x `3M-february-2018.txt` whitelist). The
analysis machine appears as `analysis-host`; its hardware is described in `PROVENANCE.md`.

The minimal end-to-end path on the PBMC 1k dataset ($D_0$) is: build the STAR indexes
(`05`), run the STARsolo reference (`10`), run the annotation-free ingest alignment (`40`), build
the archive and replay it (commands in `PROVENANCE.md` §10), then compare against the reference
with the fidelity scripts below.

## Where each part of the paper comes from

Section numbers refer to the manuscript; "Note" refers to the supplementary notes.

**Building references, oracles, and archives**

- `05_build_anno_indices.sh`, `10_oracle_starsolo.sh`: one STAR index per GENCODE release and
  the STARsolo runs that every replay is compared against.
- `40_ingest_align.sh`, `60_junction_catalogue.sh`: the annotation-free, barcode-aware ingest
  alignment and the cross-sample junction catalogue.
- `90_d1_reverify.sh`, `96_d2_campaign.sh`, `97_d2p_campaign.sh`, `123`--`125`,
  `197_build_sez_full_archives.sh`: the same pipeline on the PBMC 5k, glioblastoma, brain
  nuclei, mouse 5′, and eight-donor SEZ datasets (Note S8 lists every dataset).

**Replay fidelity and what explains the residual (Sections 3.1, 3.2; Notes S1, S4)**

- `30_gateb_compare.py`, `31_delta_capture.py`: how much a count matrix changes between
  GENCODE releases and whether replay captures that change.
- `91_velo_gate.py`, `94_velo2_strata.py`: velocity replay per component and the attribution
  of the residual to the two-representative read reduction.
- `216_d1_geometry_fidelity.sh`: the `--geometry-fidelity` archive that stores each accepted
  read geometry once, used to measure how much of the deviation the reduction itself causes.
- `217_annotation_junction_stability.py`, `218`--`219`: stability of transcript and junction
  sets across GENCODE releases, and the junction-seeded ingest alignment experiment.
- `126`--`132`: downstream stability of clusters, markers, and velocity under replay,
  with the clustering settings selected on the reference matrices alone.

**Storage and runtime (Section 3.3; Note S4)**

- `50_gatee_baselines.sh`, `110`--`120`: capability-matched storage baselines (BAM, CRAM,
  post-correction molecule BAM/CRAM) and the randomized, cache-controlled timing comparisons
  between fresh STARsolo and native replay.
- `215_d2p_genefull_timing.sh`: the GeneFull timing comparison on brain nuclei.
- `121`--`122`: bounded-memory streaming replay against eager replay.
- `142`--`143`, `171`--`176`: hot-path, compiled-annotation, batched-query, and cell-group
  query cost measurements.

**Federation, junction-by-shape indexing, and cohort discovery (Section 3.4; Notes S2, S3, S4)**

- `103_atlas_build.sh`, `206`--`207`: building and benchmarking the multi-archive collection
  index.
- `209`--`214`, `archive_root_*.py`, `jshape_routing_gate_common.py`: content-addressed
  archive roots and junction-by-shape routing, with their exactness checks and adversarial
  fixtures.
- `177`--`181`, `205`: the event engine, genome-wide event federation across archives, and
  the sparse cohort output.
- `220`--`224`: cohort-wide cassette discovery over the SEZ collection and the naive
  per-archive discovery-and-merge baseline it is compared with.

**Biological analyses (Section 3.5; Notes S5, S6)**

- `179`--`185`: the splice-event screen between T cells and monocytes and the FYB1 control.
- `182`--`193`, `194`--`195`: the FNBP1 exon-inclusion result, its confirmation in an
  independent normal-brain dataset, localization to broad cell classes, the eight-donor SEZ
  validation, and the molecular splice graphs.
- `70`--`82`, `92_discr_gate.py`, `98_prospective.py`, `136`--`141`: withheld-annotation
  panels, recovery of genes added in a later annotation, and discovery recall.
- `102_pool_extension.py`: reference optimization as replay (Pool et al.).
- `80`--`81`, `93_apa_diff.py`, `100_apa_replicate.py`, `105_ptprc.py`, `196`--`204`,
  `208_aggregate_sparse_terminal_tail.py`: terminal fragment-boundary redistribution, the
  PolyASite mixture, the NTRK2 terminal groups, and non-templated poly(A)-tail evidence.

**Inference on multi-gene evidence (Section 3.6; Note S7)**

- `95_em0_gate.py`, `99_em_af_concordance.py`, `104_em_de.py`: pooled EM against STARsolo's
  multimapper EM and against alevin-fry, and its effect on differential expression.
- `133`--`135`, `144`--`170`: the packed cell-sharded EM, masked recovery and leakage controls,
  candidate-normalized partial pooling, depth-dependent combination, and parameter selection on
  a development dataset with confirmation on untouched ones.

Every script's header comment states what it measures; `PROVENANCE.md` §8 gives the command
lines, inputs, and results for each experiment in the same order.

## Large artifacts

No large artifact is tracked by Git. `manifests/publication-artifacts.tsv` records the
publication archives and CRAM baselines with SHA-256 checksums; `manifests/excluded-artifacts.tsv`
records regenerable BAMs, public FASTQ bundles, and run trees by identity and byte size; and
`manifests/input-digests.sha256` holds the checksums of the reads, whitelist, genome, and
compressed GENCODE inputs.

## License

See `LICENSE`.
