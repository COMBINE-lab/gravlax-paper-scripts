# Follow-up experiments for the bioRxiv preprint, 2026-09-11

Both run with the released 0.2.2 binary on analysis-host; records under
`$P/runs/paper-followup-20260911/` (`$P` = the analysis project root).
Compact records (summaries, `/usr/bin/time -v` outputs, ingest reports, run logs) are archived in
`results/paper-followup-20260911/`; the large per-run trees (STAR alignments, archives, discovery
JSON) are excluded under the large-artifact policy. The 0.2.2 and 0.2.3 binaries are the released
Gravlax tags v0.2.2 and v0.2.3.

## D2' GeneFull timing (`scripts/215_d2p_genefull_timing.sh`)
Protocol as gravlax-paper-scripts/120: cold arm per method, then five warm
blocks, alternating, STAR first in odd blocks, 24 threads on CPUs 0-23.
STARsolo: GeneFull only, --outSAMtype None. Replay: `replay-rows --gene-full`
with the v49 GTF. All 12 runs exact against their references.

| arm | warm median wall | warm peak RSS |
|---|---:|---:|
| STARsolo GeneFull | 361.46 s | 32.55 GiB |
| replay GeneFull | 9.27 s | 5.79 GiB |
| ratio | 38.99x | 5.6x |

## D1 geometry fidelity (`scripts/216_d1_geometry_fidelity.sh`)
`ingest-archive --geometry-fidelity --zstd-level 19` on the D1 ingest BAM:
7:26 wall, 40.2 GB peak RSS. Sizes: compact 544,022,039 B (11.34 bits/read),
fidelity 692,731,515 B (14.44 bits/read).

| metric (5,038 called cells, v49) | compact | fidelity |
|---|---:|---:|
| Gene moved UMI mass | 0.307% | 0.253% |
| velocity spliced rel L1 | 2.672% | 2.354% |
| velocity unspliced rel L1 | 0.859% | 0.852% |
| velocity ambiguous rel L1 | 6.132% | 4.669% |

Compact replays under 0.2.2 reproduce runs/replay/d1-aie and d1-velocity
byte for byte. Under fidelity every entry is in the complete stratum.

## Annotation stability (`scripts/217_annotation_junction_stability.py`)
Summary in `$P/runs/paper-followup-20260911/annotation-stability/summary.md`.
v32->v49: transcripts by id Jaccard 0.44, by intron chain 0.42, junctions 0.61;
v32 junctions retained in v49: 379,610/382,880 (99.1%). Observed junction reads
(unique, >=3 reads) on v49-not-v32 junctions: D0 0.24%, D1 0.32%, D2' 0.73%.

## Junction-seeded ingest, D0 (`scripts/218_junction_seeded_ingest.sh`, `219_junction_seeded_v49_onepass_arm.sh`)
| arm | Gene v49 | Gene v32 | spliced | unspliced | ambiguous | bytes |
|---|---:|---:|---:|---:|---:|---:|
| unseeded two-pass | 0.244% | 0.232% | 1.808% | 0.520% | 3.059% | 111,087,381 |
| v32 two-pass | 0.228% | 0.214% | 1.750% | 0.509% | 2.989% | 112,286,136 |
| v49 two-pass | 0.213% | 0.223% | 1.576% | 0.502% | 2.652% | 112,452,111 |
| v32 one-pass | 0.160% | 0.140% | 1.319% | 0.042% | 2.045% | 110,403,101 |
| v49 one-pass | 0.143% | 0.154% | 1.088% | 0.032% | 1.506% | 111,457,143 |
Seeds inserted with --sjdbFileChrStartEnd --sjdbOverhang 90 on ref/star-base; aie 0.2.2
ingest with --junction-discovery/--junction-catalogue/--alignment-annotation provenance.

## Cohort-wide event discovery, SEZ collection (`scripts/220_collection_discovery.sh`, `221_collection_discovery_rerun.sh`, `222_naive_discovery_merge.py`)
Outputs in `$P/runs/paper-followup-20260911/collection-discovery/`. Eight SEZ archives
(1,326,392,531 B), route-root collection (39,286,118 B) and coordinate-only fallback-root
(8,931,406 B); groups from post-v1-sez-transcript-end-feasibility donor-*.tsv (9 groups);
`collection find-events --kinds cassette --min-support 2 --min-samples 4 --min-donors 4
--min-umi-classes 8 --min-side-umi-classes 2 --min-group-umi-classes 2
--require-groups astro_nsc,mature_neuron`, 24 threads pinned, `/usr/bin/time -v`.

| arm | annotation filter | wall (time -v) | total_seconds | peak RSS | archive bytes read | candidate defs | candidate entities | retained |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| routed | novel vs v32 | 1:28.97 | 83.6 | 12.44 GB | 1,017,509,893 | 5,413,686 | 2,976,274 | 205,735 |
| coordinate-only | novel vs v32 | 1:33.36 | 88.0 | 12.30 GB | 1,017,509,893 | 5,413,686 | 2,976,274 | 205,735 |
| routed | none | 1:31.48 | 85.9 | 12.34 GB | 1,017,509,893 | 5,413,686 | 2,976,274 | 209,570 |

Pooled catalogue: 1,344,270 junctions; exact_match_attempts 324,154,860; annotation
comparisons 44,294,758; sidecar bytes read 9,605,128. Genome-wide discovery is a scan:
the junction-by-shape routes change neither bytes read nor time (they pay off for panels).
Novel-vs-v32 retained entities by primary gap class: strand 125,005; boundary 76,172;
missing_junction 3,442 (donors: 8: 2,408; 7: 583; 6: 272; 5: 122; 4: 57); overlap 1,116.
Distinct overlapping genes: 5,083 overall, 2,843 in the missing_junction class.
FNBP1 cassette `cassette:chr9:-:129908999-129915965,129919242-129923843,129908999-129923843`:
exact 44 UMI classes, 44 cells, 6 samples/6 donors, gap class missing_junction,
catalogue donor route upper bound 6.

Naive baseline (`query events --kind cassette` per archive per chromosome, 25 chromosomes,
same thresholds where they apply): 200 runs, 22.9 s total wall (A 1.6, B 1.5, C 1.0, D 5.5,
E 1.2, F 1.4, G 9.4, H 1.2 s); merge 8.7 s; union 552,167 events; 1,196 pass the cohort rule
(>=4 samples, >=8 classes, >=2 per side, both required groups with >=2). The FNBP1 cassette
appears only in D's and G's per-archive output, because per-archive candidate generation
needs all three components in that archive's own catalogue; the merge therefore sees two
donors and rejects it, whereas the collection counts exact support in six.

### Rerun with the released 0.2.3 binary (`scripts/223_collection_discovery_v023.sh`, `224_naive_discovery_v023.sh`)
Outputs in `$P/runs/paper-followup-20260911/collection-discovery-0.2.3/`; binary built from tag v0.2.3
(`$P/runs/gravlax-release-0.2.3-20260912/aie-0.2.3`). Output identical to 0.2.2 after stripping timing fields.
| arm | wall (time -v) | peak RSS | total_seconds | stage: discover / exact / annotation |
|---|---:|---:|---:|---|
| routed novel-v32 | 0:08.82 | 2.98 GB | 6.4 | 1.41 / 3.86 / 0.41 s |
| coordinate-only novel-v32 | 0:08.84 | 2.97 GB | 6.4 | 1.38 / 3.86 / 0.40 s |
| routed all | 0:08.97 | 2.96 GB | 6.5 | 1.41 / 3.95 / 0.41 s |
Naive baseline with 0.2.3: 200 runs 23.9 s (A 1.6, B 1.5, C 1.0, D 5.9, E 1.2, F 1.4, G 10.0, H 1.2), merge 8.9 s,
union 552,167, recurrent 1,196, FNBP1 absent (unchanged).
