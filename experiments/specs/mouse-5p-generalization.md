# Mouse 5′ generalization gate

Status: accepted Gene pilot on 2026-08-31; velocity is marginal, and GEM-X v3 replication awaits
public raw evidence. The compact record is `results/post-v1-mouse-5p-pilot.json`.

## Outcome

- The same M39 replay from `.aie`, the annotation-free ingest BAM, and the exported
  post-correction molecule BAM is byte-identical. Archive- and BAM-sourced velocity artifacts are
  also byte-identical.
- Against fresh annotation-aware STARsolo on 1,129 current-M39 called cells, Gene union-L1/2 is
  0.751% of oracle UMI mass (PASS); replay retains 5,513,815 versus 5,522,226 oracle UMIs, and only
  2/1,129 cells change total count by more than 1%.
- The deliberately incomplete one-pass ingest moved 13.892% of UMI mass. Restoring the established
  annotation-free two-pass protocol reduced this to 0.751%; retain the failed pilot because it
  demonstrates that junction discovery, not archive serialization, is the binding evidence step.
- The wrong `Forward` orientation retains only 16.68% of the correct reverse-strand filtered-cell
  UMI mass, validating chemistry orientation as necessary replay configuration.
- Velocity is PASS for unspliced (0.939% relative L1) but MARGINAL for spliced (2.907%) and
  ambiguous (4.624%). This supports mouse 5′ Gene replay, not an exact or fully passing fresh
  velocity claim.
- The 62.40 MB archive is 14.82 bits/read pair. A function-matched post-correction CRAM is
  131.25 MB (2.103× larger), and five exact streamed-CRAM replays take a median 11.89 s versus
  1.79 s native replay (6.64× slower).

## Hypothesis

The existing evidence abstraction should support a public mouse 10x 5′ experiment without a
format change: genomic placement strand is already archived, 10 bp UMIs fit the existing packed
representation, and Gene/Velocyto replay needs only an explicit STARsolo strand relationship.
After selecting `Reverse`, archive-sourced and direct-BAM Gravlax replay must remain byte-identical.
Against a fresh annotation-aware STARsolo run on the same GRCm39/GENCODE M39 reference, Gene UMI
mass movement on oracle-called cells should be at most 1%. Velocity retains its predeclared
per-component relative-L1 threshold of 1.5%.

## Dataset choice

The resource-efficient primary gate is the official 10x 1k C57BL/6 mouse splenocyte 5′ v2
dataset (33,677,140 read pairs; one lane). It tests two previously unvalidated axes together:
mouse gene structure and the v2 10-base UMI. The producer reports 1,119 cells, but used
mm10-2020-A; that matrix is an external sanity check, not a numerical oracle for the current-M39
experiment.

The preferred replication is the official 5k Mouse PBMC GEM-X 5′ v3 dataset produced by Cell
Ranger 9.0.0. Its filtered matrix and molecule-info file are public, but a public raw-read or
alignment URL was not resolvable on 2026-08-30. Do not substitute its processed matrix for raw
evidence. Replicate when the producer exposes FASTQ/BAM access or supplies an accession.

## Design and controls

1. Build one annotation-free STAR 2.7.11b GRCm39 base index (no GTF at index time).
2. Run the fresh oracle with M39 supplied at mapping time, `--soloStrand Reverse`, the documented
   16+10 barcode geometry, CellRanger4 adapter clipping, and Gene/GeneFull/SJ/Velocyto features.
3. Run annotation-free two-pass ingest from the same reads/base index with no GTF and only `SJ`
   enabled, retaining secondary placements and raw barcode/UMI tags. Two-pass discovery is part of
   the established Gravlax ingest protocol, not a post-hoc mouse-specific adjustment.
4. Build `.aie` from the ingest BAM with Gravlax commit recorded in the manifest.
5. Replay M39 with `--solo-strand reverse` from both `.aie` and the ingest BAM. Require
   byte-identical matrix, feature, and barcode artifacts.
6. Compare archive replay with fresh STARsolo only on the oracle's filtered cells, matching genes
   by unversioned Ensembl identifier. Report union-gene L1/2 as moved UMI mass and decompose net
   gain from reassignment.
7. Replay Velocyto with the same strand choice and apply the frozen VELO-1 component thresholds.
8. Run the default-forward replay as a negative control. It should lose most Gene mass, showing
   that chemistry orientation is necessary analysis configuration rather than an innocuous label.

All timed commands use 24 threads on CPUs 0–23. STAR coordinate sorting requires a child-only
soft file-descriptor limit of 65,535 on this host; the default 1,024 limit failed before FASTQ
processing and that failed log is retained. Record wall, CPU, RSS, file-system I/O, exact command
lines, input/output hashes, archive/BAM bytes, archive bits per read pair, cell count, and matrix
mass.

## Decision gates and interpretation

- **STOP:** archive and direct-BAM Gravlax artifacts differ. This is an implementation/source-path
  defect, independent of comparison with STARsolo.
- **PASS Gene:** filtered-cell union moved mass ≤1%; **MARGINAL:** ≤5%; **STOP:** >10%.
- **PASS velocity:** each component relative L1 ≤1.5% and total mass delta ≤5%; retain the frozen
  marginal/stop interpretation from VELO-1.
- A failure caused by reverse-strand assignment is a tool-capability defect; a residual fresh-STAR
  difference after exact fixed-alignment replay is alignment/reduction deviation and must be
  described that way.
- Passing this pilot supports “mouse 5′ v2” rather than “all organisms and chemistries.” The 5k
  GEM-X v3 replication is required before making a broad 5′ claim.

## Paper value

A pass removes the manuscript's largest generalization caveat with a compact public experiment
and demonstrates that the store retained the molecule evidence required by a genuinely different
library orientation. The wrong-strand negative control also sharpens the conditional-sufficiency
story: archived geometry is sufficient only together with the chemistry's assignment semantics.
