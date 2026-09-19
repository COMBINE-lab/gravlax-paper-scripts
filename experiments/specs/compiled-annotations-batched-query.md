# Compiled annotations and batched query plans

**Locked 2026-08-31 before implementation or candidate measurement.** Existing GTF and
single-query behavior are compatibility references.

## Compiled annotation

`aie compile-annotation annotation.gtf --out annotation.aic` will serialize exactly the model
already constructed by the `anno` crate: first-occurrence gene/chromosome identifiers,
transcripts, merged half-open exons, and the augmented per-chromosome overlap arrays. The
artifact has a fixed magic, version, payload length, and BLAKE3 payload checksum. Every existing
GTF argument will detect and load either representation; the CLI option names remain unchanged
for compatibility. The compiled artifact is an optimization of annotation loading, not a new
assignment model, and the source GTF remains authoritative.

Two independent compilations must be byte-identical. GTF and AIC paths must produce identical D0
Gene replay files, discovery output, and junction-annotation flags. Tests must reject bad magic,
future version, truncation, invalid indices, checksum mismatch, and trailing bytes. On five warm
D0 annotation-marked junction queries, the AIC median wall time must be no more than 25% of the
GTF path, peak RSS at most 1 GiB, and the artifact at most 5% of the 3.1-GiB GTF.

## Batched query plan

The v1 plan is strict TSV with header `id<TAB>kind<TAB>locus`; `kind` is `region` or `junction`,
and loci retain current 0-based half-open syntax. The executor opens archive dictionaries and
junction metadata once, unions required chunk postings/ranges, decodes each required chunk once,
and evaluates every predicate assigned to that chunk. It returns one JSON object in plan order;
absent junctions are result rows rather than fatal errors. Region results preserve the existing
anchor semantics, and both kinds preserve per-cell, per-class UMI deduplication.

The fixed 32-query D0 plan is deliberately a high-locality gene-panel workload: 16 distinct
windows and 16 exact junctions in one chromosome neighborhood. Every aggregate and per-cell
answer must equal the corresponding independent command. Seven byte-identical batch repetitions
are required. Relative to separate processes, unique chunk decodes must fall by at least 75%,
median wall time to at most 35%, and peak RSS must remain at most 1 GiB. These thresholds apply
only to the registered locality profile; dispersed plans still avoid repeated archive open and
metadata work but may share few chunks.

Region tracks, APA, genome-wide APA testing, discovery, and junction enumeration are deferred:
their outputs require distinct merge policies and should not be silently approximated by the v1
planner.

