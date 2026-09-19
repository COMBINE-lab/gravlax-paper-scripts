# Content-addressed local-shape routing: Gate-B result

## Verdict

**PASS.** A derived local-shape route accelerates exact point-junction and junction-set queries
without changing their molecular answers. The confirmatory r3 campaign matches routed and exact
fallback execution for all 101 fixed and prospective scientific cases on both a fresh eight-archive
root and an immutable A--D base plus E--H extension. All 26 positive and malformed-input predicates
also behave as frozen.

The route is a candidate index, never a count store. Within each exact intron-span bucket it retains
sorted unique `(local_shape_id, donor_offset)` pairs; the acceptor offset is recovered by checked
addition of the span. A routed hit is verified against the source molecule's shape and both genomic
coordinates, and the source molecule chunks remain authoritative for UMI-class reduction.

## Representation and size

The eight source archives occupy 1,326,036,691 bytes. The fresh route-bearing collection is
39,286,118 bytes, of which 30,002,228 bytes are compressed route blocks. Relative to the
8,931,406-byte route-free collection, the route premium is 30,354,712 bytes, or **2.289131%** of
the source archives. The immutable chain is 39,823,041 bytes versus 9,470,591 bytes route-free, a
30,352,450-byte (**2.288960%**) premium. Complete root and chain collections are 2.962672% and
3.003163% of source bytes. Compressed routes are 0.885093 times the prospectively frozen
33,897,246 compressed `shapes` bytes. Every frozen size gate passes.

The fresh root, repeated build, and reversed-input build are byte-identical within each arm. Route
bindings retain the immutable layer-local archive ordinal, source archive root, compressed
`shapes` digest, source shape count, and ordered block descriptors. Full inspection independently
reconstructs all 1,335,446 exact-span rows and 8,092,537 route pairs from the eight source
dictionaries.

## Exact scientific behavior

The 96 panel junctions were selected and committed before route implementation. Together with the
five fixed controls they give 101 exact routed-versus-fallback cases across root and chain forms.
The fixed totals are unchanged:

| query | exact result |
|---|---:|
| dense FNBP1 junction | 182 UMIs, 180 cells |
| sparse neighboring junction | 4 UMIs, 4 cells |
| absent neighboring junction | 0 UMIs, 0 cells |
| include-dense/exclude-sparse junction set | 182 include-only, 4 exclude-only, 0 both |
| NTRK2 region control | 74,848 molecules, 74,352 UMIs, 20,019 sample-cells |

This is an engineering exactness result over already archived evidence. It adds no biological
replication, does not improve the evidence retained by an archive, and does not compare archived
placements with a fresh annotation-aware aligner.

## Construction and query resources

All performance comparisons use one pinned binary, 24 threads, and eight alternating warm paired
blocks on the same host.

| operation | fallback median | route median | route/fallback |
|---|---:|---:|---:|
| collection build wall | 0.705 s | 1.340 s | 1.900709 |
| collection build peak RSS | 478,126 KiB | 607,132 KiB | 1.269816 |
| dense junction wall | 0.130 s | 0.040 s | 0.307692 |
| sparse junction wall | 0.120 s | 0.035 s | 0.291667 |
| junction-set wall | 0.130 s | 0.040 s | 0.307692 |
| region-control wall | 0.130 s | 0.135 s | 1.038462 |

Dense, sparse, and junction-set source-byte ratios are 0.105145, 0.129262, and 0.105145; including
the route sidecar, their total-logical-byte ratios are 0.135978, 0.198297, and 0.142038. Routed
peak RSS is 76,262 KiB for dense, 39,392 KiB for sparse, and 77,010 KiB for the junction set,
respectively 0.189193, 0.220777, and 0.190820 times fallback. The route-neutral region control has
source ratio 1.0, total-logical-byte ratio 1.010634, wall ratio 1.038462, and RSS ratio 0.996925.

Across the prospective 96-query panel, source-byte ratios have median 0.119379 and nearest-rank
p95 0.332421; total-logical-byte ratios have median 0.174608 and p95 0.528042. The maximum panel
source ratio is 0.589488. Aggregate panel wall time is 0.386838 times fallback. Thus every frozen
construction, I/O, wall-time, and memory threshold passes. The region result is a non-regression
control, not evidence that local-shape routes accelerate region queries.

## Integrity and failure behavior

The final fixture hook covers 11 accepted and 15 rejected cases. Accepted cases exercise
single-block shapes, multiple introns, repeated spans and offsets, distinct shapes sharing a span,
both strands, both retained representatives, multimappers, cross-chunk UMI classes, absent
coordinates, and exact fallback when a route is absent. Rejected cases cover unknown codecs,
wrong source roots or `shapes` digests, out-of-range shapes, reversed or inconsistent offsets,
unsorted or duplicate pairs, wrong archive ordinals or span buckets, count and checked-addition
overflow, truncation, checksum mismatch, and trailing bytes. Full route reconstruction and the
ordinary source-coordinate predicate are required independently of container checksums.

## Locked provenance and run history

The sole confirmatory run is `runs/post-v1/jshape-routing-gate-b-r3`. It records clean Gravlax
commit `7d997c6066b5a3d1efde9c466d5dce66d9fb8743`, clean reproducibility commit
`0de74cd901c39ce50898f3a0a23af05ecd23ef4b`, and the 101,954,648-byte release binary with SHA-256
`356d7699c163d31802e89d3ee54d7b6fccd176701a56bb99dbfd7f747088061a`. Its complete artifact
manifest binds 8,223 files and 647,307,809 bytes with SHA-256
`70b48020f62a99437027332cd24ceb493c9bdd8c4debbbc8819b4693cefef52e`. The compact result is
`results/post-v1-jshape-routing-gate-b.json` (SHA-256
`e3c7430ebed30d67711a22e70e6bc05a6d22726218cc4da927f83b646d372848`); its 37-row excluded
artifact inventory is `manifests/jshape-routing-gate-b-artifacts.tsv` (SHA-256
`b4e8d11ccbd64442bc2b6f921195998b7db66eb0fec41d64c61c6ed5cd4ffcf4`). An independent second
reduction to temporary under-project outputs returned PASS and those temporary outputs were
removed.

Earlier runs are retained only as development provenance:

- r1 is a preflight-only empty run skeleton. It contains no files, protocol, measurement,
  completion record, or artifact manifest and supplies no estimate.
- r2 completed the workload with the then-current implementation and was used diagnostically.
  Its route construction wall ratio was 2.528169, route construction RSS was 871,494 KiB, and
  region wall ratio was 1.153846, failing the frozen 2.0, 786,432-KiB, and 1.05 limits. This
  motivated implementation optimization, so r2 is not a confirmatory arm. Its manifest binds
  8,223 files and 647,324,075 bytes with SHA-256
  `ce2076e52dacfd085d4fd639cdc83c3b2e877934ff61a475e48eeb81c69f42de`.

## Decision boundary

Gate B licenses the optional content-addressed, local-shape-routed collection for exact
point-junction and junction-set acceleration over rooted archives. It does not license replay or
region acceleration, a general cold-filesystem latency claim, a biological-validity claim, or use
of routes as authoritative molecular counts.
