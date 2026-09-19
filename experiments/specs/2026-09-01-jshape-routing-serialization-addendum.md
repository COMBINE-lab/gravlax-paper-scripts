# Frozen Gate-B local-shape serialization addendum

**Status: FROZEN AFTER GATE-A PASS AND BEFORE GATE-B IMPLEMENTATION OR MEASUREMENT.** This addendum
is frozen together with the 96-row prospective panel. It changes no acceptance threshold, stop
rule, input, control, or query.

## Freeze record

Gate A is recorded by reproducibility commit
`9cba7af125dfce2d0b0c0b14ea56e92e00215a5d`. The committed Gate-A result
`results/post-v1-archive-root-gate-a.json` has SHA-256
`3f75fb1f7d5cc8b726405017a320c3698b7384c5da1f833f6d1f0a29453429db` and status `PASS`.

The prospective Gate-B artifacts frozen with this addendum are:

- `experiments/designs/archive-root-jshape-routing-panel.tsv`: 96 rows, 11,578 bytes, SHA-256
  `55980706cd6fbf0c07a1e652b330e5aca6946b1365010c8ef4e8cb6d6c4f1ca1`;
- `experiments/designs/archive-root-jshape-routing-panel.sha256`: SHA-256
  `db723b3883fdd9c8a57a4f74090aa3ddf68a018c3d6feccd54c42309180092b1`; and
- `experiments/designs/archive-root-jshape-routing-panel.metadata.json`: SHA-256
  `c3a0cd6d7cd53dd5e734f425f3dbc5c106911a0bda3ab4627dcee1dad2031d37`.

The panel contains exactly 24 rows in each frozen presence stratum. The metadata records 1,139,256
eligible one-archive junctions, 114,543 eligible two-to-three-archive junctions, 68,436 eligible
four-to-seven-archive junctions, and 22,035 eligible eight-archive junctions. No Gate-B route code
or result existed when these artifacts were materialized.

## Compact exact-span bucket

The original prospective specification writes a route tuple as
`(local_shape_id, donor_offset, acceptor_offset)` inside an exact intron-span bucket. The compact
serialization fixed by this addendum stores only the sorted unique pair

```text
(local_shape_id, donor_offset)
```

and reconstructs

```text
acceptor_offset = donor_offset + exact_bucket_span
```

with checked unsigned addition. Construction rejects an overflow, a shape id outside the bound
source dictionary, an offset that does not identify the stated intron in that shape, duplicate or
unsorted pairs, a pair in the wrong span bucket, and any mismatch found by full reconstruction
from the bound source `shapes` section. The on-disk payload of each exact-span bucket is therefore
the sorted unique sequence of `(local_shape_id, donor_offset)` pairs; `acceptor_offset` is never
stored. Query execution still verifies both genomic coordinates against the source representative;
the route never supplies a count.

## Equivalence to the frozen triples

Fix a bucket span `h`. Every valid frozen triple `(s, d, a)` obeys `a - d = h`, so define
`F_h(s, d, a) = (s, d)`. Conversely define `G_h(s, d) = (s, d, d + h)`, where the addition must
succeed and the source shape must contain that intron at those offsets.

- `G_h(F_h(s, d, a)) = (s, d, d + h) = (s, d, a)` because bucket membership gives `a = d + h`.
- `F_h(G_h(s, d)) = (s, d)` directly.

Thus `F_h` and `G_h` are mutual inverses on valid route records. The pair representation is
bijective with the previously specified triple representation within each exact span bucket.
Repeated equal spans at different locations in one shape remain distinct through `donor_offset`;
equal spans in different shapes remain distinct through `local_shape_id`. Sorting and uniqueness
are preserved under the bijection. Candidate membership and both final coordinate checks are
therefore unchanged.

## Threshold and projection guard

All thresholds in `2026-09-01-archive-root-and-jshape-routing.md` remain unchanged. In particular,
the 3.0% route premium, 4.0% complete-collection fraction, 1.25× compressed-shapes ceiling, source
and total-I/O ratios, runtime ratios, RSS limits, exactness requirements, and integrity stop rules
are not relaxed.

The previously cited 33,897,246 compressed `shapes` bytes, 36,588,673 dense-query source bytes,
14,100,036 sparse-query source bytes, predicted route size, and predicted I/O reduction are
**pre-implementation audit measurements or projections**, not Gate-B results. They motivate the
prospective gate but cannot be reported as achieved route performance. Only results from the
committed 96-row panel and the subsequently locked validator may support a Gate-B claim.
