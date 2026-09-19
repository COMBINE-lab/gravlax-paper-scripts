# Bounded-memory PolyASite cohort execution

## Verdict

PASS. Gravlax commit `655a7de923e71d9eaf0ef5b40f2a0afb3c3f3e55` reproduces the frozen
eight-donor PolyASite result in 28.21 seconds at 2,262,980 KiB maximum RSS. The registered limits
were 30 seconds and 2.25 GiB (2,359,296 KiB). `sites.tsv`, `genes.tsv`, and
`fragment-kernel.tsv` are byte-identical to the frozen primary result.

Relative to the prior 22.02-second hot path, maximum RSS falls from 3,660,468 to 2,262,980 KiB,
a 38.2% reduction, for a 1.28x wall-time cost. The new path remains 1.32x faster than the
37.16-second pre-hot-path implementation. The acceptance binary was built from a clean detached
checkout of the recorded commit; its SHA-256 is
`b7af8ad99d78a260cf74df5c73b211238587a3c9f3a9732dc0a131453302ddb9`.

## Implementation

- Endpoint records are collapsed after each chromosome rather than retained for the whole
  archive panel.
- Coordinate counts use one sorted packed array per sample instead of a hash table followed by a
  duplicate sample-by-gene hierarchy. Gene fits borrow contiguous ranges from that array.
- Decoded molecule batches contain at most eight chunks. Cell-of-class blocks are prefetched,
  read through an immutable lookup during parallel classification, and released between
  chromosomes.
- At most two adjacent archives are active, and only when their combined compressed size is at
  most 256 MiB. Larger archives are reduced alone. This preserves useful concurrency without
  restoring an atlas-wide frontier.
- Non-monotone or out-of-range chromosome chunk indexes fail before reduction.

This is a measured working-set bound, not a claim of constant memory at arbitrary sequencing
depth: the active chromosome's selected molecule classes still scale with evidence in that
chromosome. A disk-backed class shard would be needed for a strict asymptotic bound.

## Reproduction

Build commit `655a7de923e71d9eaf0ef5b40f2a0afb3c3f3e55`, then run
`scripts/203_benchmark_polyasite_mixture_hotpath.sh` against
`runs/post-v1/sez-polyasite-mixture-crossfit-controls-r1/primary-gap24`. Summarize fail-closed with
`scripts/204_summarize_bounded_polyasite.py`. The compact result is
`results/post-v1-polyasite-bounded-memory.json`; large tables and timing records remain ignored.
