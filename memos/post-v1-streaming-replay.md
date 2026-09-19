# Bounded streaming replay verdict

**Date:** 2026-08-30 · **Datasets:** D0 PBMC 1k and D1 PBMC 5k · **Threads:** 24 pinned

## Verdict

**PASS.** Ordinary archive-backed Gene replay can stream bounded molecule batches without making
UMI aggregation local or approximate. Across a balanced five-block same-binary crossover, all 20
timed outputs were byte-identical to the frozen matrices, feature tables, and barcode tables.

| dataset | streaming wall | eager wall | wall ratio | streaming RSS | eager RSS | RSS reduction |
|---|---:|---:|---:|---:|---:|---:|
| D0 PBMC 1k | 2.77 s | 2.38 s | 1.164× | 1.36 GB | 4.22 GB | 3.10× |
| D1 PBMC 5k | 6.90 s | 5.57 s | 1.239× | 5.32 GB | 17.62 GB | 3.31× |

RSS values are median `/usr/bin/time` maximum RSS in decimal GB; the largest individual streaming
observations were 1.37 GiB on D0 and 5.15 GiB on D1. Thus every replicate clears the preregistered
<2 GiB and <6 GiB thresholds, while both median wall ratios clear the ≤1.25× threshold.

## What changed

- The eager `read_archive` path is retained behind `--eager` as a performance-faithful reference.
- Default Gene replay decodes at most two archive chunks per Rayon worker, splits density-skewed
  genomic chunks into fine reducer tasks, and releases each decoded batch immediately.
- Compact `(cell, UMI class, gene, weight)` tuples accumulate in cell shards. MultiGeneUMI_CR,
  edge-based 1MM collapse, and matrix emission still run globally after the scan.
- GTF compilation is also bounded in 192 MiB line-aligned batches. The former whole-file compiler
  alone made the 3.1 GB GENCODE v49 file an impossible floor for the D0 RSS gate.
- The CLI golden fixture now puts one UMI class on both sides of a physical archive boundary and
  checks streaming, eager, BAM, and post-correction BAM replays for identical bytes.

The initially attractive implementation was rejected during tuning: it appeared to pass the
runtime gate only because its `--eager` comparator had itself been slowed by artificial reducer
batching. Restoring the original parallel eager shard assembly exposed a 1.75× D1 penalty. Finer
work stealing and a bounded two-unit decode window produced the accepted result above.

## Interpretation

This is a memory-scalability result, not a new quantification claim. It changes the default Gene
replay execution path while preserving the exact archive-replay contract. Velocity, EM, audit,
from-BAM, and post-correction-BAM modes still use eager molecule materialization and should not be
described as streaming.

The next binding generalization gate is mouse 5′ data and a chemistry manifest; controlled depth
scaling remains deferred because it is less likely to change the paper than testing a new species,
strand geometry, UMI length, and assay end.
