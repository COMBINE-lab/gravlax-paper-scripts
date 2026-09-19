# Bounded streaming replay gate

## Hypothesis

An archive reader that decodes bounded genomic chunk batches and immediately reduces them into
global cell shards can preserve exact Gene replay while cutting peak RSS below 2 GB on D0 and
6 GB on D1. Because final aggregation remains global by `(cell, UMI class, gene)`, UMI classes
spanning archive chunks must behave identically to eager replay. Median wall time must be no more
than 1.25 times the same-binary eager reference.

## Implementation and controls

- `aie replay-rows` streams by default for ordinary `.aie` Gene replay; `--eager` retains the
  full-materialization reference arm.
- Independently compressed chunks are decoded in batches of at most two chunks per Rayon worker.
  Each decoded chunk is subdivided into uniform reducer tasks to avoid genomic-density imbalance.
- Only compact `(cell, class, gene, weight)` assignments survive a batch. The existing global
  cell-shard sort, MultiGeneUMI_CR choice, UMI-edge collapse, and deterministic matrix emitter are
  unchanged.
- The GTF compiler uses bounded 192 MiB, line-aligned batches; otherwise the 3.1 GB uncompressed
  GENCODE v49 file alone makes the D0 RSS threshold impossible.
- `scripts/121_streaming_replay_gate.sh` runs a five-block alternating same-binary crossover on a
  fixed CPU set. Every matrix, feature table, and barcode table is byte-compared with the frozen
  reference. `scripts/122_summarize_streaming_replay.py` creates the compact result record.

## Datasets and gates

- D0: 10x PBMC 1k 3′ v3, GENCODE v49; peak RSS <2 GB.
- D1: 10x PBMC 5k 3′ v3, GENCODE v49; peak RSS <6 GB.
- Both: five exact streaming and five exact eager outputs; median streaming wall / eager wall
  ≤1.25.

Failure of exactness rejects the implementation. Passing memory but failing time requires either
better within-batch load balancing or a documented fallback to eager mode; it does not support a
default-mode change.
