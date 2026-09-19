# Locked gate: packed, cell-sharded EM

**Locked:** 2026-08-31, before implementation or measurement.  The v1 implementation and its
published outputs remain recoverable at tag `gravlax-v1-8-30-2026`; `aie em --eager` is retained
as the within-binary semantic reference during this gate.

## Motivation and design constraints

The optimized D1 EM takes about 64 s but peaks at 38.1 GB.  The peak is not inherent to EM: the
implementation simultaneously retains the eager molecule table, a flat row table, a heap-owning
candidate vector per retained row, a hash table of heap-owning class support vectors, and a heap
vector per target's responsibilities.

The replacement follows the packed equivalence-class layout used by piscem-infer and Salmon:

- candidate labels are flat `u32` arrays with a CSR offset array;
- per-row support during construction is two packed `u32` words (class and gene/evidence flags),
  with no per-row allocation;
- cells are assigned deterministically to 64 shards; each shard is sorted and reduced into a
  packed class/target representation;
- per-cell expression state is a sorted packed key/value table, not a hash table;
- responsibilities are computed into bounded shard-local storage and are not retained as one
  heap allocation per target;
- global reductions occur in fixed shard order, so repeated runs are deterministic.

The masked holdout is selected by a keyed hash of `(seed, cell, class)`, rather than by iteration
order of a hash table.  This deliberately makes the experimental sample independent of allocator
and hash-table layout.  It is a post-v1 protocol improvement, so masked scores need statistical
agreement with the v1 reference rather than byte identity.

## Inputs and arms

- D0: `runs/archive/d0.aie`, GENCODE v49.
- D1: `runs/archive/d1.aie`, GENCODE v49.
- Packed arm: default `aie em`.
- Reference arm: `aie em --eager` from the same binary.
- Impact arm: `--mask 0 --emit`, compared at the emitted MatrixMarket layer.
- Masked arm: `--mask 0.2 --seed 7 --alpha 20`.
- 24 Rayon threads; warm page cache; `/usr/bin/time -v`; three repetitions per timed arm.

## Gates

Correctness:

1. For `--mask 0`, class partition totals (single, mixed, multi-only), recovered target count,
   confident fraction, and emitted nonzero coordinates must equal the eager reference.
2. Emitted D0 and D1 `em.mtx`, `features.tsv`, and `barcodes.tsv` must be byte-identical to the
   v1 frozen outputs.  If float accumulation order alone changes a rounded entry, the gate may be
   relaxed only after reporting maximum absolute error, maximum relative error, relative L1,
   moved mass, and the number of changed rounded entries; scientific equivalence requires
   relative L1 <= 1e-10 and maximum absolute error <= 1e-8.
3. The keyed masked arm must be repeatable byte-for-byte.  Its best-mode top-1 and expected
   accuracy must be within 0.5 percentage points of the v1 eager experiment, and all original
   EM-1 PASS thresholds remain binding.
4. Invalid masks (outside `[0,1]` or non-finite) and invalid alpha (negative or non-finite) fail
   before archive decoding.

Memory and runtime:

- D0 packed peak RSS <= 4 GiB.
- D1 packed peak RSS <= 12 GiB (primary binding gate; >=3.17x below the 38.1 GB reference).
- Median packed wall time <= 1.5x the optimized eager reference on each dataset.
- No output-bearing result may depend on worker count or scheduling.

Interpretation:

- PASS supports describing EM memory as bounded by packed evidence plus shard-local working state.
- A correctness PASS with D1 RSS between 12 and 18 GiB is MARGINAL and triggers an out-of-core
  support-spill follow-up.
- Any scientific-equivalence failure is STOP regardless of memory or speed.

