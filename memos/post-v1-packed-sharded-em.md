# Packed, cell-sharded EM: PASS

**Date:** 2026-08-31  
**Code:** Gravlax `3d7dcd5d739fc9d1bbee0a889569ff30d1e1f913`  
**Locked protocol:** `experiments/specs/packed-sharded-em.md`  
**Results:** `results/post-v1-packed-sharded-em-d0.json`,
`results/post-v1-packed-sharded-em-d1.json`

## Outcome

The high EM peak was representational, not intrinsic to inference. The v1 path retained the
eager molecule table, flattened rows, a candidate `SmallVec` per retained row, a hash table of
heap-owning class vectors, and one heap vector of responsibilities per target. The new default
streams archive batches into 8-byte `(class, gene/evidence-flags)` records, reduces 64
deterministic cell shards, and runs over flat `u32` labels with `u64` CSR offsets and packed
cell--gene state. `aie em --eager` retains the v1 implementation as a reference.

| dataset / arm | packed median wall | packed median RSS | v1 wall | v1 RSS | speedup | RSS reduction |
|---|---:|---:|---:|---:|---:|---:|
| D0 masked | 3.10 s | 1.221 GB | 12.93 s | 8.14 GB | 4.17x | |
| D0 impact + emit | 3.64 s | 1.243 GB | 13.27 s | 8.14 GB | 3.65x | **6.55x** worst-arm |
| D1 masked | 10.06 s | 4.559 GB | 64.1 s | 38.09 GB | 6.37x | |
| D1 impact + emit | 10.01 s | 4.506 GB | 73.3 s | 38.09 GB | 7.32x | **8.35x** worst-arm |

All figures are medians of three 24-thread runs. D0 and D1 clear the pre-registered 4/12 GiB
ceilings and the 1.5x-runtime ceiling by wide margins. Three repeated masked outputs per dataset
are byte-identical.

## Scientific equivalence

Class partitions are unchanged: D0 has 10,084,709 single, 886,655 mixed, and 688,558 multi-only
classes; D1 has 46,853,364 / 4,600,483 / 3,084,135. The keyed holdout deliberately selects a
different, hash-layout-independent evaluation sample. Pooled/blend top-1 remains 98.2/98.3% on
D0 and 98.3/98.3% on D1; expected accuracy is 97.0--97.7%, within 0.1 percentage point of v1 and
well above every original EM-1 PASS threshold.

At four-decimal MatrixMarket precision every impact-layer value and total mass equals the frozen
v1 reference: relative L1, maximum absolute error, and moved mass are all zero. Packed reduction
places nine D0 and three D1 sub-rounding values just above the historical `1e-6` storage cutoff,
so those files contain that many additional explicit `0.0000` coordinates. There are no added or
missing nonzero coordinates. We report this rather than tuning a data-dependent epsilon to force
byte identity.

## Reference implementations inspected

- piscem-infer `7182fce639d3a0994fa99bf621fa76d66e92d286`,
  `src/utils/eq_maps.rs`: flat `eq_labels`, `eq_label_starts`, and counts.
- Salmon rewrite `d3be8b31d31242683c01fc5510c48fbd8bc917fb`,
  `crates/salmon-infer/src/packed.rs`: pre-sized flat labels/weights, `u64` incidence offsets,
  and omission of arrays unused by the chosen inference mode.

The Gravlax adaptation additionally shards by cell and recomputes responsibilities from packed
state instead of retaining a vector per target. An out-of-core spill is not justified: D1 is
already below 4.7 GB in every observed run.

