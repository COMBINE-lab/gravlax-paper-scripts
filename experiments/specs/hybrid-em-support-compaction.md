# Exact batch-local support compaction

**Locked:** 2026-08-31 after serial finalization failed and before this code change or its memory
measurement. D3 scores are already known; this is an engineering test only.

The streamed builder currently appends one 8-byte support for every classified representative,
although finalization consumes only whether each `(class,gene)` received a unique and/or
multimapper flag. Within each bounded decode batch and cell shard, sort supports by class and gene
and replace duplicate records by one record whose flags are their bitwise OR. Retain the final
global sort/reduction to combine keys spanning batches, and remove the `n_molecules`-sized support
preallocation. This is algebraically exact and changes neither candidates nor EM arithmetic.

The exact D3 shuffled-group seed-7 hybrid-only command at `D=8,p=8` must reproduce the reference
stdout and JSON byte-for-byte, pass all Rust tests, use at most 4,194,304 KiB peak RSS, and finish
within 60 seconds. Failure triggers disk-backed support shards.

## Result

**FAIL on the 4-GiB cap; retained as an exact partial improvement.** Commit `7c1ef42` reproduced
stdout and JSON byte-for-byte, reduced peak RSS from 7,653,508 to 6,562,508 KiB (14.24%), and ran
in 16.71 s versus 17.16 s. Redundant representatives are material but do not explain most of the
peak. The registered next step is disk-backed compact-support shards.

