# Disk-backed compact-support shards

**Locked:** 2026-08-31 after batch compaction failed the 4-GiB cap and before disk-shard code or
measurement. This remains an engineering remediation; the D3 statistical result is not rerun or
reinterpreted.

For archives whose `n_molecules * sizeof(EmSupport)` exceeds 512 MiB, write batch-compacted
supports as two little-endian `u32` words into one temporary stream per deterministic cell shard.
After the single archive scan, load, sort, reduce, and release one shard at a time. Small archives
retain the current in-memory path. Temporary storage lives below the standard temporary directory
and must be removed on success and error. The candidate order and all EM floating-point operations
remain unchanged.

The exact D3 shuffled-group seed-7 `D=8,p=8` hybrid-only run must reproduce stdout and JSON
byte-for-byte, use at most 4,194,304 KiB RSS, finish within 60 seconds, and leave no temporary
artifacts. All Rust tests and a D0 in-memory byte-identity regression must pass. If it fails, the
next choices are an external sort or an explicitly measured larger operational envelope.

