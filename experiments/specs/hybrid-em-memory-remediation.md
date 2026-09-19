# D3 hybrid-EM memory remediation

**Locked:** 2026-08-31 after the D3 formal failure and before any remediation code change or
measurement. This is an engineering gate, not a second statistical confirmation.

The D3 hybrid passed every scientific criterion but used 7.30--7.45 GiB peak RSS in hybrid-only
runs. Code inspection identifies a concrete transient: all 64 support shards are finalized in
parallel, so their sorting and temporary unique/key event arrays coexist even though the final
EM executes shards in fixed order. The first remediation bounds shard finalization concurrency
without changing the packed persistent layout, estimator, candidate order, or arithmetic order.

Run the exact D3 shuffled-group seed-7 hybrid-only command at `D=8,p=8`. Its `metrics.json` and
stdout must be byte-identical to the archived reference, all Rust tests must pass, peak RSS must
be at most 4,194,304 KiB, and wall time must be at most 60 seconds (3.5x the 17.16-s reference).
Passing shows that the operational cap has been remediated, but the historical D3 gate remains a
formal FAIL. Failure triggers true out-of-core target-shard storage rather than another model or
threshold change.

## Result

**FAIL.** Serial finalization preserved `metrics.json` byte-for-byte but used 7,639,532 KiB peak
RSS and 28.69 s. This rejects the proposed root cause: the peak precedes final shard sorting and
arises while raw per-representative support records accumulate. The serial patch is discarded.
A separately locked exact batch-compaction gate precedes a disk-backed implementation because it
can eliminate the same redundant support records before they become resident.
