# Post-v1 junction enumeration and machine-readable output

**Date:** 2026-08-31  
**Protocol lock:** `52dca7a87942db4b44a7d5aea8feeb2c679c7511`  
**Implementation:** `0efabed0d3e92a72c432b4c815c2e0b70c7155e4`  
**Verdict:** correctness, latency, memory, and no-format-change gates **PASS**

## Delivered surface

`aie query ARCHIVE junctions chrom:start-end` now enumerates junctions without requiring the user
to know an exact coordinate first. The default is an explicit 0-based half-open window containing
both endpoints; `--either` selects alternative donors/acceptors crossing a window boundary.
`--min-support` filters the index's supporting-child count. The output deliberately calls that
quantity `supporting_children`: it is neither a read count nor a UMI count.

`--with-cells` decodes the union of selected posting chunks once and adds exact, class-deduplicated
UMI and cell counts; JSON also carries every barcode/count row. `--min-cells` implies that path.
`--gtf` marks exact and individual endpoint membership. The v1 junction catalogue has no strand
bit, so annotation flags are explicitly either-strand. No archive byte or version changed.

The old point query now has clean `--tsv` and `--json` modes, and `--top 0` returns every cell.
Machine stdout contains no summary/timing text; diagnostics go to stderr.

## D0 correctness

The fixed CD44 interval `chr11:35138870-35232402`, filtered at 20 supporting children, returns:

| junction | index children | exact UMIs | cells | v49 annotated |
|---|---:|---:|---:|---:|
| chr11:35186900-35189834 | 26 | 17 | 15 | yes |
| chr11:35221732-35229128 | 77 | 73 | 66 | yes |

For each row, a separate `query junction --top 0 --json` matches the interval result in
coordinates, index support, posting chunks, total UMIs, cells, and the complete ordered barcode
count list. TSV and JSON index rows agree; five output repetitions are byte-identical. Fixed
tests exercise the support and cell filters plus the contained/`--either` boundary distinction.
The full workspace now passes 93 Rust tests.

## Performance

Five warm 24-thread repetitions:

| path | median wall | maximum RSS | gate |
|---|---:|---:|---:|
| catalogue/postings metadata only | **0.01 s** | **28,388 KB** | ≤0.20 s / ≤0.50 GiB — PASS |
| union-postings decode + exact cells | **0.06 s** | **108,988 KB** | ≤2.0 s / ≤2 GiB — PASS |

Supplying the 3.3 GB uncompressed v49 GTF is intentionally outside the index-only timing; it adds
about 1.2 s on this host. The archive's catalogue and postings streams are currently parsed from
the beginning because v1 stores no per-junction byte offsets. At 0.01 s on D0 this is not a reason
for a format change.

## Scientific and product value

This closes the largest credibility gap in the “queryable molecular index” story: a user can now
ask what junctions exist in a locus, filter them, annotate them at query time, and recover exact
cell-level evidence. It also supplies the structure needed for the next high-value features:
group/cell scoping, per-cell junction-set/PSI queries, and overlap-aware discovery models.

## Reproduction

```sh
AIE_BIN=/path/to/gravlax/target/release/aie \
  bash scripts/138_query_junctions_gate.sh
```

