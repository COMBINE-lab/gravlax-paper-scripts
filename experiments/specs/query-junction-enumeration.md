# Junction enumeration and machine-readable query gate

**Locked:** 2026-08-31, before implementation or full-archive timing. **Scope:** add the first
query-surface sprint without changing the archive format.

## User-facing contract

Add:

```text
aie query ARCHIVE junctions chrom:start-end \
  [--either] [--min-support N] [--with-cells] [--min-cells N] \
  [--gtf annotation.gtf] [--tsv|--json]
```

- Coordinates are documented as 0-based half-open. By default both donor and acceptor must lie
  within `[start,end)`; `--either` selects a junction when either endpoint lies in the window.
- Rows are sorted `(chrom, donor, acceptor)` and contain `supporting_children` and
  `posting_chunks`. These are index-construction quantities and must not be mislabeled reads or
  UMIs.
- `--with-cells` additionally decodes only the selected postings and reports exact
  class-deduplicated `umis` and `cells`, using the existing point-junction semantics.
  `--min-cells` implies this path.
- `--gtf` adds `annotated`, `donor_annotated`, and `acceptor_annotated`; because the v1 junction
  catalogue has no strand bit, these mean present on either annotation strand.
- `--tsv` emits a header and rows only to stdout; `--json` emits one valid JSON object only.
  Diagnostics/timing go to stderr. The two flags conflict.

Also add `--tsv` and `--json` to the existing point `junction` query. `--top 0` means all cell
rows. The legacy no-format output remains unchanged.

## Correctness gates

1. Unit fixtures cover catalogue delta decoding, interval boundary behavior, postings decoding,
   support filtering, and annotation endpoint flags.
2. On D0, every row from an interval `--with-cells --json` query equals a separate point
   `junction --json` in support, UMI, cell, and per-cell counts.
3. TSV and JSON carry the same ordered rows and parse without log stripping; repeated invocations
   are byte-identical after excluding stderr timings.
4. `--min-support`, `--min-cells`, default contained semantics, and `--either` are each exercised
   by an automated gate.
5. The complete Rust workspace test suite passes.

## Performance and memory gates

Benchmark D0 at 24 threads over a fixed expressed locus selected from the existing catalogue:

- index-only TSV enumeration, five warm repetitions: median wall <=0.20 s and peak RSS <=0.50 GB;
- `--with-cells`, five warm repetitions: median wall <=2.0 s and peak RSS <=2.0 GB;
- index-only implementation performs no chunk reads or molecule decode;
- no archive-format bytes are added.

Failure narrows the release to correctness-only and triggers a postings-offset design study; it
does not justify a v2 format change without profiling evidence.

## Outputs

- code tests plus `scripts/138_query_junctions_gate.sh`;
- compact result `results/post-v1-query-junctions.json`;
- memo `memos/post-v1-query-junctions.md`.

