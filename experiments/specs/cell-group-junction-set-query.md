# Cell/group-scoped junction-set query

**Locked 2026-08-31 before implementation or candidate measurement.** The reference is Gravlax
`2690905c6d4d5d47fbd3949692ee8c90e46465e6`. Existing unscoped output is a compatibility
contract.

## Shared scope

`query region`, `junction`, `junctions`, `batch`, and the new `jset` accept mutually exclusive
`--cells` and `--groups` files plus `--agg auto|cell|group|bulk`. Cell files are headerless,
one archive barcode per line. Group files are headerless `barcode<TAB>group`; group order is
first occurrence. `auto` means cell rows without a group file and group rows with one. A group
file is also a cell filter. Parsing is strict: malformed rows, duplicate barcodes, unknown
archive barcodes, empty scopes, and group aggregation without groups fail before chunk decode.

Default invocations must retain their current output schemas and values. Scoped totals count
only selected archive cells. Group totals are sums of already class-deduplicated cell totals;
classes are never deduplicated across cells. Every group is represented, including zero-support
groups, in deterministic first-occurrence order.

## Junction-set semantics

`aie query A jset --include chr:d-a [--include ...] --exclude chr:d-a [--exclude ...]` asks one
class-level set question over the union of all catalogue postings. Within a cell, a UMI class is
assigned to exactly one category:

- `include_only`: supports at least one inclusion junction and no exclusion junction;
- `exclude_only`: supports at least one exclusion junction and no inclusion junction;
- `both`: supports both sides and is not forced into either numerator.

`informative_umis = include_only + exclude_only`, and `usage_fraction = include_only /
informative_umis`. A zero informative denominator is JSON `null` and TSV `NA`. The `both`
category is always reported and excluded from the denominator. This is called junction-set usage,
not transcript PSI: a PSI label is licensed only for a separately specified event geometry and
coverage model.

Duplicate loci within a side and loci appearing on both sides are errors. Catalogue-absent loci
are reported as `present: false` without aborting. Selected chunks are decoded once; each class's
side bitmask is OR-reduced across molecules and physical chunks before its archive cell is read.

## Gates

The synthetic CLI fixture must exercise union deduplication, `both`, null denominator, absence,
and strict scope errors. On the fixed D0 pair, the independent inclusion point count must equal
`include_only + both` for every cell, and the exclusion count must equal `exclude_only + both`.
Scoped cells, grouped rows, and bulk totals must conserve exactly on D0 and D1. Seven JSON runs
must be byte-identical.

At 24 threads, seven D0 runs compare one jset process with two standalone point-junction
processes. Median jset wall time must be at most 75% of the two-process arm; grouped jset at most
1.5 times the unscoped jset, with peak RSS at most 1 GiB. These are execution gates, not a claim
that the registered pair is a biological splicing event.
