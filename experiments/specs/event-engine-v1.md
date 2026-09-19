# Event engine v1: frozen semantic contract

Locked 2026-08-31 before implementation and before inspecting event-level D4 output.

## Event catalogue

The catalogue is derived only from junction coordinates and catalogue support in the
annotation-independent archive. It does not use a transcript annotation to decide which
events exist.

- `alt_acceptor`: two observed junctions share a donor. The lower acceptor is the
  deterministic include side and the higher acceptor is the exclude side.
- `alt_donor`: two observed junctions share an acceptor. The lower donor is the
  deterministic include side and the higher donor is the exclude side.
- `cassette`: observed junctions `(d1,a2)`, `(d2,a3)`, and `(d1,a3)` satisfy
  `d1 < a2 <= d2 < a3`. The two flanks form the include side and the skipping junction
  forms the exclude side.

Alternative-site names are genomic, not transcript-directional. An optional GTF or AIC may
label the resulting coordinate-defined events, but may not create or remove candidates.
Event IDs, ordering, and side assignment are deterministic. Every required component must
meet `--min-support`. A hard `--max-events` limit must fail rather than silently truncate.

## Reduction

All component junctions for all events in a request are planned as one union. Each selected
archive chunk is decoded at most once. Within each event, every molecule class is reduced to
exactly one of `include_only`, `exclude_only`, or `both`, even if the class contains several
include junctions or is observed in several chunks. `informative = include_only +
exclude_only`, and `usage = include_only / informative`; `both` is reported but excluded
from the conservative denominator.

Cell and group scopes use the already frozen strict scope grammar. Bulk, group, and cell
rows are exact reductions of the same class-level result. Events below `--min-informative`
after scope selection are omitted.

## Cohort federation

The cohort surface accepts named `ID=ARCHIVE` samples and optional `ID=GROUPS.tsv` group maps.
It forms a coordinate-keyed union of per-archive catalogues, retains events present in at
least `--min-samples` catalogues, and executes each archive locally with the same exact
reducer. Missing components yield `present=false`, not an error or imputed zero evidence.
Outputs are ordered by event ID and caller-supplied sample order. No p-value is reported in
v1; cross-sample comparisons are descriptive because sample, donor, chemistry, and depth are
not interchangeable replicates.

## Annotation labels

Optional annotation labels report overlapping gene identifiers/names, consistent strand
when identifiable, and whether all component junctions are annotated. Labels cannot alter
counts or catalogue membership.

## Compatibility and failure semantics

Existing query and federation output bytes are compatibility gates. Malformed loci,
duplicate sample IDs, unmatched group IDs, unknown barcodes, invalid thresholds,
combinatorial catalogue overflow, and incompatible archive references must fail clearly.
JSON schemas are versioned as `gravlax.query.events.v1` and `gravlax.cohort.events.v1`.

