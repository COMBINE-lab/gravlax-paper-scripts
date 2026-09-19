# Atlas collection and evidence-v2 campaign

## Ordered objective

This campaign begins only after the protocol-aware PolyASite/NTRK2 result. Its first objective is
an optional federated collection index, followed by bounded-memory cohort execution. Archive-v2
payload changes are pilots and may be rejected without weakening the collection work. Parallel
implementation is allowed, but promotion and manuscript claims follow this order.

## Collection sidecar

The collection artifact is derived state. Every `.aie` remains independently readable and is the
only authoritative molecular evidence. The sidecar records exact archive identities and reference
signatures, then merges only existing index metadata: genomic chunk extents, junction coordinates,
per-junction total-support upper bounds, and local posting routes. Building it must not read any
`cN` molecule section. A global coordinate maps first to samples, then to each sample's existing
local chunks. Sample absence is represented explicitly by the global sample dimension.

The first production surface covers exact region, junction, and junction-set federation. A query
plan reports candidate samples, pruned samples, local chunks, support-bound pruning, estimated and
actual bytes read, and phase timings. Total support may prove that a row cannot pass a gate; it may
never substitute for group- or cell-resolved counting. Group maps and annotations remain query-time
inputs. A cache derived from either must carry their full digests and live outside `.aie`.

Construction is deterministic and segment-friendly. Existing immutable segments are not rewritten
when archives are appended. A compacted clean build and the logical union of segments must expose
the same sorted catalogue. Duplicate archive identities, duplicate sample identifiers, changed
archive bytes, mixed reference signatures, invalid local postings, integer overflow, and truncated
segments are hard failures.

Acceptance compares the indexed path to the current naive federation on adversarial synthetic
archives and frozen real queries, including absence, repeated representatives, cross-chunk UMI
classes, reverse strand, and incompatible references. Two identical builds must hash identically.
The sidecar must remain under 5% of aggregate archive size. The six-archive targeted query may not
regress by more than 10%; sparse queries must demonstrate that irrelevant archives and chunks are
not opened or decoded.

## Bounded cohort execution and output

The SEZ PolyASite mixture is the binding workload because it combines 45,705 recurrent sites,
35.7 million expected UMIs, leave-one-donor-out kernels, and per-gene EM. The implementation may
use a two-pass scan or an ignored temporary spool. It must preserve the fitted kernel, floating-
point operation order within every site mixture, recurrence decisions, paired inference, and final
row order. The three numerical TSVs must remain byte-identical.

Peak memory must be bounded by an active chromosome or gene shard and fall from 3.49 GiB to at
most 2.25 GiB on the frozen eight-donor run, while wall time remains at most 30 seconds. Large
machine output uses a versioned sparse fact table plus explicit sample/group/site dimensions;
missing fact rows are logical zeros, never missing observations. TSV/JSON remain available for
small outputs. A columnar representation is promoted only if it streams, round-trips exactly, and
does not impose disproportionate core binary or build cost. The existing pushed-down federated
event result is the output-size gate and must shrink at least fourfold without dropping fields.

## Adaptive evidence frontier

The pilot generalizes the fixed most-contained/most-extended pair only where additional read
geometries are nondominated under an annotation-independent containment relation. Ordering and tie
rules are deterministic; exact chain multiplicity is retained. The current representation remains
the default and must be byte-identical when the option is absent.

The pilot reports representatives per chain, extra bytes by section, replay discrepancy against
the full BAM path, and spliced/unspliced/ambiguous velocity discrepancy, stratified by chains whose
middle reads were previously omitted. Promotion requires at most 5% archive growth and either a
20% relative improvement in the principal incomplete-chain velocity discrepancy or a 25% relative
improvement in Gene replay discrepancy, with no component worsening by more than 0.25 percentage
point. A written result must specify the exact assignment-predicate class preserved by the
frontier and give counterexamples for policies outside that class.

## Sparse terminal-tail evidence

The pilot inspects oriented terminal soft clips during ingest and stores only positive compact
tail evidence: genomic anchor, direction, capped untemplated-tail length, and a confidence tier.
It never stores arbitrary bases or qualities. Disabled mode must leave archive bytes and replay
unchanged. No-tail is missing evidence, not a negative call.

Synthetic reads fix forward/reverse, left/right clipping, short/long tail, mismatch, and wholly
templated cases. Real evaluation is performed independently in at least three 10x 3′ datasets.
Tail-positive molecules are compared with matched-depth and matched-terminal-region controls,
external PolyASite candidates, and genomic internal-priming candidates. Promotion requires a
maximum 2% archive premium (target below 0.5 bit/read), enrichment in at least three datasets, and
positive held-out discrimination or calibration gain in at least three. Otherwise the prototype
is retained only as an audit result.

## Interpretation

Failure of a v2 pilot is useful evidence against paying a permanent storage cost. Collection
speedups never imply biological validity. Tail support never converts all fragment endpoints into
cleavage calls. No result may relax a frozen gate after outcome inspection.
