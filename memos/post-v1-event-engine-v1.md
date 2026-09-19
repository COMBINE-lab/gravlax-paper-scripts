# Coordinate-defined event and cohort engine

## Verdict

The query capability passes; the locked biological vignette does not. Gravlax now discovers
alternative-acceptor, alternative-donor, and cassette structures from archive junction
coordinates/support alone, then reduces all event sides in one decode plan. Named cohort shards
return exact local bulk or group reductions without treating cells or samples as replicates.

## Exactness and scale

- The synthetic integration fixture exercises every event type, multi-flank deduplication,
  include-only/exclude-only/both categories, annotation-neutral labels, cohort identity, strict
  input rejection, and hard catalogue overflow.
- All 115 Rust tests pass.
- The first 50 retained D0 chromosome-1 events equal standalone `jset` totals and every group
  row. All 2,198 cohort rows that clear the independent per-sample threshold equal separate
  `query events` results. Repeated D0 JSON is byte-identical.
- Outputs from the pre-feature binary and implementation commit are byte-identical for region,
  point junction, junction enumeration, junction set, and batch JSON; federate bytes are equal
  after normalizing its printed elapsed time.
- D0 chromosome 1: 1,905 candidates, 1,211 retained at ten informative molecules, 4,263
  independent chunk decodes reduced to 57, 0.201 s median and 336,472 KiB peak RSS.
- D1 chromosome 1: 17,951 candidates, 9,242 retained, 39,526 decodes reduced to 58, 1.64 s and
  1,469,384 KiB.
- D0/D1 shared cohort: 1,311 candidates, 1,291 retained, 0.80 s and 1,171,404 KiB.
- Returning all 1,211 D0 events takes 6.86% of the median wall time required to run only the
  first 50 events as standalone `jset` processes.

## Biological gate and correction

The earlier post-hoc D0 AKR1A1 junction set was not a complete cassette: one inclusion flank is
absent from the D0 catalogue. The event engine correctly does not infer an unobserved component
from annotation. The locked event is complete in untouched D4 and has 110 informative
molecules, passing the support threshold. Only g3 and g4 have at least 20 informative molecules;
their usages are 0.86957 and 0.82759, a 0.04198 range below the frozen 0.05 threshold. This opens
the predeclared `FAIL_held_out_only` branch: publish the query engine, not a biological
generalization claim.

## The next ambitious step

The strongest successor is a replicate-aware molecular splice-graph layer, not another local
event spelling. Build gene/locus graphs whose edges are archive junctions and whose hyperedges
are molecule-class path fragments, then query mutually exclusive paths, multi-exon co-occurrence,
terminal-exon choice, and intron retention across a sample/design table. A donor-aware
beta-binomial or Dirichlet-multinomial contrast should sit above exact per-shard counts and require
genuine biological replicates; cells must not become pseudoreplicates. The binding experiment is
an external multi-donor PBMC dataset with a frozen design, truth/replication denominator, and
comparison to established differential-splicing methods. This would turn the current local-event
engine into a defensible population splicing atlas rather than merely adding more query syntax.

