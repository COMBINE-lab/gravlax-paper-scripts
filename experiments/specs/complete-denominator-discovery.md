# Complete-denominator annotation-history discovery gate

**Locked:** 2026-08-31, before generating D1/v32 discovery output or computing this gate's
statistics. **Scope:** human 3′ PBMC D0 (primary) with D1 as an independent-library support arm.
This is a prospective annotation-history test, not a general de-novo gene-calling benchmark.

## Question and fixed inputs

Can the v1.5 archive, queried using only GENCODE v32, rediscover genes that first appear in
GENCODE v49, when every sufficiently expressed v49-added gene is included in the recall
denominator—including genes overlapping an older gene model?

- evidence archives: `runs/archive/d0.aie` and `runs/archive/d1.aie`;
- old and future annotations: GENCODE v32 and v49, with version suffixes removed from gene IDs;
- D0/D1 expression: fresh v49 STARsolo `Gene/raw`, summed only over each dataset's fixed v49
  `Gene/filtered/barcodes.tsv` cell set;
- discovery calls: `aie query <archive> discover --gtf gencode.v32.annotation.gtf --tsv`, with
  the command defaults frozen (`merge-gap=1000`, `min-umis=10`);
- orthogonal/supporting observations: D1 replication, PolyASite 2.0, and fresh v49 STAR splice
  junctions. GENCODE v49 defines annotation-history truth and is not described as orthogonal.

All coordinates are converted to 0-based half-open intervals. A discovery locus matches a truth
gene only when chromosome and strand agree and their intervals overlap. Ties and duplicate calls
are resolved by deterministic maximum-cardinality bipartite matching; permissive any-overlap
recall is reported only as a sensitivity analysis.

## Complete truth denominator

The primary truth set contains **every** v49 gene ID absent from v32 with at least 50 D0 Gene
UMIs in the fixed D0 called-cell set. No overlap filter is allowed. Each truth gene is assigned to
exactly one predeclared stratum:

1. `clean`: overlaps no v32 gene span;
2. `antisense_only`: overlaps one or more v32 genes, all on the opposite strand;
3. `same_strand`: overlaps at least one v32 gene on the same strand (including mixed
   same/opposite-strand overlap).

Report numerator and denominator for every stratum, all unmatched truth genes, and both macro and
UMI-weighted recall. Also report thresholds 10 and 100 UMI as sensitivity analyses, without using
them to replace the primary 50-UMI result.

## Ranked event set and annotation-history precision

The ranked event table consists of every emitted D0 v32-discovery locus plus one zero-score row
for each unmatched eligible truth gene. A matched emitted locus is positive; extra loci—including
duplicate calls that cannot receive a one-to-one truth match—are negative for this narrowly named
**annotation-history confirmation** endpoint. Rank by D0 UMI support, use deterministic coordinate
order for ties, and report average precision, precision/recall at the fixed 10-UMI emission floor,
and precision at top 100/500/1000. This endpoint must not be called biological precision: a locus
unsupported by a later annotation is unresolved, not proven false.

## Independent/support evidence

For every primary truth gene and every emitted D0 candidate, record these flags independently of
the v49 match label:

- `D1_expression`: at least 10 D1 Gene UMIs over D1's fixed called cells;
- `D1_discovery`: same-strand overlap with a D1 v32-discovery locus having at least 10 UMIs;
- `PolyASite`: a same-strand PolyASite 2.0 site lies inside the interval or within 500 bp on its
  transcript-oriented downstream side;
- `STAR_splice`: a fresh D0 v49 STAR junction with at least 3 uniquely mapping reads lies wholly
  inside the interval on a compatible strand.

For candidate replication, interval overlap with a same-strand D1 discovery locus is the primary
test. Report support combinations, not a post-hoc composite truth label. The STAR and PolyASite
flags are expected to be insensitive for unspliced/single-exon genes, so absence is not a failure.

## Gates and interpretation

- accounting: all primary truth genes occur in exactly one stratum and exactly one of
  matched/unmatched; all D0 candidates occur in exactly one annotation-history class;
- clean-stratum recall: **PASS ≥90%**, marginal 75–90%, fail <75%;
- all-strata macro recall: **PASS ≥60%**, marginal 35–60%, fail <35%;
- same-strand-overlap recall is reported with no minimum; it tests the expected conservative
  failure mode rather than silently removing it;
- D0-candidate independent-library recurrence: report locus and UMI-weighted fractions overall
  and by v49-match status; no minimum is imposed because locus-definition sensitivity is itself a
  result;
- rerunning the evaluator must reproduce all non-timing fields byte-for-byte.

A clean-pass/overall-fail result narrows the paper claim to conservative rediscovery in
previously unoccupied loci. No result licenses “complete discovery,” “de novo gene calling,” or
precision without the `annotation-history` qualifier.

## Outputs

- compact summary: `results/post-v1-complete-denominator-discovery.json`;
- eligible truth audit: `results/post-v1-complete-denominator-truth.tsv`;
- D0 candidate audit: `results/post-v1-complete-denominator-candidates.tsv`;
- protocol/evaluator: `scripts/136_complete_denominator_discovery.sh` and
  `scripts/137_complete_denominator_discovery.py`;
- interpretation memo: `memos/post-v1-complete-denominator-discovery.md`.

