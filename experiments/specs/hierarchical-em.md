# Locked gate: query-time hierarchical EM

**Locked:** 2026-08-31, before implementation or measurement.  This is a post-v1 experiment;
the frozen release remains recoverable at tag `gravlax-v1-8-30-2026`.  The packed/sharded EM at
Gravlax `edec25f` is the implementation and performance reference.

## Question and estimand

The experiment asks whether ambiguous molecules are assigned more accurately when expression is
shared primarily among transcriptionally related cells rather than across the entire sample.
The hierarchy is query-time metadata and does not change the archive.  It is evaluated only on
called cells with supplied group labels; unlabelled archive cells retain the current global
fallback and are not included in the head-to-head scores.

For candidate gene `g`, cell `c`, and supplied group `h(c)`, the new model uses

```text
w(c,g) = local(c,g)
       + alpha_group  * group(h(c),g) / total_group(h(c))
       + alpha_global * global(g)     / total_global
       + epsilon,
```

then normalizes `w` over the target class.  A `group` arm uses only the normalized group term.
Existing `uniform`, `cell`, `pooled`, and `blend` arms retain their exact definitions.  Ten EM
iterations, 64 deterministic cell shards, fixed shard reduction order, and epsilon `1e-9` are
unchanged.  Production `--mask 0 --emit` output remains the pooled arm regardless of group input.

## Inputs and leakage control

- D0 PBMC 1k is development only; D1 PBMC 5k and D2-prime brain nuclei are confirmation.
- Archives are `runs/archive/{d0,d1,d2p}.aie`, replayed against GENCODE v49.
- Targets use the existing keyed 20% mixed-class mask with seeds 7, 17, and 29.  Selection depends
  only on `(seed, cell, class)` and the target's unique contribution is removed from fitted counts.
- Candidate genes appearing in any masked target for a dataset are exported before grouping and
  excluded from group construction.  If fewer than 500 eligible genes or fewer than 20 cells in
  any selected group remain, that dataset is STOP rather than silently relaxing the exclusion.
- Groups are constructed from the base Gravlax Gene matrix before any EM layer, restricted to the
  externally fixed called-cell barcode list.  Normalization is library-size 10,000 plus `log1p`;
  replay-derived HVGs (up to 2,000) and PCA are fit on replay alone.  Dimensions and resolutions
  are the already locked robust macro-partitions: D0 40/0.60, D1 30/0.40, D2-prime 40/0.25.
  Twenty Leiden seeds (0--19) are run and the within-arm ARI medoid labels are used.  No fresh
  STARsolo expression values or cluster labels are model inputs.

This mixed-class test is biologically grounded but not a complete causal validation: target
classes are selected from observed mixed evidence.  Results must be described as masked recovery,
not as ground-truth quantification accuracy.  Independent ambiguity corruption remains a distinct
follow-up and cannot be silently pooled with this gate.

## Controls and selection

For each dataset, evaluate:

1. the real replay-derived labels;
2. a deterministic size-preserving permutation of those labels;
3. all archive cells collapsed to one group while retaining the called-cell scoring mask (this
   must agree numerically with pooled sharing on those same scored targets).

On D0 only, evaluate `alpha_group in {5,20,80}` crossed with `alpha_global in {0,5,20}` and choose
the pair minimizing mean hierarchical negative log loss across the three seeds.  Exact ties are
broken by smaller `alpha_group`, then smaller `alpha_global`.  Carry that pair unchanged to D1 and
D2-prime.  The existing blend alpha remains 20.

## Metrics and acceptance

For every arm and seed report target count, top-1 accuracy, mean probability on the truth,
negative log loss, multiclass Brier score, calibration counts, wall time, and peak RSS.  Preserve
raw stdout, `/usr/bin/time -v`, group/candidate digests, command lines, code commit, binary digest,
and a compact JSON summary.  Real and control arms must score identical target identities.

PASS requires all of the following:

- combined D1+D2-prime hierarchical negative log loss improves at least 5% over the better of
  pooled and blend;
- top-1 falls by no more than 0.5 percentage point on either confirmation dataset;
- on D2-prime, hierarchical top-1 improves at least 1 point or log loss improves at least 10%;
- real labels achieve at least half of the hierarchy-over-baseline log-loss gain beyond the
  size-preserving shuffled labels;
- the one-group arm agrees with pooled probabilities to `1e-10` in aggregate metrics;
- median wall time is at most 1.5 times packed EM and peak RSS is below 6 GiB on D1.

MARGINAL means a repeatable proper-score improvement fails one hard-assignment or magnitude gate;
retain it as analysis but not a default.  STOP means confirmation fails, shuffled labels explain
the gain, collapse invariance fails, or a group is harmed enough to violate the aggregate top-1
guardrail.  No outcome changes the default production emission without a separate decision.
