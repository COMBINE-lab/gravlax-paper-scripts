# Next goals after bounded discovery and hot-path optimization

**Date:** 2026-08-31  
**Status:** recommendation draft, not a frozen experiment specification  
**Current code:** Gravlax `edec25f` (documentation) / `c1ef9c4` (last implementation change)

## Executive recommendation

The next paper-changing goal should be **cross-fitted hierarchical EM**, ahead of another storage
codec pass or a larger query surface. Packed, cell-sharded EM removed the engineering blocker:
D1 now runs in 10.06 s at 4.56 GB rather than 64 s at 38.1 GB. The remaining weakness is
scientific, not computational. The current `cell`, `pooled`, and `blend` comparison shows that
sharing helps, especially on sparse nuclei, but `blend` has only cell and whole-sample levels and
the masked truth is constructed from the same mixed-class structure the model consumes.

A group-aware model can answer the more interesting question: *does borrowing primarily from
related cells improve ambiguous-molecule assignment without allowing dominant populations to
wash out rare ones?* A negative result is still useful: it would show that global pooling already
captures nearly all available signal and prevent hierarchy language from entering the paper.

## Why this is now first

1. **Large sparse-cell headroom.** On D2′ brain nuclei, per-cell top-1 is 59.4% while pooled is
   90.9%; D0/D1 are near ceiling. A hierarchy has a plausible place to improve calibration and
   rare-population behavior even if aggregate D0 top-1 barely moves.
2. **No format change.** The archive already stores the needed cell, class, and candidate-gene
   relations. Group labels are query-time metadata.
3. **The packed layout is ready.** Adding a `u32` group per cell/target and compact group-gene
   states is small beside current CSR targets. A dense 10-group by ~80k-gene `f64` table is only
   about 6 MB; sparse group-gene CSR is available if labels become numerous.
4. **It strengthens an existing central claim.** This tests when cross-cell sharing helps rather
   than adding another isolated command.

## Proposed model surface

Keep the current modes unchanged and add two evaluation modes when `--groups barcode-group.tsv`
is supplied:

- **group:** responsibility proportional to the expression estimate within the supplied group;
- **hierarchical:** per-cell evidence plus fixed pseudo-count mass from the group distribution
  plus a smaller fixed mass from the sample-wide distribution.

For candidate gene `g` in cell `c`, group `h`, use the interpretable shrinkage form

```text
w(c,g) = local(c,g)
       + alpha_group  * group(h,g)  / total_group(h)
       + alpha_global * global(g)   / total_global
       + epsilon
```

and normalize `w` over the target's candidate set. This extends the existing `blend` model rather
than replacing its semantics. Production emission must remain pooled until hierarchy passes the
evaluation; the first implementation is evaluation-only.

Implementation changes are bounded:

- parse and validate a two-column barcode/group map against the archive cell dictionary;
- retain target class IDs (one additional `u32` per target) for folds and audit;
- accumulate deterministic group-gene bases and iterates alongside cell/global arrays;
- add log loss, Brier score, and group/depth strata to the machine-readable evaluator;
- preserve the packed streaming builder, 64 cell shards, fixed iteration order, and `--eager`
  reference.

## Evaluation must be stronger than the current mask

The hierarchy should not be judged only by the existing 20% mixed-class mask. Use two target
constructions and report them separately:

1. **Cross-fitted mixed classes:** retain the present biologically grounded targets, but assign
   folds by a keyed `(seed, cell, class)` hash, remove each target's unique evidence before fitting,
   and choose shrinkage on development folds only.
2. **Independent ambiguity corruption:** select single-gene classes by keyed hash, remove their
   unique count, and attach a candidate set sampled deterministically from a different real
   ambiguous class containing the truth gene. This keeps real paralog-set geometry while making
   target selection independent of the fitted ambiguity evidence. It is a stress test, not a
   claim that the synthetic class occurred in the reads.

Group labels must come from the **base Gravlax Gene matrix, before any EM layer**, using the
already locked stable macro-partition settings: D0 40 PCs/resolution 0.60 (7 groups), D1
30/0.40 (9), and D2′ 40/0.25 (9). To avoid direct target-gene leakage, recompute these partitions
after excluding genes appearing in evaluation candidate sets, and record label digests. Do not
use fresh-STARsolo cluster labels as model inputs.

Required controls:

- deterministically permute group labels while preserving group sizes;
- collapse all labels to one group (must reduce to the global model within numerical tolerance);
- stratify by cell UMI quartile, group size, candidate-set size, and truth-gene abundance;
- report D0 development separately from untouched D1 and D2′ confirmation;
- repeat keyed folds/seeds; observable output must not depend on hash-table or worker order.

## Draft selection and acceptance gate

This section should be copied into a frozen spec before code or result generation. A sensible
small development grid is `alpha_group in {5,20,80}` and `alpha_global in {0,5,20}`. Select once
on D0 by mean held-out negative log likelihood, breaking ties toward less shrinkage. Carry the
chosen pair unchanged into D1 and D2′.

Recommended PASS conditions:

- on combined untouched D1+D2′ targets, hierarchical negative log likelihood improves at least
  5% over the best current non-hierarchical mode (`pooled` or `blend`);
- aggregate top-1 does not fall by more than 0.5 percentage point on either dataset;
- on D2′'s lowest-UMI cell quartile, top-1 improves by at least 1 point **or** negative log
  likelihood improves by at least 10%;
- real groups outperform size-preserving shuffled groups by at least half of the hierarchy's gain
  over the non-hierarchical baseline;
- calibration error and abundance-weighted error do not worsen materially;
- median wall time is no more than 1.5× packed EM and peak RSS remains below 6 GiB on D1.

MARGINAL means hierarchy improves proper scoring but not hard assignments or sparse-cell strata;
keep it as analysis, not the production default. STOP means shuffled labels perform similarly,
confirmatory log loss does not improve, or rare groups worsen; retain pooled emission and report
that sample-wide sharing was sufficient.

## Ranked goals after hierarchical EM

2. **Cell/group-scoped replay and queries, then junction-set/PSI.** One shared barcode/group file
   should filter region, junction, APA, discovery, and replay consistently. This has more product
   value than another single query command and supplies the grouping contract hierarchical EM
   needs.
3. **Compiled annotations and batched query plans.** Cache a digest-keyed compiled GTF and share
   chunk decoding across query sets. This is the next plausible speed win after inner-loop replay
   rearrangements failed.
4. **Unblock GEM-X 5′ v3 replication.** The current mouse 5′ v2 Gene result is good, but a broad
   5′ claim remains blocked by the unresolved raw-data source and marginal velocity components.
5. **Sparse E2/E3 only behind a biological query.** Prototype edits/residual sequence only after
   selecting an allele/edit/fusion question with orthogonal truth and a <=25% storage-premium
   gate. Do not begin with a container redesign.

Not next: another entropy codec, GPU EM, atlas-shared dictionaries, or a v2 archive rewrite. None
has a measured opportunity comparable to the hierarchical scientific gate.
