# Post-v1 complete-denominator discovery

**Date:** 2026-08-31  
**Protocol lock:** `522bad43daf8d22c5141d7c5b4eb3f8fe9ca6563`  
**Result:** clean-locus **PASS**; all-strata **MARGINAL**; accounting and determinism **PASS**

## Executive result

The earlier 52/52 prospective result was real but conditional: it excluded every future gene
whose span overlapped a v32 gene. Under the locked complete denominator of all 367 v49-added
genes with at least 50 D0 Gene UMIs, deterministic one-to-one recovery is **164/367 (44.7%)**;
UMI-weighted recall is **42.6%**. This is inside the pre-registered marginal band, not a pass.

The failure is sharply mechanistic:

| v32-overlap stratum | recovered | recall | UMI-weighted recall |
|---|---:|---:|---:|
| no span overlap (`clean`) | 51 / 52 | **98.1%** | 97.0% |
| opposite-strand overlap only | 32 / 60 | 53.3% | 60.7% |
| at least one same-strand overlap | 81 / 255 | **31.8%** | 31.5% |
| **all** | **164 / 367** | **44.7%** | **42.6%** |

The clean result differs from the old 52/52 only because two truth genes compete for one emitted
locus under the new one-to-one event rule; permissive any-overlap recovery remains 52/52. At 10
and 100 D0 UMIs, total recall is 47.5% and 45.8%, respectively, so the conclusion is not an
artifact of the primary 50-UMI cutoff.

## Why the denominator members are credible

All 367 primary truth genes are independently expressed in D1 at at least 10 UMIs. PolyASite 2.0
supports 314/367 (85.6%), and a compatible D0 STAR splice junction supports 148/367 (40.3%);
the latter is expected to miss single-exon genes. The 203 unmatched genes are not merely a weak
tail: all are expressed in D1, 166 (81.8%) have PolyASite support, and 87 (42.9%) have STAR splice
support. Only 19/203 are emitted by the independent D1 v32 discovery, confirming that the old
annotation's claiming rule—not D0 sampling alone—is the main limiter.

## Candidate-side result

D0 emits 12,389 v32-discovery loci at the fixed 10-UMI floor. One-to-one matches to the primary
truth set are deliberately only 164/12,389 (1.32%; ranked-event AP 0.049). This is an
**annotation-history confirmation endpoint**, not biological precision: another 2,145 loci
overlap lower-expression or otherwise non-primary v49-added genes, and absence from a later
annotation does not make an expressed locus false.

The orthogonal observation is recurrence. **12,204/12,389 (98.5%) of D0 loci and 99.4% of D0 UMI
mass overlap a same-strand D1 discovery locus**. Even the 8,455 later-annotation-unresolved loci
recur at 98.3% (99.3% of mass). This is unusually high because the second PBMC library is deeper
and the 1 kb merge rule yields broad, shared transcription/readthrough intervals. It supports
reproducibility, not novel-gene identity, and should be stated that way.

## Claim correction and next design implication

The defensible claim is now:

> Gravlax conservatively rediscovers essentially every sufficiently expressed future gene in
> previously unoccupied loci, but only 45% across a complete denominator; same-strand overlap with
> an existing annotation is the dominant blind spot.

Do not claim complete discovery, de-novo gene calling, or biological precision from this test.
The result gives a concrete capability target: add a strand-aware `discover --allow-overlap`
mode that can retain unclaimed antisense and local splice/3′-end evidence inside old gene spans,
then rerun this *unchanged* denominator. It also strengthens the rationale for junction
enumeration: overlap-aware candidate models need direct access to splice structure rather than
only merged locus spans.

## Reproduction

```sh
AIE_BIN=/path/to/gravlax/target/release/aie \
  bash scripts/136_complete_denominator_discovery.sh
```

The runner regenerates D0 and D1 v32 calls, evaluates twice, and byte-compares the JSON and both
audit TSVs. All three non-timing outputs reproduced byte-for-byte.

