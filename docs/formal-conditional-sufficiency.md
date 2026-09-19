# Conditional sufficiency and exactness taxonomy

Status: post-v1 specification, 2026-08-30. This document sharpens the claim language; it does not
retroactively alter the `gravlax-v1-8-30-2026` snapshot.

## Objects and scope

Let:

- `R` be the sequenced reads and raw barcode/UMI observations;
- `A(R; H, P)` be genome placements produced against reference `H` with mapper parameters `P`;
- `B(R; W, K)` be the barcode correction map produced with whitelist `W` and correction policy `K`;
- `g` be a GTF whose contig identities and coordinates are compatible with `H`;
- `C(A, B, g; S)` be a named assignment-and-collapse consumer with semantics `S`;
- `T_E1(A, B)` retain blocks, junctions, strand, multimapping structure, corrected cell identity,
  and the UMI equality/adjacency relations consumed by `S`.

The target family is the implemented STARsolo Gene/GeneFull/Velocyto-style consumer under frozen
alignment, barcode correction, strandedness, UMI, and ambiguity policies. “Any annotation” means
any `g` in this compatibility domain. It does not mean another genome build, a different mapping
policy, a different whitelist/corrector, or a consumer that reads sequence, base qualities, edit
positions, or alignment scores.

## Conditional result to prove and test

If `C` reads only fields retained by `T_E1`, then for every compatible `g`:

```text
C(A, B, g; S) = C(T_E1(A, B), g; S).
```

The proof obligation is a dependency trace from every branch in assignment, molecule ambiguity
resolution, and UMI collapse to a retained field. The implementation test is a read-level golden
path over adversarial annotations, not a comparison to a newly generated alignment.

The statement is conditional sufficiency, not global minimality. A minimality result would require
an equivalence relation over all inputs and a lower-bound proof that no strictly coarser statistic
preserves the full consumer family. Gravlax v1 does not provide that proof.

## Necessary-evidence counterexamples

For each field, construct two inputs that become indistinguishable when that field is removed but
produce different consumer output:

| removed evidence | counterexample family | changed decision |
|---|---|---|
| strand | overlapping antisense genes with identical spans | gene identity |
| exon blocks | gapped versus contiguous placement over an intron | exon/intron compatibility |
| junction endpoints | annotated versus one-base-shifted splice boundary | concordance |
| alternative placements | alternatives falling in one versus two genes | unique-gene union |
| cell identity | same UMI observed in two cells | molecule namespace |
| UMI equality | duplicate observations versus distinct UMIs | molecule count |
| 1-mismatch adjacency | abundant/rare neighboring UMIs versus non-neighbors | correction/collapse |
| read multiplicity | equal geometry with unequal support | best-gene and EM weight |

These establish necessity for the named consumer through counterexamples. They do not prove that
the current encoding is byte-minimal.

## The deployed coarsening

The v1 archive stores `Q(T_E1)`: within a molecule and junction chain, all reads are replaced by
up to two span-extreme representatives plus a weight. `Q` is not claimed to preserve `C` for every
compatible GTF. Report three comparisons separately:

1. **Serialization identity:** archive-sourced and BAM-sourced Gravlax replay apply the same `Q`
   and must emit byte-identical matrices.
2. **Reduction error:** `C(Q(T_E1), g)` versus Gravlax replay retaining every read placement.
3. **Fresh-pipeline deviation:** `C(Q(T_E1(A_no_anno)), g)` versus fresh annotation-aware STARsolo,
   which also changes `A` and may change other upstream decisions.

Only comparison 1 is “byte-identical” or “exact.” Comparisons 2 and 3 are deviations with a named
denominator and metric. Component-wise velocity deviations must not be collapsed into the Gene
number.

## Release gates

- A claim table maps every manuscript use of “exact,” “sufficient,” “minimal,” and “any
  annotation” to one of the definitions above.
- Adversarial golden cases cover every counterexample row.
- The BAM-to-archive-to-matrix CLI test proves serialization identity on unique, repeated, and
  multimapping evidence.
- Fresh STARsolo comparisons are always labeled end-to-end and report nonzero deviation.
- Any structural v2 coarsening must report its incremental effect against full-read Gravlax replay,
  separately from alignment-induced deviation.

