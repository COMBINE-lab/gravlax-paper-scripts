# Frozen pilot: adaptive unique-read geometry frontier

**Frozen:** 2026-09-01, before implementing the frontier or opening any frontier result.
**Base:** `e378516a0fb6602e63d626d0c49fb11e58f463c8` (`post-v1-week2`).
**Status:** exploratory archive-v2 pilot; the default two-representative encoding remains unchanged.

## Object being tested

Within one corrected `(cell, UMI, locus, strand, absolute junction chain)`, represent a unique-read
geometry by its half-open outer interval `[s,e)`.  Duplicate intervals are one geometry because,
conditional on the fixed absolute junction chain, they reconstruct the same aligned blocks used by
E1 assignment.

Interval `a` **contains** interval `b` when `a.s <= b.s && a.e >= b.e`.  The proposed frontier is
the union of:

- the containment-maximal geometries (not strictly contained by another observed geometry); and
- the containment-minimal geometries (not strictly containing another observed geometry).

This is the unique smallest observed subset preserving every existential predicate that is either
upward-closed or downward-closed under containment.  Every retained maximal (minimal) element can
be isolated by an upward-closed (downward-closed) principal set, so removing it breaks that class
of predicates.  The raw chain read count is retained exactly.  The frontier is annotation-free;
no GTF, gene, transcript, cell type, or group label participates.

This guarantee does **not** cover predicates confined to an interior endpoint band, counts of reads
satisfying a predicate (geometry multiplicities are not stored), velocity state transitions at an
arbitrary exon boundary, sequence/edit/quality-dependent consumers, or complete isoform phasing.
The implementation must include executable counterexamples for at least the first two limitations.

## Arms and immutable inputs

- **Control:** current two lexicographic span extremes, written by the same pilot binary with no
  new option.
- **Treatment:** `--unique-representatives containment-frontier`.
- **D0:** `pbmc_1k_v3`, annotation-free BAM, existing v49 Gene and Velocyto oracles and called-cell
  list.
- **D2':** `Brain_3p` nuclei, attempted if the D0 implementation and roundtrip gates pass.
- Same BAM, whitelist, GTF, barcode universe, chunk size, zstd level, threads, executable, and host
  within each paired comparison.

## Correctness gates (all mandatory)

1. The option is opt-in; the absent/default option emits a byte-identical archive to the current
   encoder on D0 (allowing only an explicitly documented pre-existing genome-meta difference).
2. Archives containing more than two representatives round-trip exactly through eager decoding;
   malformed representative counts fail closed.  All existing Rust tests and new frontier,
   roundtrip, and counterexample tests pass.
3. D0 archive replay remains byte-identical between eager and streaming paths.  Frontier
   representative order and bytes are deterministic across two ingests.
4. Report chain-weight conservation and the full representative-count distribution.  Any chain
   loses the exact read count => STOP.

## Measurable-value gate

The pilot is promoted only if the correctness gates pass, archive size is no more than **3.0%**
above the paired control, and replay wall time is no more than **1.15x** the paired control, plus at
least one of:

- **Gene:** D0 called-cell UMI-mass movement versus fresh STARsolo improves by at least **0.015
  percentage points** (half the previously measured 0.031-point gap between two representatives
  and the full-read ceiling), with no tested annotation worsening by more than 0.005 points; or
- **Velocity:** on D2' (preferred) or D0, at least one component's called-cell mass-weighted
  relative L1 improves by **0.50 percentage points**, or the mean across S/U/A improves by **0.25
  points**, with no component worsening by more than 0.10 points.

If storage exceeds 3%, or neither accuracy threshold is met, retain the implementation and results
as a negative design experiment but do not make it the archive default.  A failure may motivate a
future annotation-conditioned disposable sidecar, but must not be rescued by tuning this frozen
gate.

## Frozen focused rescue: junction-only frontier

**Frozen:** 2026-09-01, after the all-chain frontier missed the storage and replay gates and before
implementing or measuring this rescue arm.

The one permitted rescue is `--unique-representatives junction-containment-frontier`.  For a fixed
absolute junction chain containing at least one splice junction (equivalently, its observed shapes
have more than one aligned block), retain the containment frontier defined above.  For an empty
junction chain, retain the archive-v1 two span extremes exactly.  The predicate is derived only
from alignment block geometry and remains independent of annotations, genes, transcripts, groups,
and cell labels.

The formal monotone-containment guarantee therefore applies unchanged to junction-containing
chains and is deliberately not claimed for empty junction chains.  Exact chain read weight remains
mandatory in both cases.  This pre-declared restriction targets the 6,387,206 D0 empty-junction
chains with more than two representatives that dominated the failed all-chain frontier's growth;
no further data-dependent thresholds or tuning are allowed.

The original correctness, **+3.0% archive-size**, **1.15x replay-time**, and Gene/Velocity accuracy
gates remain unchanged.  Failure of either resource gate, loss of exact weight, or failure to meet
an accuracy threshold is a final **STOP** for adaptive in-archive geometry frontiers in this line of
work; the gates will not be diluted.
