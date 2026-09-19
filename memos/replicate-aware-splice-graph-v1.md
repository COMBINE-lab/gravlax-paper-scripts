# Replicate-aware molecular splice graph v1

## Verdict

The population layer passes its frozen engineering and synthetic-inference gate.
Gravlax commit `173d9256b48eb4d729e73732b8a0b7e0cf01c3ff` adds
`aie cohort splice-graph`. The command defines a recurrent coordinate-edge
catalogue, reduces exact molecular path fragments independently in every
archive, and emits complete zero-explicit sample × path and sample × edge
matrices. With a strict two-condition design it fits path-versus-rest
beta-binomial contrasts whose only replicate is one unique biological-sample
archive row. Cells and molecules never become replicates.

The real eight-donor FNBP1 run is deliberately counts-only. It passes as an
engineering check and emits no p-value. The prospectively registered
`STOP_UNDERPOWERED` primary result and `STOP_DESCRIPTIVE` microglial result are
unchanged.

## Model and failure guards

The strict design header is
`sample<TAB>condition<TAB>archive<TAB>cells`; paths are relative to the design
file, `.` selects all archive cells, and both sample IDs and resolved archive
paths must be unique. Inference requires exactly the two conditions named by an
ordered `A:B` contrast and reports fitted usage in B minus A. Counts-only mode
allows a neutral single-condition panel and cannot be combined with a contrast.

A coordinate enters the common edge catalogue when local support is at least
`min_support` in at least `min_edge_samples` archives. Selection does not turn a
below-threshold local catalogue into an imputed zero: every archive is reduced
against the common edges, and genuine observed zeros are explicit in the
complete matrix. A sample enters a strand-specific contrast only when its total
common-graph path-fragment depth reaches `min_sample_umis`; ineligible rows stay
visible but are not fitted as zeros.

For path $j$ and eligible sample $i$, the response is $(k_ij,n_i)$: exact path
fragments versus all common-graph path fragments on the same strand. The null
has one mean and one beta-binomial concentration. The alternative has one mean
per condition and a shared concentration. Deterministic bounded profile maximum
likelihood yields a one-degree-of-freedom asymptotic likelihood-ratio p-value;
BH correction covers every tested path in the locus. The JSON retains raw
sample tuples, both likelihoods, fitted means/concentrations, effect, p, q, and
explicit reasons for skipped paths. This marginal path-versus-rest model is not
a joint Dirichlet-multinomial graph model.

The integration gate covers exact standalone/cohort agreement, repeated
representative deduplication through the shared primitive, complete zeros,
per-sample edge/path conservation, low-depth suppression, counts-only output,
duplicate archives, malformed designs, missing or mismatched contrasts, unknown
cells, all zero thresholds, union overflow, and incompatible genome signatures.
A balanced null is nonsignificant; a balanced four-versus-four usage shift has
the expected effect signs and smaller p-values. All 124 workspace tests pass,
and the documentation site builds successfully.

## Eight-donor FNBP1 counts-only result

The frozen design contains donors A--H from GSE234790 under one neutral
`cohort` condition. The eight tiny called-cell locus archives share genome
digest `2817ea0b…472`. At local support ≥1 and recurrence in at least two
archives, 44 union coordinates yield eight recurrent coordinates, ten
strand-specific edges, 15 nodes, and 12 exact path patterns. The complete
matrix contains 96 sample-path rows (63 zeros) and 80 sample-edge rows (49
zeros), demonstrating that absent observations remain explicit.

The paths carry 86 total UMI-class fragments: 75 on ten one-edge patterns and
11 on two two-edge patterns. One two-edge reverse-strand pattern has nine UMIs
across donors B, D, and G; the other has two UMIs across D and G. The strongest
single-edge pattern has 29 UMIs across seven donors. At the default ten-fragment
strand gate, only D, E, and G are eligible on the reverse strand and none is
eligible on the forward strand. This sparsity reinforces rather than rescinds
the registered donor-level STOP.

Twenty warm release-process invocations are byte-identical and take 0.006347 s
median (0.005930--0.007258 s), with 14,844 KiB maximum peak RSS. Every aggregate
path count and supporting-sample count reconciles with sample rows; every
sample edge count equals the sum of its path fragments containing that edge.

## Interpretation and next gate

This completes the requested query architecture, not the biological validation
of a population-splicing estimator. The synthetic beta-binomial gate establishes
correct plumbing and effect direction, but asymptotic calibration with two to
four low-depth replicates is not yet a paper claim. Promotion requires an
external, genuinely replicated two-condition dataset locked before locus
inspection, comparison with an established differential-splicing method, and
null calibration by sample-label permutation. A joint Dirichlet-multinomial or
hierarchical Dirichlet graph model should be considered only if that experiment
shows path covariance or unstable path-versus-rest multiplicity to be the
binding limitation.
