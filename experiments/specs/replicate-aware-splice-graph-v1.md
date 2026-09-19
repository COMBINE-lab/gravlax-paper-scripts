# Replicate-aware molecular splice graph v1 gate

## Objective

Extend the exact per-archive molecular splice graph with a population layer
whose inferential unit is an explicitly named biological sample or donor. The
command must retain one exact numerator/denominator tuple per sample before
model fitting. Cells, UMIs, archive chunks, and technical lanes are never
replicates.

## Design input

The input is a strict TSV with the exact header
`sample<TAB>condition<TAB>archive<TAB>cells`. `sample` and `condition` use the
existing restricted ASCII identifier grammar. `archive` is one `.aie` file;
`cells` is either `.` for all archive cells or a strict one-barcode-per-line
scope file. Paths are resolved relative to the design file. Sample IDs and
resolved archive paths must both be unique, preventing one archive from being
split into artificial replicates.

Inference requires an explicit ordered `A:B` contrast and exactly those two
conditions in the design. `B-A` is the reported effect direction. Counts-only
mode permits one or more conditions and emits no test.

## Common graph and exact rows

A coordinate enters the common edge set when its archive catalogue support is
at least `min_support` in at least `min_edge_samples` archives. Once selected,
all archived evidence for that coordinate is reduced in every sample; a local
catalogue below the discovery threshold contributes observed evidence rather
than an imputed zero. Stamped archives must share a reference digest and may
not be mixed with unstamped archives.

For each sample and strand, a path fragment is the selected junction set
co-supported by one archive UMI class. The output contains the union of path
and strand-aware edge definitions plus a complete zero-explicit sample matrix
for both. Repeated representatives and chunks cannot count a class twice.
Every edge UMI count must equal the sum of sample path counts containing that
edge. `max_paths` is a hard union limit, never a truncation option.

## Eligibility and beta-binomial contrast

Each tested path is contrasted against every other observed path fragment on
the same strand. For sample `i`, `k_i` is the path count and `n_i` is the total
strand-path UMI count. A sample is eligible on that strand only when
`n_i >= min_sample_umis`; ineligible samples remain in the output but are not
fit as zeros. Each condition must retain at least `min_replicates` eligible
samples. A tested path must meet both `min_path_umis` across eligible samples
and `min_path_samples` eligible samples with nonzero support.

The beta-binomial null has one mean usage and one shared concentration. The
alternative has a separate mean for each condition and the same shared
concentration. Parameters are bounded away from degenerate zero/infinity and
fit by deterministic profile maximum likelihood. The likelihood-ratio test
has one asymptotic degree of freedom. BH correction is applied across all
tested paths in the locus. The output must preserve the raw per-sample
`(k_i,n_i)` rows, null/alternative likelihoods, fitted means and
concentrations, LRT, p-value, q-value, and explicit `B-A` effect.

This is a path-versus-rest marginal model, not an omnibus graph model. A full
Dirichlet-multinomial or hierarchical Dirichlet model remains the candidate
extension when joint path covariance is scientifically necessary.

## Acceptance

- Synthetic repeated representatives, opposite strands, and multi-edge paths
  match standalone `splice-graph` counts exactly.
- The complete sample path and edge matrices are deterministic and explicitly
  contain zeros.
- Duplicate samples, duplicate resolved archives, malformed design rows,
  unknown cells, reference mismatch, a missing/ambiguous contrast, zero
  thresholds, and union overflow fail closed.
- A balanced null simulation is not significant and a balanced strong usage
  shift has the expected effect sign and a smaller p-value.
- Removing eligibility from one condition suppresses the test rather than
  manufacturing zero-count replicates.
- Counts-only mode emits exact rows and no p-values.
- All legacy cohort and query output remains unchanged.

## Real-data interpretation gate

The eight GSE234790 FNBP1 donor archives are an engineering counts-only check.
All design rows use one neutral condition and no contrast. Their previously
registered `STOP_UNDERPOWERED` and `STOP_DESCRIPTIVE` outcomes are immutable;
graph rows cannot retroactively relax those thresholds or create a case/control
comparison.
