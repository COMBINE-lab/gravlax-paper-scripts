# Discovery-recall and streaming-hotpath optimization gate

**Locked:** 2026-08-31, before running either new discovery claiming mode on D0 or D1.
**Baseline code:** `0efabed0d3e92a72c432b4c815c2e0b70c7155e4`.

This is a post-v1 engineering gate. It does not alter the complete truth denominator, matching
rule, thresholds, or interpretation in `complete-denominator-discovery.md`.

## Discovery recall ablations

Run GENCODE-v32 discovery on D0 and D1 with three predeclared claiming rules:

1. `span` (compatibility/default): a molecule is claimed by any overlapping v32 transcript span;
2. `strand-span`: a molecule is claimed only by an overlapping transcript on the library-aware
   strand;
3. `compatible`: a molecule is claimed only when at least one stored placement is exon/junction
   concordant with a v32 transcript on the library-aware strand.

All other discovery settings stay fixed (`merge-gap=1000`, `min-umis=10`). Evaluate each mode
with the unchanged 367-gene primary denominator and deterministic one-to-one matcher. The baseline
must reproduce 164/367. Report recall and UMI-weighted recall overall and by clean,
antisense-only, and same-strand strata; candidate count; D1 candidate recurrence; runtime; and
peak RSS.

The preferred mode is the highest-recall mode satisfying both safeguards: no more than 3x the
baseline D0 candidate count and at least 95% D0-candidate same-strand D1 recurrence. An all-strata
gain of at least 15 percentage points is practically meaningful. If no mode satisfies the
safeguards, retain `span` as the default and report the ablations as failure analysis rather than
calling the larger candidate set improved discovery. Later-annotation confirmation remains a
narrow endpoint, never biological precision.

### Predeclared follow-up after the claiming ablations

If `compatible` residuals are joined into broad transitive components, test one structural remedy
before abandoning them: `residual-sites`. This mode is the union of the unchanged `span` calls and
evidence that span-claiming suppresses but transcript compatibility does not explain. Cluster the
second channel by the molecule's transcript-oriented 3' evidence coordinate; a site cluster is
bounded to `merge-gap` from its first coordinate rather than transitively chained. Require the
same 10 distinct UMI classes per channel and report candidate width (median, p95, maximum) so broad
intervals cannot game interval-overlap recall. Apply the same 3x-volume, 95%-D1-recurrence, and
+15-point meaningful-gain safeguards above. This follow-up was locked after observing the three
claiming-mode results but before implementing or running `residual-sites`.

For avoidance of ambiguity, residual-site candidate geometry is the half-open extent of the
observed transcript-terminal bases (minimum terminal base through maximum terminal base plus one),
not the full alignment or splice span. Thus a 1-kb bounded site cluster can never create a
megabase-scale truth overlap. This clarification was committed after the required width audit
rejected the first implementation's broad intervals and before rerunning the denominator with
corrected geometry.

The 10-UMI residual-site arm is diagnostic if it exceeds the volume guard. In that case, sweep a
small fixed residual-channel support grid of 25, 50, 75, and 100 distinct UMI classes while always
retaining the unchanged 10-UMI `span` channel. Select the highest-recall point that satisfies the
already locked 3x candidate-volume and 95% D1-recurrence safeguards; do not interpolate or add a
cutoff after seeing recall. Report the full grid, including failures. This support-grid extension
was locked after establishing that the 10-UMI residual channel contains the complete denominator,
but before calculating recall for any filtered point.

## Performance pass

Use the same binary, inputs, 24 pinned threads, and warm randomized/block-balanced runs where a
command lasts long enough for scheduling noise to matter.

- `replay-rows`: D0 and D1 GENCODE-v49 Gene replay. Output matrices must be byte-identical to the
  baseline commit. A useful change lowers median wall time by at least 10%, does not increase peak
  RSS by more than 5%, and preserves bounded streaming.
- `query discover`: D0 and D1 v32 `span` output must be byte-identical to baseline. A useful change
  lowers median wall time by at least 25% without increasing peak RSS by more than 10%.
- indexed queries: profile a fixed region query, exact-junction query, junction enumeration
  index-only query, and junction enumeration with exact cells. Outputs must be byte-identical.
  Optimize only a measured hot path; sub-100-ms commands are reported but are not paper claims.

Every reported timing includes command construction and output emission. Record separately the
timings already exposed by the binary (open/dictionary, annotation compile, decode/reduce, and
total) so an optimization is attributed to the correct stage. Do not trade correctness or archive
validation for speed.
