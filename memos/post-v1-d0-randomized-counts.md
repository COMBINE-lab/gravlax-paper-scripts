# D0 randomized warm and specific-file cold timing

Date: 2026-08-30  
Result: `results/post-v1-d0-randomized-counts.json`  
Protocol: `scripts/116_d0_randomized_counts_comparison.sh`  
Gravlax: `ef2126684378d8ba0dad5605ff81eed635bcab6a`

## Verdict

Randomized interleaving does not explain the D0 speed result. In five balanced two-arm crossover
blocks, Gene-only STARsolo with alignment output disabled took 82.44--84.60 s (median 83.71 s),
and Gravlax replay took 2.41--2.45 s (median 2.44 s): **34.31× by the ratio of medians**. Median
peak RSS was 34.93 GB versus 6.94 GB (**5.03×**). All five outputs from each arm reproduced their
separate frozen references exactly.

This is almost identical to the earlier sequential pilot's 34.58×, but with lower Gravlax
dispersion and the principal order confound removed. The descriptive first/second-position means
are 83.07/84.20 s for STARsolo and 2.430/2.443 s for Gravlax; five blocks are not enough to treat
those small differences as an inferential order-effect result.

## Cache control

Host-wide cache dropping is inappropriate on the shared server. The committed helper instead
opened only each arm's named inputs and issued `POSIX_FADV_DONTNEED` immediately before one cold
observation. The kernel accepted all calls: 3 Gravlax files totaling 3.55 GB and 19 STAR/FASTQ
files totaling 35.50 GB. The timed commands then reported nonzero file-system input counters,
unlike every warm replicate.

The cold observations were 2.75 s for Gravlax and 85.82 s for STARsolo (**31.21×**). Relative to
warm medians, cold overhead was 12.7% and 2.5%, respectively. The larger proportional Gravlax
effect is expected because I/O is a larger share of its short run. These are one-observation
diagnostics, and advisory eviction cannot prevent another user from recaching shared files, so no
distributional claim is made.

## Paper decision

Use 34.3× and 5.0× as the controlled D0 Gene-only warm result; retain 31.2× as a parenthetical
specific-file-cold observation. Remove the sequential-order caveat. Continue to state that this
compares different pipelines whose outputs have the separately measured 0.22% UMI-mass deviation;
"equal output capability" means Gene matrices, not identical matrices. D1 and D2′ remain the
binding generalization gate before a multi-dataset speed headline.
