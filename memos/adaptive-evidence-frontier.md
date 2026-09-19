# Adaptive evidence-frontier pilot: result

**Verdict: STOP.**  The all-chain containment frontier improves D0 Gene and Velocity fidelity, but
misses both frozen resource gates.  The one pre-declared junction-only rescue clears the resource
gates but removes essentially all of the fidelity gain.  Neither mode should be promoted or merged
into the production archive.

## Provenance and frozen gates

- Gravlax base: `e378516a0fb6602e63d626d0c49fb11e58f463c8`
  (`post-v1-week2`).
- Worktree/branch: `work/gravlax-frontier-pilot`, `post-v1-frontier-pilot`.
- Initial specification SHA-256, frozen before implementation:
  `3f0a037aec4d393010cc79be1f4f34294b09ca3e6f85bb48bdedfc75f7c88505`.
- Specification plus focused-rescue addendum SHA-256, frozen before rescue measurement:
  `ba3113855652bbc88e83e3dbe43835f0237d0ac5fc88f59c02bb83ea6113eff3`.
- Promotion required archive growth at most 3.0%, replay wall time at most 1.15x, and either a
  0.015-percentage-point D0 Gene improvement or a 0.50-point single-component / 0.25-point mean
  Velocity improvement, with the frozen non-worsening limits.
- Measurements used `rustc 1.98.0 (88d9e12ae 2026-08-18)` on the same 256-logical-CPU Linux host
  within each paired comparison.

The source remains deliberately uncommitted and has not been pushed.

## Design and preserved predicate class

For the unique reads of one corrected `(cell, UMI, locus, strand, absolute junction chain)`, each
distinct geometry is the outer interval `[s,e)`.  Under interval containment, the treatment stores
the union of the observed containment-maximal and containment-minimal elements and stores the raw
chain read weight exactly.

This union is the unique smallest observed subset preserving all existential predicates that are
upward-closed or downward-closed under containment.  If an upward-closed predicate contains an
observed element, it contains a maximal element above it; the dual argument uses a minimal element
for a downward-closed predicate.  Conversely, each retained maximal or minimal element is isolated
among the observations by its corresponding principal set, so none can be removed while preserving
the full predicate class.  For a fixed absolute junction chain, this covers questions such as
whether any observed geometry contains a requested interval or is contained by one.  A fixed-
annotation concordance predicate is covered only when it is genuinely containment-monotone and the
consumer depends on existence, not representative multiplicity.

It does not preserve interior endpoint-band predicates, geometry multiplicities, arbitrary
Velocity state changes at exon boundaries, sequence/edit/quality predicates, or complete isoform
phasing.  Executable tests include an interior-band counterexample and a multiplicity
counterexample.

The focused rescue applies the identical frontier only to chains with at least one alignment splice
junction; empty-junction chains retain the archive-v1 two extremes exactly.  This predicate is
annotation-independent, but its formal guarantee is correspondingly restricted to
junction-containing chains.

## Implementation

- `aie ingest-archive --unique-representatives containment-frontier` selects the all-chain pilot.
- `--unique-representatives junction-containment-frontier` selects the focused rescue.
- With no option, `two-extremes` remains the default and emits the exact archive-v1 bytes.
- Chains with one or two representatives retain the existing `(weight << 1) | has_second` token.
  Wider pilot chains use a bit-48 escape plus exact `u32` weight and `u16` representative count;
  reserved bits and counts below three fail closed.  This is experimental reader support, not a
  production format-version proposal.
- Ingest prints exact weight and the complete representative-count histogram.  Debug output splits
  the distribution into empty-junction and junction-containing chains.

## D0 (`pbmc_1k_v3`) result

Called-cell fidelity is against fresh STARsolo v49.  Gene replay time is the median of three paired
sequential runs; percentages for Velocity are mass-weighted relative L1.

| arm | archive bytes | growth | representatives | chains >2 / max | Gene replay | Gene movement | Velocity S / U / A |
|---|---:|---:|---:|---:|---:|---:|---:|
| two extremes | 111,079,672 | baseline | 33,118,838 | 0 / 2 | 3.13 s (all-chain pair); 3.17 s (rescue pair) | 0.243904% | 1.807734 / 0.519740 / 3.058809% |
| all-chain frontier | 132,487,429 | **+19.272% FAIL** | 48,640,255 | 7,051,867 / 80 | 3.92 s (**1.252x FAIL**) | 0.215574% | 1.588026 / 0.516139 / 2.042336% |
| junction-only rescue | 114,201,153 | +2.810% PASS | 34,130,033 | 664,661 / 15 | 3.25 s (1.025x PASS) | **0.244068% FAIL** | **1.807649 / 0.519740 / 3.058679% FAIL** |

The all-chain arm improves Gene movement by 0.028330 percentage points and mean Velocity error by
0.413260 points, satisfying the accuracy alternative, but its resource failures are decisive.  The
rescue worsens Gene movement by 0.000164 points and improves its best Velocity component by only
0.000130 points (mean improvement 0.000072 points), orders of magnitude below the gate.
Median Gene replay peak RSS is 1,370,936 KiB for the original paired control, 1,740,264 KiB for the
all-chain arm (1.269x), and 1,407,104 KiB for the rescue (1.026x against its paired control).

The all-chain distribution localizes the cost:

- empty-junction: 16,756,542 chains, 48,142,141 exact reads, 41,075,791 representatives,
  6,387,206 chains above two, maximum 80;
- junction-containing: 4,844,966 chains, 8,838,888 exact reads, 7,564,464 representatives,
  664,661 chains above two, maximum 15.

Thus, the endpoint alternatives responsible for the measurable fidelity benefit are overwhelmingly
in the empty-junction evidence that the focused rescue must omit to remain small.

## Brain nuclei D2' result

| arm | archive bytes | growth | representatives | chains >2 / max | Velocity replay | Velocity S / U / A |
|---|---:|---:|---:|---:|---:|---:|
| two extremes | 595,246,459 | baseline | 173,703,807 | 0 / 2 | 177.52 s | 2.788983 / 0.574331 / 3.946257% |
| all-chain frontier | 640,092,167 | **+7.534% FAIL** | 206,669,195 | 20,614,311 / 122 | 193.70 s (1.091x PASS) | **2.647599 / 0.572579 / 3.711359% FAIL** |

The largest D2' component improvement is 0.234898 percentage points and the mean improvement is
0.126011 points, so the all-chain arm misses both Velocity accuracy alternatives in addition to its
storage failure.  The pre-declared rescue cannot exceed the all-chain evidence set and was not run
on D2' after it failed the D0 accuracy gate.  Velocity replay peak RSS rises from 41,479,468 KiB to
44,078,308 KiB (1.063x).

## Correctness and determinism

- Exact D0 chain weight is 56,981,029 in every arm; exact D2' chain weight is 231,711,883.
- The current pilot binary's default D0 archive is byte-identical to the control:
  `a664edd16cecc86347341bff177f184f67ad53007d3c194ed8b41f0f4d49c82c`.
- The D0 all-chain treatment is
  `611e073550a0184f3c6c6cd6ff37d1c864375caa1eade4d5ba3b1c1c29f61f25`.
- Two independent D0 rescue ingests are byte-identical:
  `337a7fe802131a98abb176ad67e8624740dd275b74e80654e896525efe1a1cd4`.
- Two independent development all-chain ingests are byte-identical:
  `108947cc7c87d4903544572129c56252521a07e21d637c03e1c1d2e4a5b8bcfe`.
- D0 rescue eager and streaming Gene matrices are byte-identical:
  `397259ebcb5889f50323cfd707f4b125c3e14b45efe96122fa83031da6d5cd17`.
- The D2' control/treatment archives are respectively
  `b493523e1f6a05df05c656ca15c8791b2d5da2e42c8fe700950030c9886dc06c` and
  `0fb6fef789f59a17ef45caa1253173f55a0531ac76dce34cf2e171d35c3b697c`.
- The full workspace test suite passes: 143 tests, including five containment/proof/counterexample
  tests and four wide-token/malformed-stream tests.  `git diff --check` passes.

## Artifact root

Raw local artifacts are under
`${GRAVLAX_PROJECT_ROOT}/runs/post-v1/adaptive-frontier-pilot-r1`:

- `d0-control`, `d0-control-current`, `d0-frontier`, `d0-junction-frontier`, and
  `d0-junction-frontier-repeat` contain archives, ingest logs, timings, and histograms;
- `d0-replay-*` and `d0-velocity-*` contain matrices, timings, and the compact Gene comparison
  JSON files;
- `d2p-frontier`, `d2p-control-debug`, and `d2p-velocity-*` contain the brain-nuclei archive,
  representative audit, matrices, and timings.

The large generated artifacts are ignored and must not be committed.

## Decision

Do not merge or promote either adaptive in-archive frontier.  The experiment identifies a real but
unfavorable trade: the unspliced endpoint diversity that improves replay is also the dominant byte
and replay-time cost.  If revisited, use an explicitly optional, query-targeted or
annotation-conditioned disposable sidecar rather than weakening these gates or expanding the core
archive.
