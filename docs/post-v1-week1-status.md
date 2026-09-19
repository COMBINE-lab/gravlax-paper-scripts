# Post-v1 week-one status

Date: 2026-08-30

This work begins after, and does not modify, the coordinated
`gravlax-v1-8-30-2026` snapshot.

## Local review commits

- Gravlax code and public docs: `b5a791e681e7f1f93d46459745eedba59ac841d6` on `post-v1-week1`
  (`work/gravlax`, based on tagged code commit
  `d9257ca2c81569aa14dfeb2899beb4b90ebb5a9d`).
- Manuscript: `1a1ef0bf3501ee57c009e37661c773eff88ff046` on `post-v1-week1`
  (`work/moleceular-evidence-store-paper`, based on tagged paper commit
  `67da17775aed087a16ab9dbd1aea91e3eca7160a`).
- Reproducibility specifications: this branch and commit.

These commits are local review commits and have not been pushed to the authoritative GitHub
repositories.

## Claim audit

| claim family | permitted statement | evidence or gate | current action |
|---|---|---|---|
| archive serialization | archive- and BAM-sourced Gravlax matrices are byte-identical after the same row reduction | new `crates/aie/tests/cli_golden.rs`; full-data regression logs | retained and scoped |
| fresh STARsolo | Gene replay differs by 0.22–0.45% of UMI mass on four evaluated datasets | `results/archive-fidelity-full-v49.json`, `results/gate1a-d1-v49.json`, `results/d2-fid.json`, `results/d2p-fid.json` | “exact” removed |
| evidence sufficiency | conditional on fixed placements, barcode correction, compatible GTF coordinates, and a named consumer that reads only retained fields | `docs/formal-conditional-sufficiency.md`; necessary-field counterexample plan | premise now explicit |
| deployed two-representative reduction | empirical coarsening, separately priced against full-read Gravlax replay | `results/chainrep-full-v49.json`, `results/tworep-full-v49.json` | no theorem claim |
| minimality | not established | would require a lower-bound proof over the consumer family | headline removed |
| storage | 11–18 bits/read and observed absolute artifact sizes; BAM/CRAM are information-richer | existing size logs and manifests | “every baseline retains less” removed; matched experiment specified |
| runtime | D0 warm Gene-only counts-only STARsolo median 83.68 s versus Gravlax 2.42 s (34.6×); five runs per arm | `results/post-v1-d0-counts-only-pilot.json` | stale 150–370× headline removed; randomized cold/warm and D1/D2' confirmation required |
| EM compatibility | implemented STARsolo-compatible update scheme; current comparison is nonzero | EM logs and results | unaudited byte-exact wording removed |

## Code hardening completed

- Reject seekable container versions newer than the reader.
- Validate directory bounds, duplicate names, inline headers, lengths, and integer conversions.
- Reject overlong varints and malformed/trailing rANS residual streams.
- Bound rANS decode counts with chunk-derived metadata before allocation.
- Replace chunk stream indexing panics and coordinate/class arithmetic wraparound with errors.
- Reject malformed exon records, zero/reversed GTF intervals, and unsupported strand values.
- Add a generated tagged-BAM → `.aie` → matrix CLI golden test covering repeated unique evidence
  and a primary/secondary multimapper; compare the output byte-for-byte with BAM-sourced replay.
- Correct public CLI examples and describe the currently reserved `replay` and `eval` crates
  honestly.

## Validation

- `cargo test --release --workspace`: **86 passed, 0 failed**. The 18 pre-existing compiler
  warnings remain a tracked quality item.
- Astro/Starlight static documentation build: succeeded under Node 22.18.0; 15 pages built.
- Typst manuscript compilation: succeeded; output written only to `/tmp` (13 pages).
- `git diff --check`: clean in all three working repositories.
- `experiments/manifests/post-v1-gates.yaml`: parsed successfully with `js-yaml`.

## D0 baseline pilot

The ingest-equivalent D0 field transform passed on 2026-08-30. Full- and reduced-BAM ingests were
byte-identical; genome stamping reproduced the frozen D0 archive; and archive/full-BAM/reduced-BAM
replays matched the frozen matrix outputs byte-for-byte. The reduced archive-mode CRAM is
758,514,044 bytes, 6.83× the 111,087,381-byte `.aie`, while still retaining more information.
See `memos/post-v1-d0-fair-baseline-pilot.md` and the compact JSON result.

The counts-only D0 arm now has five warm observations per method. The paper reports the 34.6×
ratio with the explicit caveat that arms were sequential, not randomized. Next obtain a safe
cold-cache observation, randomized interleaving, and D1/D2' confirmation, then run the accepted
storage transform on D1 and D2'.
