# Coordinated Gravlax v1 snapshot

**Snapshot tag:** `gravlax-v1-8-30-2026`

**Snapshot date:** 2026-08-30

**Archive format:** `.aie` v1.5 (`meta.codec = "rans2"`, ten chunk streams, 65,536-class
cell-of-class blocks)

## Authoritative revisions

| Component | Repository | Commit |
|---|---|---|
| Code | `COMBINE-lab/gravlax` | `d9257ca2c81569aa14dfeb2899beb4b90ebb5a9d` |
| Paper | `COMBINE-lab/moleceular-evidence-store-paper` | `67da17775aed087a16ab9dbd1aea91e3eca7160a` |
| Experiments | `COMBINE-lab/gravlax-paper-scripts` | the commit addressed by this snapshot tag |

The analysis workspace's `src/` tree was compared recursively with the tagged Gravlax checkout;
the only checkout-only paths were repository metadata, documentation, CI, assets, license, and
README. Source, workspace manifests, and lockfile were identical.

The tagged paper source compiled successfully with Typst 0.15.1 to an 11-page, 272,867-byte
temporary PDF. The `main.pdf` already tracked at the specified paper commit is an earlier 10-page
render; `main.typ` and the tagged figure sources are authoritative for this snapshot.

## Snapshot validation contract

- The Gravlax workspace test suite must pass at the code revision above.
- The paper must compile from the paper revision above without modifying its tracked PDF.
- All archived shell and Python drivers must parse successfully.
- All archived JSON outputs must be valid JSON.
- Code/paper/archive replay identity is defined precisely in `PROVENANCE.md`: archive-sourced and
  BAM-sourced Gravlax matrices are byte-identical; comparisons to fresh annotation-aware STARsolo
  are high-fidelity but not byte-identical.

Validation performed on 2026-08-30:

- `cargo test --release --workspace`: **77 passed, 0 failed** using
  `CARGO_TARGET_DIR=/tmp/gravlax-snapshot-target` (the build emitted 18 non-fatal compiler
  warnings, retained as a software-quality item for the red-team review).
- Typst compilation: succeeded to `/tmp/gravlax-v1-8-30-2026-paper.pdf`.
- Shell syntax, Python AST parsing, and JSON parsing: all passed.
- Credential-pattern scan: no matches; no committed snapshot file exceeds 5 MB.

## Scope boundary

This tag freezes the weekend prototype and the complete paper/results state immediately before the
independent red-team and v2 work. Later commits may deliberately break archive compatibility; this
tag remains the recovery point for the v1 claims and artifacts.
