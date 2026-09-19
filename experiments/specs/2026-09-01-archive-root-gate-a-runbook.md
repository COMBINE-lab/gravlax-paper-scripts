# Gate-A execution and validation runbook

This runbook implements Gate A of `2026-09-01-archive-root-and-jshape-routing.md`; the frozen
specification and machine-readable thresholds remain authoritative. No command below opens Gate B.

## Preconditions

- Gravlax and this reproducibility repository are clean, committed checkouts.
- The final binary implements `seal-archive`, `inspect-archive`, rooted source identity, and
  collection-build schema `gravlax.collection.build.v1`.
- Inputs match `manifests/archive-root-gate-a-inputs.tsv`.
- The destination run directory and compact result do not exist.
- `scripts/archive_root_adversarial.py` is passed as a fixture hook. A final run without the full
  accepted-boundary and rejected-corruption panel is non-promotable.

## Driver

From the project root, initialize the project-local toolchain, disable the iteration-only
chromosome setting, and build the release binary from the clean Gravlax commit:

```bash
source "$PWD/gravlax-paper-scripts/environment/setup.sh"
unset AIE_CHROM
CARGO_TARGET_DIR="$PWD/env/cargo-target-archive-root" \
  cargo build --manifest-path "$PWD/work/gravlax/Cargo.toml" --release -p aie
```

The confirmatory r3 invocation was:

```bash
python3 gravlax-paper-scripts/scripts/209_benchmark_archive_root.py \
  --project-root "$PWD" \
  --run-dir "$PWD/runs/post-v1/archive-root-gate-a-r3" \
  --aie "$PWD/env/cargo-target-archive-root/release/aie" \
  --code-repo "$PWD/work/gravlax" \
  --scripts-repo "$PWD/gravlax-paper-scripts" \
  --interface-status final \
  --fixture-hook "$PWD/gravlax-paper-scripts/scripts/archive_root_adversarial.py"
```

The driver refuses to overwrite the run, records exact argv/resource streams, and writes
`artifact-manifest.tsv` last. A reproduction must use a new empty run path; development runs cannot
be promoted.

## Fail-closed reduction

The confirmatory r3 reduction used the exact clean commits recorded by the driver:

```bash
python3 gravlax-paper-scripts/scripts/210_validate_archive_root.py \
  --run-dir "$PWD/runs/post-v1/archive-root-gate-a-r3" \
  --gravlax-commit fbd48e536b39729632f189a74e7ece891166ecd1 \
  --paper-scripts-commit df1e363d397c740300c97ec2b1f04fe48cd9eb6b \
  --artifact-inventory-out "$PWD/gravlax-paper-scripts/manifests/archive-root-gate-a-artifacts.tsv" \
  --out "$PWD/gravlax-paper-scripts/results/post-v1-archive-root-gate-a.json"
```

The reducer freezes the complete run-file set, independently recomputes every archive root,
checks the independently generated payload commitments, validates exact Gene/Velocity replay and
all four collection-query arms, enforces source-section and byte accounting, revalidates the
corruption panel, and applies every size/runtime/RSS threshold. It writes no result on failure.
On PASS it also writes the tracked excluded-artifact inventory for the complete run tree, binary,
nine rooted archives, all collection sidecars, and all replay matrix trees. Commit the compact
result and inventory, but do not commit any of the large artifacts they identify.

## Gate-B freeze boundary

Only after the compact Gate-A PASS and its artifact inventory are committed may the 96-row panel
be materialized:

```bash
python3 gravlax-paper-scripts/scripts/211_materialize_jshape_panel.py \
  --gate-a-result "$PWD/gravlax-paper-scripts/results/post-v1-archive-root-gate-a.json" \
  --out "$PWD/gravlax-paper-scripts/experiments/designs/archive-root-jshape-routing-panel.tsv" \
  --digest-out "$PWD/gravlax-paper-scripts/experiments/designs/archive-root-jshape-routing-panel.sha256" \
  --metadata-out "$PWD/gravlax-paper-scripts/experiments/designs/archive-root-jshape-routing-panel.metadata.json"
```

Commit the panel, digest, metadata, and the then-frozen serialization addendum before Gate-B code
or measurement begins. The materializer refuses to run without a complete Gate-A PASS and refuses
to overwrite outputs.
