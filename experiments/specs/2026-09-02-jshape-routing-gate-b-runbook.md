# Gate-B local-shape-routing execution and validation runbook

**Status: HARNESS FROZEN; NO CONFIRMATORY GATE-B MEASUREMENT OR RESULT EXISTS.** The panel,
serialization, controls, schedules, metrics, thresholds, and stop rules remain those frozen in
`2026-09-01-archive-root-and-jshape-routing.md` and its serialization addendum.

## Preconditions and fixed interfaces

- Gravlax and this reproducibility repository are clean, committed checkouts. The exact commits
  passed to the reducer must equal those recorded by the driver.
- The release `aie` binary is built once from that clean Gravlax commit and is used unchanged for
  both fallback and route arms.
- The eight A--H root-committed archives are the files bound by the Gate-A PASS, and the Gate-A
  compact result has SHA-256
  `3f75fb1f7d5cc8b726405017a320c3698b7384c5da1f833f6d1f0a29453429db`.
- The 96-row prospective panel has SHA-256
  `55980706cd6fbf0c07a1e652b330e5aca6946b1365010c8ef4e8cb6d6c4f1ca1` and remains unchanged.
- Python provides the `zstandard` module. The independent pure-Python BLAKE3 implementation must
  pass its known empty-input and `abc` vectors before any command runs.
- `/usr/bin/time` is available. The frozen execution uses 24 threads and eight alternating warm
  paired blocks; neither value is configurable in a promotable run.
- Script `214_jshape_routing_adversarial.py` is supplied as the fixture hook. Its exact named Rust
  tests provide positive geometry, representative, multimapper, and cross-chunk UMI predicates;
  independent collection mutations exercise malformed format, binding, context, reconstruction,
  checksum, truncation, and trailing-byte rejection through the release binary.

The source-controlled adapter fixes the following contracts:

| Operation | Command/interface | JSON schema |
|---|---|---|
| fallback build | `aie collection build ... --json` | `gravlax.collection.build.v2` |
| route build | `aie collection build ... --shape-routes --json` | `gravlax.collection.build.v2` |
| point query | automatic route use when present; exact fallback otherwise | `gravlax.collection.junction.v2` |
| junction set | automatic route use when present; exact fallback otherwise | `gravlax.collection.jset.v2` |
| region | route-neutral control | `gravlax.collection.region.v2` |
| inspect | `aie collection inspect COLLECTION` | `gravlax.collection.v4` |
| full route reconstruction | `aie collection inspect COLLECTION --verify-routes` | `gravlax.collection.v4` |

Build records separately expose the exact compressed source `shapes` payload bytes, attributable
source structure plus `shapes` bytes, and total construction source I/O. Query planning separately
exposes source identity/execution bytes, collection-sidecar bytes, compressed route-block payload
bytes, their total logical I/O, route/fallback archive counts, blocks loaded, and four phase times.

## Frozen driver

After the final code and harness commits are pushed, build one clean release binary and invoke the
driver from the project root with a new run directory:

```bash
python3 gravlax-paper-scripts/scripts/212_benchmark_jshape_routing.py \
  --project-root "$PWD" \
  --run-dir "$PWD/runs/post-v1/jshape-routing-gate-b-r1" \
  --aie "$PWD/work/gravlax/target/release/aie" \
  --interface-status final \
  --fixture-hook "$PWD/gravlax-paper-scripts/scripts/214_jshape_routing_adversarial.py"
```

The defaults resolve the two repositories, committed Gate-A result, frozen panel and metadata,
and A--H archives recorded by Gate A. Explicit `--archive SAMPLE=PATH` arguments may only replace
paths with byte-identical files. The driver refuses an existing run directory, dirty checkout,
changed input, non-final adapter, missing fixture hook, wrong thread/block count, or different
binary between arms. Its default `--interface-status draft` is an intentional second lock; a
confirmatory run must explicitly select `final`.

The workload builds route-free and routed collections in fresh, repeated, reversed-input, A--D
base, and E--H extension forms. It checks deterministic identities; inspects route-free, lazy
routed, and fully reconstructed routed forms; runs the known dense, sparse, absent, junction-set,
and region controls plus every fixed panel junction on both fresh-root and extension chains; and
executes eight alternating warm paired build and query blocks. It records command, stdout, stderr,
GNU-time, phase, source, sidecar, and total-I/O evidence. The artifact manifest is written last.
The driver can emit only `COMPLETE` or `DEVELOPMENT_ONLY`; it never emits `PASS`.

## Fail-closed reduction

Only after a complete clean run, substitute the exact clean commits recorded in `protocol.json`:

```bash
python3 gravlax-paper-scripts/scripts/213_validate_jshape_routing.py \
  --run-dir "$PWD/runs/post-v1/jshape-routing-gate-b-r1" \
  --gravlax-commit GRAVLAX_COMMIT \
  --paper-scripts-commit PAPER_SCRIPTS_COMMIT \
  --out "$PWD/gravlax-paper-scripts/results/post-v1-jshape-routing-gate-b.json" \
  --artifact-inventory-out "$PWD/gravlax-paper-scripts/manifests/jshape-routing-gate-b-artifacts.tsv"
```

The reducer never executes Gravlax. It validates the exact run-file manifest; reconstructs source
archive roots and committed `shapes` identities; validates every v2 build, v4 inspection binding
and canonical route descriptor; checks all scientific JSON after removing only frozen planning
fields; recomputes source/sidecar/total byte accounting, phase timing, process resources, paired
summaries, and every frozen size, construction, I/O, runtime, and RSS threshold; and requires all
26 named fixture predicates. The 26-case set intentionally excludes the formerly drafted
`maximum_safe_addition` label: checked-addition overflow is a required rejection, while ordinary
valid additions are already exercised by exact reconstruction.

On any discrepancy the reducer exits without either tracked output. Only after all checks pass
does it atomically expose the excluded-artifact inventory and compact `PASS`, writing `PASS` last.
The run tree, collections, binary, and rooted archives remain ignored large artifacts; only the
compact result and inventory are eligible to be committed.

## Current boundary

This runbook and scripts `212`--`214` freeze the measurement protocol after the route interface
was finalized. They contain no Gate-B treatment observation. No size, I/O, runtime, memory, or
integrity claim may be promoted until the single clean run is reduced to a committed PASS.
