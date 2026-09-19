# Optional federated collection index v2

## Verdict

**PASS.** A deterministic 8,931,229-byte derived sidecar routes exact junction, region, and
junction-set queries across eight SEZ archives totaling 1,326,036,691 bytes. The sidecar is
**0.674%** of the source archives (148.5-fold smaller). A true immutable four-archive base plus
four-archive extension occupies 9,469,771 bytes, **0.714%** of the sources. Building the extension
does not change one byte of the base, and building the fresh root with reversed input order produces
the same SHA-256, `9b33caf1...05245`.

The final benchmark is `runs/post-v1/collection-index-v2-r5`, built from clean Gravlax commit
`6b92503cffbd48da8312bdd75c1f1f70d5ff2a27` with binary SHA-256
`3c5faaf38d8b9ea14001a8561b78b026452804d454b16a5f2e11dbc089605a4d`. Runs r1 and r2 were failed
driver attempts; r3 and r4 were successful dry runs used to harden validation. None is a scientific
comparison arm. R5 alone is the locked result and adds matched aggregate timings for the eight
ordinary per-archive queries.

## Exactness and routing result

The fail-closed validator checks every per-sample result against an ordinary query of its source
archive: eight samples times dense junction, sparse junction, absent junction, region, and
junction-set equals **40 exact comparisons**. It also requires fresh-root and base-plus-extension
answers to be logically identical for all five queries. The fixed scientific totals are:

| query | exact result | root wall / RSS | eight-query naive wall / RSS | opened / pruned archives | actual archive bytes read |
|---|---:|---:|---:|---:|---:|
| dense FNBP1 junction | 182 UMIs, 180 cells | 0.13 s / 402,700 KiB | 0.45 s / 156,320 KiB | 8 / 0 | 36,588,673 |
| sparse neighbor | 4 UMIs, 4 cells (A and G only) | 0.12 s / 178,528 KiB | 0.24 s / 156,312 KiB | 2 / 6 | 14,100,036 |
| absent neighbor | 0 UMIs, 0 cells | <0.01 s / 12,448 KiB | 0.11 s / 18,676 KiB | 0 / 8 | 0 |
| NTRK2 region | 74,848 molecules, 74,352 UMIs, 20,019 sample-cells | 0.14 s / 111,116 KiB | 0.18 s / 123,156 KiB | 8 / 0 | 60,464,172 |
| joint include/exclude set | 182 include-only, 4 exclude-only, 0 both | 0.13 s / 402,872 KiB | 0.52 s / 203,380 KiB | 8 / 0 | 36,588,673 |

Thus the matched single-run wall ratios are 3.46x, 2.00x, at least one timer quantum for the
absent query, 1.29x, and 4.00x, respectively. These timings describe this fixed warm-host run; they
are not a multi-replicate performance distribution. The collection path also uses more peak RSS
than the sequential naive baseline for dense junction and junction-set output. Its primary gains
are routing, one command over the collection, and avoiding irrelevant archive reads—not a general
memory reduction.

The routed reads are instrumented at the source archive reader rather than inferred from result
size. Dense/junction-set, sparse, and region operations read 2.759%, 1.063%, and 4.560% of aggregate
source bytes. The joint query decodes eight unique chunks rather than ten independent component
decodes, a 20% reduction. Raising the dense-junction minimum support from 0 to 211, just above its
catalogue upper bound of 210, prunes all eight archives and reads zero archive bytes while retaining
the explicit catalogue-presence state.

## Integrity and construction

The collection is an optional, rebuildable index, not a replacement archive or a change to the
v1 evidence contract. Version 2 authenticates its fixed header and section directory, stores a
digest for every decompressed section, binds archives to a stamped reference, and rejects duplicate
paths, files, and content identities. Inspection decodes and semantically audits every sidecar
section. Ordinary queries recheck the opened source file's size, nanosecond mtime/ctime, device, and
inode. `--verify-content` additionally rehashes source content.

The verified dense query accounts for exactly 1,326,036,691 identity bytes—the complete eight-file
source size once—and returns the unchanged 182-UMI result. It takes 0.20 s and 404,720 KiB versus
0.13 s without content verification. Full semantic inspection takes 0.15 s and 20,300 KiB; adding
content verification takes 0.22 s and 50,140 KiB. These very high hash rates reflect warm local
storage and must not be generalized to cold or remote filesystems.

The fresh eight-archive root builds in 0.76 s at 508,752 KiB. The four-archive base and extension
build in 0.45 s / 272,616 KiB and 0.47 s / 269,628 KiB. Because the current `.aie` format does not
carry a reusable whole-archive root commitment, exact source identity still requires reading all
compressed archive bytes at build time, although no molecule section is decompressed. Embedding an
archive root in a future format would remove that scan; the present result must not be described as
metadata-only construction.

## Reproduction and retained evidence

From the project root, regenerate the excluded run tree with
`gravlax-paper-scripts/scripts/206_benchmark_collection_index.sh`, then reduce it with:

```bash
python3 gravlax-paper-scripts/scripts/207_validate_collection_index.py \
  --run-dir runs/post-v1/collection-index-v2-r5 \
  --binary env/cargo-target-atlas-v2/release/aie \
  --gravlax-commit 6b92503cffbd48da8312bdd75c1f1f70d5ff2a27 \
  --out gravlax-paper-scripts/results/post-v1-collection-index-v2.json
```

The validator fixes the complete 114-file run layout, validates all JSON schemas and internal
accounting, checks the executable and four sidecar SHA-256 records, checks source sizes and recorded
BLAKE3 identities, parses all 24 GNU-time records, and commits a canonical digest over relative
path, byte size, and SHA-256 for every run file. Large archives, sidecars, query JSON, and timing
records remain excluded; their identities are in the compact result and excluded-artifact manifest.
