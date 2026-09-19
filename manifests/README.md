# Manifest conventions

- Paths are relative to the historical analysis root
  `${GRAVLAX_PROJECT_ROOT}` unless marked as URLs.
- `input-digests.sha256` is the original input checksum file created during the experiments.
- `publication-artifacts.tsv` contains SHA-256 checksums computed at snapshot time for the `.aie`
  files, function-preserving CRAM baselines, and generated annotations most directly supporting
  manuscript claims.
- `excluded-artifacts.tsv` records large, regenerable inputs and intermediates that are not stored
  in Git. A blank checksum means that the snapshot records identity and size but deliberately does
  not spend I/O hashing a public download or regenerable intermediate.
- `adaptive-evidence-frontier-artifacts.tsv` expands the rejected frontier pilot's excluded run
  tree into exact logical paths, byte sizes, SHA-256 identities, and experimental roles.
- The collection-index entries use the canonical relative-path/size/SHA-256 file-set digest defined
  in `results/post-v1-collection-index-v2.json`; source `.aie` identities are the BLAKE3 digests
  verified by the collection builder and recorded per archive in that result.
- `sparse-terminal-tail-artifacts.tsv` pins every input and excluded output used by the eight-donor
  tail pilot, the NTRK2 secondary query, and the final-binary donor-A reproducibility check. It also
  records the exact experimental executable; the experimental Gravlax source itself is identified
  by commits in the result record and is not duplicated here.
- `jshape-routing-gate-b-artifacts.tsv` pins the confirmatory Gate-B release binary, frozen panel,
  eight rooted input archives, complete r3 run tree, and canonical and repeated routed/fallback
  collection files. Gate-B r1, r2, and r3 are separately labeled preflight, diagnostic, and
  confirmatory in `excluded-artifacts.tsv`; only r3 supports the compact PASS.
- `software.tsv` records the tools used on the snapshot host. Rust crate versions are frozen by
  the tagged code repository's `Cargo.lock`.
