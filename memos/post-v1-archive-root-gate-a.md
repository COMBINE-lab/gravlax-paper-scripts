# Content-addressed archive identity: Gate-A result

## Verdict

**PASS.** Gravlax archive v2 commits to the exact encoded section content with a
domain-separated BLAKE3 archive root. Across the eight SEZ archives, the commitment adds 355,840
bytes to 1,326,036,691 source bytes (**0.026835%**) and eliminates the complete-file identity scan
when constructing a federated collection. The rooted build reads zero identity-content bytes and
10,341,234 total source bytes (**0.779860%**) for the declared six metadata/index sections, versus
1,326,036,691 identity-content bytes and 1,336,022,077 total bytes for the v1 control.

This licenses a content-addressed archive identity and metadata/index-only collection construction
after one-time sealing. It does not license a publisher-authentication claim, zero-I/O collection
construction, or whole-file payload verification on every ordinary query.

## Representation and integrity result

Each physical-order directory entry carries a 32-byte BLAKE3 digest of its compressed payload; the
footer carries a 32-byte root over the versioned header, directory offset, and exact directory
bytes. The eight SEZ archives contain 11,112 sections, so the observed premium is exactly 32 bytes
per section plus 32 bytes per archive. D0, used only as the replay guard, adds 36,096 bytes over
1,127 sections.

The locked reducer independently reproduced every root and compressed-payload commitment and
required byte-identical compressed payloads before and after sealing. Repeated seals were
byte-identical. Fresh eight-archive and reversed-input collection builds were byte-identical within
each identity arm. The 30-case adversarial panel accepted all six legal boundary cases and rejected
all 24 malformed/corrupted cases without partial scientific stdout, including changed roots or
directory fields, swapped/corrupt payloads, overlap/gaps/truncation, allocation and compression
bombs, selected-payload corruption on an ordinary read, and unselected corruption under full
verification.

Ordinary reads validate the directory/root and each selected compressed payload. An unselected
corrupt payload is intentionally outside that lazy guarantee; `--verify-content` checks every
payload. The root is a content commitment, not a signature or provenance attestation.

## Exact scientific behavior

Gene and Velocity replay are exact across v1/v2 and streaming/eager execution. Fresh-root and
A--D-base-plus-E--H-extension collections agree exactly, as do v1 and v2, for all fixed queries:

| query | exact result |
|---|---:|
| dense FNBP1 junction | 182 UMIs, 180 cells |
| sparse neighboring junction | 4 UMIs, 4 cells |
| absent neighboring junction | 0 UMIs, 0 cells |
| NTRK2 region | 74,848 molecules, 74,352 UMIs, 20,019 sample-cells |
| include-dense/exclude-sparse junction set | 182 include-only, 4 exclude-only, 0 both |

Thus the representation and identity change does not change any tested molecular, count, or matrix
answer. The Gate-A result is an engineering/integrity result; it adds no new biological evidence.

## Resource result

All measurements use eight alternating warm paired blocks with 24 threads on the same host.

| operation | v1 median | v2 median | v2/v1 |
|---|---:|---:|---:|
| eight-archive collection build wall | 0.77 s | 0.71 s | 0.9221 |
| collection build peak RSS | 498,510 KiB | 477,336 KiB | 0.9575 |
| D0 streaming replay wall | 3.260 s | 3.245 s | 0.9954 |
| D0 streaming replay peak RSS | 1,315,786 KiB | 1,332,758 KiB | 1.0129 |
| dense query wall / RSS | 0.13 s / 401,148 KiB | 0.14 s / 403,390 KiB | 1.0769 / 1.0056 |
| sparse query wall / RSS | 0.12 s / 178,628 KiB | 0.12 s / 180,628 KiB | 1.0000 / 1.0112 |
| junction-set wall / RSS | 0.13 s / 403,032 KiB | 0.14 s / 401,974 KiB | 1.0769 / 0.9974 |

The rooted build performs 129.19-fold less total source I/O than the legacy build and is 7.8%
faster in this warm comparison. Replay is effectively neutral; its 1.29% RSS increase is 16,972
KiB and remains inside the frozen 1.02 ratio. Dense and junction-set wall differences are one
10-ms timer quantum and remain inside the 1.10 non-regression limit. These results support a build
and identity claim, not a query-speed claim.

Primary sealing of the eight archives took 3.93 s in aggregate on the warm local host, with
202,324 KiB maximum RSS. Sealing copies and verifies compressed frames without recompression. It is
a one-time construction cost and is excluded from replay/query comparisons.

## Locked provenance and run history

The confirmatory run is `runs/post-v1/archive-root-gate-a-r3`, using clean Gravlax commit
`fbd48e536b39729632f189a74e7ece891166ecd1`, clean reproducibility commit
`df1e363d397c740300c97ec2b1f04fe48cd9eb6b`, and release binary SHA-256
`f0fb4c0d6a12bc9b3877bbab09168e5754939fe97ad4a363db51fa5c91c06211`. Its complete manifest
binds 813 files and 6,896,522,383 bytes with SHA-256
`d64a953320a56d53e0caaa11d81c91bffbb55a2067027466201f161a49def64d`.
An independent second reduction produced the same result after normalizing only the requested
inventory-output path, and its 59-row excluded-artifact inventory was byte-identical (SHA-256
`60ac744b559e42918ebab40570bde981600cfd01736e2f744203958f6e87dada`).

Two earlier runs are retained only as harness-development provenance:

- r1 stopped at the first paired collection build because its output parent was not precreated. It
  has no completion record or driver artifact manifest and supplies no confirmatory measurement.
- r2 completed the workload, but was used diagnostically to find and correct reducer assumptions
  about manifest path ordering, intentional replay diagnostics, repeated replay resource parsing,
  and identity-dependent collection metadata size. Its legacy manifest binds 813 files and
  6,896,522,424 bytes (`e6d5758a...ba6f7`), but r2 is not the reported arm.

The large r1/r2/r3 trees remain excluded. The compact PASS is
`results/post-v1-archive-root-gate-a.json`; exact r3 artifact identities are in
`manifests/archive-root-gate-a-artifacts.tsv`.

## Decision boundary

Archive v1 remains readable and can be sealed without recompression. Gate A supports making the
root-committed representation the authoritative identity for new archives and federated
collections. Gate B remains unopened in this commit: its 96-junction panel must be materialized,
hashed, and committed together with the serialization addendum before any local-shape-routing code
or measurement begins.
