# Authenticated archive root, then local-shape routing

**Frozen:** 2026-09-01, before either implementation or candidate measurement.

**Gravlax base:** `1c1b4fd0c5058d06081c6c3ea3db7c453fa6ef97`.

**Reproducibility base:** `7a6818377476b2fad1b3b4505bcd82017faca922`.

**Order:** Gate A must pass every correctness and integrity requirement before Gate B begins.

Archive v1 remains readable for the tagged snapshot, but only a root-committed v2 archive may use
the fast identity path. Here, *authenticated* means a domain-separated cryptographic content
commitment. It is not a publisher signature and does not establish provenance by itself.

## Gate A: archive-internal authenticated root commitment

### Hypothesis and fixed representation

A v2 `.aie` can expose an exact encoded-content identity from its header and directory alone. This
should eliminate the complete compressed-file identity scan during collection construction with
negligible storage and no replay or query regression.

Each physical-order v2 directory entry adds a 32-byte BLAKE3 digest of the compressed section
payload to the existing canonical fields `(name, offset, raw_len, compressed_len)`. The fixed
44-byte footer is:

```text
directory_offset: u64 little endian
archive_root:     32 bytes
magic:            "AIED"
```

The archive root is:

```text
BLAKE3(domain_tag || "AIE0" || version_le || directory_offset_le || exact_directory_bytes)
```

The domain tag is a fixed byte string in the format implementation. It may not be changed after
this gate is frozen. The root commits to the exact directory bytes, and the directory commits to
every compressed payload. It therefore identifies exact encoded section bytes, not merely
decompressed semantics. Writers calculate each payload digest while writing; a v1-to-v2 migration
copies compressed frames without recompression and writes a new directory/footer atomically.

The v2 reader requires canonical unique names, directory order equal to physical section order,
strictly increasing nonoverlapping contiguous extents, exactly one terminator, the directory
immediately after the section area, the footer at end of file, and no gaps or trailing bytes.
Counts and lengths are bounded before allocation. A selected section's compressed digest is
verified before decompression. Collection manifests bind the archive root, byte length, reference
signature, and required directory records; path and filesystem-stat guards remain.

Normal queries validate the root/directory and every payload they select. An unselected corrupted
payload need not fail an ordinary query, but must fail `--verify-content`. This lazy guarantee must
not be described as whole-file verification.

### Immutable inputs and baselines

- The eight SEZ donor archives A--H under
  `runs/post-v1/sez-transcript-end-atlas-r1`, totaling 1,326,036,691 bytes. Their paths, BLAKE3
  identities, chunk counts, and shared reference digest are frozen in
  `results/post-v1-collection-index-v2.json`.
- D0 `pbmc_1k_v3`, with the existing v49 Gene and Velocyto replay inputs and called-cell list.
- Synthetic archives covering zero, one, and many sections, an unknown optional section, and the
  maximum legal name and declared lengths.
- Control: v1 plus a complete-file BLAKE3 identity scan.
- Treatment: v2 plus root-based identity, using the same binary, host, threads, compression level,
  and inputs within each comparison.

Construct the fresh eight-sample collection and the existing A--D base plus E--H immutable
extension. Performance uses eight warm paired blocks, with control/treatment order alternating by
block. A separately labeled first-touch observation may be retained, but is descriptive and is
not a cold-cache claim.

### Exactness and integrity gates

All are mandatory.

1. Every migrated v2 archive has the same ordered section names, raw lengths, compressed lengths,
   and compressed-payload SHA-256 and BLAKE3 values as its v1 source.
2. Eager molecule decoding, streaming replay, and D0 v49 Gene and Velocyto outputs are exact
   between v1 and v2.
3. The following collection queries are exact after removing timing and plan fields:
   - dense junction `chr9:129925157-129927194`;
   - sparse junction `chr9:129861033-129861542`;
   - absent junction `chr9:129861034-129861542`;
   - region `chr9:84800000-85040000`; and
   - junction set with dense inclusion and sparse exclusion.
4. Fresh-root and base-plus-extension answers agree exactly. Two migrations and reversed-input
   collection builds are byte-identical. An independent root implementation reproduces every
   emitted root.
5. Tests reject changes to magic, version, root, any directory byte, entry order/name/offset/
   length/digest, inline header, terminator, footer, payload, or decompressed length. They also
   reject swapped payloads, overlap, gaps, truncation, trailing bytes, future versions, and
   allocation or compression bombs.
6. Corruption of a selected payload fails an ordinary read. Corruption of any payload fails full
   content verification. All failures occur without returning partial scientific output.

### Resource gates

- Added bytes are at most 0.10% of the aggregate archives and at most 64 bytes per section plus
  256 bytes per archive.
- Root-based collection construction reports zero `identity_content_bytes_read`; all source bytes
  read are at most 2.0% of aggregate archive bytes, with no `cN`, `coc.*`, `patterns`, or `edges`
  payload read.
- Across the paired collection builds, treatment median wall time is at most 1.10 times control,
  and peak RSS is at most 1.05 times control and at most 600 MiB.
- D0 streaming-replay median wall time is at most 1.05 times control and peak RSS at most 1.02
  times control.
- Dense, sparse, and junction-set collection commands each have median wall time and peak RSS at
  most 1.10 times their v1-rootless controls.
- Migration wall time, CPU, peak RSS, bytes read, and bytes written are recorded separately. They
  are one-time construction costs, not query costs.

### Decision and stop rules

An identity collision for unequal encoded section bytes, any exactness failure, a missed required
corruption, acceptance of a noncanonical layout, or unbounded allocation is a hard **STOP**; Gate B
does not begin. If integrity passes but the size or I/O gate fails, retain v1 reading and do not
make v2 the writer default. Lack of a warm-cache speedup is not failure when the I/O gate passes,
but a regression beyond a frozen limit stops promotion.

A pass licenses *metadata/index-only collection identity after one-time archive construction*.
It does not license *publisher-authenticated*, *whole-file verified on every query*, or
*zero-I/O collection construction*.

## Gate B: derived junction-to-local-shape routing

### Freeze boundary, hypothesis, and fixed representation

The 96-row query panel described below is materialized and its SHA-256 is committed **after Gate A
passes but before any Gate B implementation or result is opened**. No panel row may be changed in
response to Gate B behavior.

The frozen eight-donor archives contain 33.897 MB of compressed `shapes` sections. The hypothesis
is that loading all of those dictionaries accounts for most excess bytes and memory in exact point
and junction-set federation. A source-root-bound derived route should remove those reads while the
molecule chunks remain authoritative.

For each source archive, derive from its `shapes` dictionary only a seekable map:

```text
intron_span -> sorted unique (local_shape_id, donor_offset, acceptor_offset)
```

For a requested `(donor, acceptor)`, load only the bucket containing
`intron_span = acceptor - donor`. A decoded representative `(position, shape_id)` matches only if
the shape id is routed and checked arithmetic proves both
`position + donor_offset == donor` and `position + acceptor_offset == acceptor`. Every stored chain
representative and multimapper placement is tested. Shapes sharing an intron span may be
over-routed, but no candidate can become a hit without the exact coordinate checks.

The map contains no cells, UMIs, molecule ordinals, annotations, groups, or biological result. It
is optional derived `.aicollection` state, partitioned into deterministic bounded sections so an
exact span never loads the complete map. Its collection commitment binds it to the source archive
root and exact `shapes`-entry digest. Construction reads no `cN` molecule payload. Collections
without the route use the existing full-shape fallback.

### Inputs, prospective panel, and controls

- Use the Gate-A root-committed SEZ A--H archives and both the fresh-root and A--D/E--H chain.
- Retain the dense, sparse, absent, region, and junction-set queries above. Region is a negative
  control because it does not consult shapes.
- Materialize 96 exact junctions from the frozen collection catalogue: 24 from each archive-
  presence stratum `{1, 2--3, 4--7, 8}`. Within a stratum, exclude malformed or zero-span rows,
  compute
  `BLAKE3("gravlax-jshape-panel-v1\0" || chrom || donor_le || acceptor_le)`, and select the 24
  smallest hashes, breaking a hash tie by `(chrom, donor, acceptor)`. Commit the resulting ordered
  TSV and SHA-256 before Gate B code or measurement.
- Control: a root-committed collection built by the Gate-B binary with routing disabled.
- Treatment: the same binary and inputs with routing enabled. Use `--json --top 0`, the same host
  and threads, and eight warm alternating control/treatment blocks. Report internal phase timings
  as well as process-level resource records.

### Exactness and integrity gates

1. Scientific JSON is exact between route and fallback for all 96 panel junctions, the three known
   point queries, the junction set, and region; root and chain agree exactly.
2. Synthetic fixtures cover a single-block shape; multiple introns; repeated equal intron spans at
   distinct offsets in one shape; distinct shapes sharing a span; forward and reverse records;
   both representatives; multimappers; a UMI class crossing chunks; an absent coordinate; and
   checked-add overflow.
3. With routing absent or disabled, behavior and answers retain the Gate-A fallback contract.
   Repeated builds and reversed archive input produce byte-identical route state.
4. Parsing and inspection reject an unknown route version; incorrect source-root or shape-digest
   binding; an out-of-range shape id; reversed or inconsistent offsets; unsorted or duplicate
   tuples; wrong archive ordinal or span bucket; count overflow; truncation; checksum mismatch;
   and trailing bytes.
5. Full route inspection independently reconstructs the route from each source `shapes` section
   and requires an exact match. A query validates each reported hit against the routed shape and
   both offsets; the route never supplies a count.

### Size, construction, I/O, runtime, and memory gates

- The route premium is at most 3.0% of aggregate source bytes; the complete collection is at most
  4.0%, retaining the earlier 5% outer ceiling. Compressed route bytes are also at most 1.25 times
  the aggregate compressed `shapes` bytes.
- Construction reads/decompresses only declared dictionaries, indexes, and `shapes`; it reads zero
  `cN`, `coc.*`, `patterns`, or `edges` payloads. It takes at most 2.0 times root-only build wall
  time and at most 768 MiB peak RSS.
- Instrument collection-sidecar and source-archive bytes separately and report their sum.
- For dense, sparse, and junction-set queries, source-archive bytes are at most 20% of their frozen
  fallbacks: 36,588,673, 14,100,036, and 36,588,673 bytes, respectively.
- Across the 96-row panel, the median source-byte ratio is at most 0.25, p95 at most 0.50, and no
  query exceeds 1.05. Including route-section reads, the median total-logical-byte ratio is at
  most 0.50 and p95 at most 0.75.
- The absent query remains at zero source-payload bytes. Region bytes and wall time may not regress
  by more than 5%.
- Median aggregate panel wall time is at most 0.75 times fallback. The known dense and junction-set
  queries are each at most 0.80 times fallback; sparse is at most 0.85 times fallback.
- Dense and junction-set peak RSS are each at most 0.65 times fallback and at most 256 MiB. Sparse
  peak RSS is at most 0.75 times fallback and at most 128 MiB.
- Exact cell/class deduplication may not be weakened to meet any resource gate.

### Decision, stop rules, and failure interpretation

Any exactness or integrity failure is a hard **STOP** and the derived route is discarded. Missing
either size or total-I/O gate prevents promotion even if source bytes fall. If I/O passes but wall
or RSS does not, retain only an exact opt-in prototype and attribute remaining time to archive
open, chunk decode, cell-of-class lookup, aggregation, and output; make no speed claim. If fixed
loci pass but the prospective 96-row panel fails, call the result locus-specific and stop.

Region neutrality is expected. A clean pass licenses exact point-junction and junction-set
acceleration only; it says nothing about replay, region-query speed, or biological validity. After
an I/O-pass/runtime-fail, the next candidate is query-local class-to-cell reduction rather than a
larger shape index.

## Required reproducibility record

Retain the frozen panel TSV and digest; benchmark drivers and fail-closed validators; compact
Gate-A and Gate-B JSON summaries; decision memos; exact code, scripts, binary, Rust, zstd, BLAKE3,
host, and thread identities; per-command stdout, stderr, and resource records; per-section byte
manifests; a corruption-fixture manifest; and excluded large v2 archive/collection inventories
with logical path, byte size, checksum, provenance, and regeneration command.

Validators must freeze the complete run-file set and its canonical artifact-manifest digest, as in
the collection-v2 campaign. Result fields are added to the machine-readable gate manifest only
after a locked validator has produced them.
