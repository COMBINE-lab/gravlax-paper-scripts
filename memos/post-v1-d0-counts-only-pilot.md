# D0 counts-only timing correction

Date: 2026-08-30  
Verdict: replace the stale speed headline; retain a replication caveat.

The v1 manuscript reported 150–370× against an end-to-end STARsolo workflow and correctly noted
that the denominator included sorted BAM output. The audit found a deeper provenance problem:
the frozen current D0 oracle log starts at 13:28:13 and finishes at 13:31:25 (about 192 s), while
the text used a roughly 16-minute D0 denominator. The compact record does not connect that older
denominator to the current prebuilt-index command. It is therefore unsuitable as the headline.

We reran STARsolo 2.7.11b at 24 threads with the same prebuilt v49 index, input reads, barcode
geometry, whitelist, correction/deduplication/filtering policies, and `Gene` feature, but with
`--outSAMtype None`. Five warm runs took 82.37–83.95 s (median 83.68 s; CV 0.75%). Every generated
file was byte-identical to the frozen oracle and no alignment file was emitted. Five warm Gravlax
replays took 2.36–2.58 s (median 2.42 s), giving a 34.6× ratio of medians. Median peak RSS was
34.93 GB versus 6.95 GB, a 5.03× reduction.

A separate no-BAM diagnostic retaining all four published STARsolo products (`Gene`, `GeneFull`,
`SJ`, and `Velocyto`) took 191.75 s. It too reproduced the frozen outputs byte-for-byte. The large
difference is mostly feature work—four-feature Solo counting took about 111 s, while Gene-only
counting took about 12 s—so future runtime tables must name the output capability explicitly.

The manuscript now uses 34.6× as the D0 preliminary counts-only result and says that the five
arms were sequential rather than randomized. A safe cold-cache observation, randomized
interleaving, and D1/D2' replication remain required before treating 34.6× as a multi-dataset
headline.

