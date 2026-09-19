# Molecular splice graph v1

## Verdict

The exact per-archive primitive passes its frozen gate. Gravlax commit
`2d4d51e0ae47bbfa402d3d39ba797f588a067e49` adds `aie query <archive>
splice-graph <locus>`. It returns separate strand-aware junction graphs and
exact UMI-class path fragments at cell, group, or bulk scope. It does not
reconstruct transcripts, infer missing exons, or emit population statistics.

The synthetic integration gate covers repeated representatives, multi-edge
co-support, opposite strands, directed reverse-strand edges, cell/group/bulk
agreement, edge/path conservation, zero thresholds, malformed coordinates,
invalid scopes, and hard overflow. All 119 workspace tests pass.

## Real-locus result

The registered real check uses the donor-G targeted FNBP1 archive from the
independent-donor experiment, its frozen seven-group map, and the locked
`chr9:129907500-129926000` window. The default thresholds select 21 catalogue
junctions. Their 21 independent postings all resolve to one archive chunk,
which is decoded once.

The query returns 37 nodes, 22 strand-specific edges, and 24 exact path
patterns representing 51 strand-path UMI counts. Twenty-one patterns have one
selected edge (45 UMIs); three have two selected edges (six UMIs). All three
multi-edge fragments are on the FNBP1-consistent reverse strand. The strongest
has four UMIs in four nuclei: two oligodendrocyte/OPC and two astrocyte. The
other two have one UMI each in microglia/immune and neuronal nuclei. These are
capability observations in one donor, not differential-splicing results.

Across 20 warm release-process invocations, median wall time is 0.004785 s
(range 0.004249--0.005645 s) and maximum peak RSS is 13,204 KiB. JSON is
byte-identical across repetitions. Both aggregate and per-group edge UMI
counts equal the sum of retained path-fragment counts containing each edge.

## Interpretation and next layer

A path fragment means only that one archive UMI class positively co-supports
the selected junction set on one strand. A one-edge fragment does not prove
that no other exon was present, and representative reduction does not retain
every original read-level co-occurrence. Multimappers follow the same archived
anchor-placement rule as existing junction/event queries. Catalogue support
is strand-combined; returned nodes, edges, paths, and counts are separated by
strand.

The next implementation should operate above this exact primitive: accept an
explicit sample/design table, retain one path/edge count row per donor or
sample, and only then evaluate replicate-aware beta-binomial or
Dirichlet-multinomial contrasts. Cells and molecules must never be substituted
for biological replicates.
