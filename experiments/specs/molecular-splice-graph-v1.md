# Molecular splice graph v1 implementation gate

The first splice-graph layer is an exact per-archive evidence primitive.  A graph edge is a
strand-aware archived junction.  A molecular path fragment is the sorted set of selected junction
edges co-supported by one archive UMI class on that strand.  Counts are class-deduplicated and may
be restricted using the existing strict cell/group scopes.

The word *fragment* is binding.  Gravlax retains representative molecule evidence, not complete
transcript sequences, so a multi-edge fragment is positive co-support while a one-edge fragment
does not imply that no additional exon was present.  The command must not reconstruct or label a
full isoform.  Multimapper anchor placements remain unresolved, matching the current junction and
event query semantics.

The initial implementation returns versioned JSON with nodes, strand-aware edges, exact edge
counts, exact path-fragment counts, per-group reductions where requested, and planning metadata.
`--max-paths` is a hard guard rather than a truncation control.  Synthetic integration tests must
cover multi-edge co-support, UMI-class deduplication, strand separation, cell/group/bulk agreement,
edge/path conservation, strict input rejection, and unchanged legacy query output.

Population inference is deliberately the next layer.  A future cohort command will read an
explicit sample/design table and preserve one count row per donor or sample before applying a
beta-binomial or Dirichlet-multinomial model.  This v1 query must emit no p-value and must never
treat cells or molecules as biological replicates.
