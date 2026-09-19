# D5 FNBP1 broad-cell localization

## Scope and lock

This is a fixed-event, single-sample localization analysis, not a differential-splicing test.  The
event and D5 sample were fixed by the preceding prospective confirmation.  The complete
machine-readable label rule is
`experiments/manifests/fnbp1-broad-cell-localization.json`.  It was written before any FNBP1
evidence was queried by cell type.

Labels use only the filtered gene-expression matrix.  They cannot use splice evidence: FNBP1 is
absent from every marker module and is also removed from the library-size denominator.  Five
broad compartments are scored from two canonical marker submodules apiece: neuronal,
oligodendrocyte/OPC, astrocyte, microglia/immune, and vascular.  An explicit margin rule leaves
weak profiles unresolved and marks conflicting or unusually high-count multi-lineage profiles as
`ambiguous_or_doublet`.  The H5 `-1` GEM-group suffix is removed only when writing the group map,
because the existing archive dictionary contains the corresponding suffix-free barcodes.

## Exact event reduction

The frozen group map is passed to the existing, unchanged executable as:

```text
aie query fnbp1.called.aie events chr9:129907500-129926000 \
  --min-support 2 --min-informative 1 --max-events 100000 \
  --gtf gencode.v49.aic --groups broad-cell-groups.archive.tsv --agg group --json
```

Only the already locked cassette event is retained in the compact result.  Because the targeted
archive contains reads selected at the FNBP1 window, a called nucleus without a local read may be
absent from its cell dictionary.  The unchanged executable's whole-chromosome region query first
enumerates that dictionary, after which the frozen map is restricted without relabeling.  Both map
digests are recorded.  Event-bearing rates use the full H5 group sizes, not the archive-selected
cell count.

## Interpretation guardrail

Report localization, sparse support, and uncertainty honestly.  One donor cannot support a
sample-level contrast, and only 74 bulk informative molecules were observed previously.  A
compartment with no informative molecule is not evidence that the event is absent there.  The
analysis may motivate a donor-resolved fixed-event validation, but it cannot establish
cell-type-specific splicing or mechanism.
