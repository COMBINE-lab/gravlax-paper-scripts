# Federated sparse cohort output validation

## Verdict

**PASS.** The versioned sparse cohort bundle reconstructs the matched dense JSON exactly for the
six-archive, called-cell chromosome-17 event query. The validator compared the complete event
catalogue, per-sample catalogue presence, every total and group count fact, explicit sample/scope
dimensions, derived informative counts and usage fractions, and the logical zeros implied by
missing sparse rows.

The validated result contains 543 events across six samples: 3,258 event–sample observations,
2,383 positive catalogue-presence facts, and 6,516 nonzero total/group count rows. No dense row was
missing, duplicated, reordered, or numerically different.

## Representation result

The dense JSON is 3,484,594 bytes. The three compressed sparse tables total 72,741 bytes and their
metadata is 4,229 bytes, for a 76,970-byte self-describing bundle: **45.27× smaller**, or 2.21% of
the dense JSON size.

The matched sparse command took 3.92 seconds and 2,473,920 KiB maximum RSS; the dense JSON command
took 3.96 seconds and 2,452,780 KiB. The 0.99× wall and 1.009× RSS ratios are effectively neutral.
The arms were sequential and cache order was not counterbalanced, so these figures document
feasibility rather than a speed claim. The main value is a compact, streamable representation
proportional to emitted nonzero facts, with unambiguous zero reconstruction.

## Exactness checks

`scripts/205_validate_federated_sparse_output.py` fails closed unless all of the following hold:

- dense and sparse commands are identical before their output-selection flag and both exit zero;
- the sparse event order and all coordinate/type/annotation fields equal the dense events;
- the presence table equals exactly the dense rows with `present=true`;
- nonzero total, group, or bulk facts equal the dense rows in deterministic order;
- omitted count facts are all and only exact logical zeros over declared event/sample/scope
  dimensions;
- `informative_umis` and `usage_fraction` recompute from the stored primitive counts;
- metadata dimensions, planning statistics, thresholds, filenames, row counts, and compressed byte
  accounting agree with the dense JSON and decoded tables; and
- the metadata printed to stdout is identical to the installed `metadata.json`.

The final bundle uses the production `annotation_genes_json` event column, so embedded delimiters
cannot make gene IDs and names pair ambiguously. The validator retains a compatibility reader for
the earlier split-column pilot layout, but the frozen result and all reported bytes use the final
JSON-annotation layout.

## Reproduction

From the project root:

```bash
python3 gravlax-paper-scripts/scripts/205_validate_federated_sparse_output.py \
  --dense-json runs/post-v1/federated-sparse-output-r2-final/six-called.chr17.json \
  --sparse-dir runs/post-v1/federated-sparse-output-r2-final/six-called.chr17 \
  --dense-time runs/post-v1/federated-sparse-output-r2-final/six-called.chr17-json.time.txt \
  --sparse-time runs/post-v1/federated-sparse-output-r2-final/six-called.chr17-sparse.time.txt \
  --metadata-stdout runs/post-v1/federated-sparse-output-r2-final/six-called.chr17.metadata.stdout.json \
  --binary env/cargo-target-atlas-v2/release/aie \
  --gravlax-commit 6b92503cffbd48da8312bdd75c1f1f70d5ff2a27 \
  --out gravlax-paper-scripts/results/post-v1-federated-sparse-output-validation.json
```

The canonical run uses clean Gravlax commit
`6b92503cffbd48da8312bdd75c1f1f70d5ff2a27` and binary SHA-256
`3c5faaf38d8b9ea14001a8561b78b026452804d454b16a5f2e11dbc089605a4d`. Large dense and sparse query
artifacts remain outside Git; their compact identities and measurements are recorded in
`results/post-v1-federated-sparse-output-validation.json`.
