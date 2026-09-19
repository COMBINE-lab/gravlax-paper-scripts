#!/usr/bin/env bash
# Run the fixed D5 FNBP1 event with the pre-locked broad-cell labels and unchanged executable.
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 ROOT CONFIRMATION_RUN LOCALIZATION_RUN" >&2
  exit 2
fi

ROOT=$(realpath "$1")
CONFIRMATION=$(realpath "$2")
OUT=$(realpath -m "$3")
mkdir -p "$OUT"

AIE=${AIE:-$ROOT/env/cargo-target-features/release/aie}
PYTHON=${PYTHON:-python3}
MANIFEST=$ROOT/gravlax-paper-scripts/experiments/manifests/fnbp1-broad-cell-localization.json

if [[ ! -x "$AIE" ]]; then
  echo "missing existing executable: $AIE" >&2
  exit 1
fi

# The group builder is deliberately completed before the event query starts.
uv run --with 'h5py==3.16.0' --with 'numpy==2.5.2' --python 3.12 \
  "$ROOT/gravlax-paper-scripts/scripts/186_build_fnbp1_broad_cell_groups.py" \
  --h5 "$CONFIRMATION/filtered_feature_bc_matrix.h5" \
  --manifest "$MANIFEST" --out "$OUT"

GROUP_SHA=$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["group_map_sha256"])' "$OUT/groups-summary.json")
echo "$GROUP_SHA  $OUT/broad-cell-groups.tsv" | sha256sum -c -

# This targeted archive contains only the 6,584 called nuclei represented by local reads.  The
# existing executable requires every scope barcode to occur in its dictionary, so enumerate that
# dictionary with a whole-chromosome region query and restrict, without relabeling, the frozen map.
"$AIE" query "$CONFIRMATION/fnbp1.called.aie" region chr9:0-138394717 \
  --top 100000 --json > "$OUT/archive-cells.json"
"$PYTHON" "$ROOT/gravlax-paper-scripts/scripts/188_filter_fnbp1_archive_groups.py" \
  --full-map "$OUT/broad-cell-groups.tsv" --archive-cells "$OUT/archive-cells.json" \
  --out "$OUT/broad-cell-groups.archive.tsv"
ARCHIVE_GROUP_SHA=$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["archive_group_map_sha256"])' "$OUT/archive-scope-summary.json")
echo "$ARCHIVE_GROUP_SHA  $OUT/broad-cell-groups.archive.tsv" | sha256sum -c -

/usr/bin/time -v -o "$OUT/query.time.txt" \
  "$AIE" query "$CONFIRMATION/fnbp1.called.aie" events chr9:129907500-129926000 \
  --min-support 2 --min-informative 1 --max-events 100000 \
  --gtf "$CONFIRMATION/gencode.v49.aic" \
  --groups "$OUT/broad-cell-groups.archive.tsv" --agg group --json \
  > "$OUT/events-by-broad-cell.json" 2> "$OUT/query.stderr.txt"

sha256sum "$AIE" "$CONFIRMATION/fnbp1.called.aie" "$OUT/broad-cell-groups.tsv" \
  "$OUT/broad-cell-groups.archive.tsv" \
  > "$OUT/execution-inputs.sha256"

"$PYTHON" "$ROOT/gravlax-paper-scripts/scripts/189_summarize_fnbp1_broad_cell_localization.py" \
  --run "$OUT" --manifest "$MANIFEST" \
  --confirmation-summary "$CONFIRMATION/summary.json" \
  --executable "$AIE" \
  --out-json "$ROOT/gravlax-paper-scripts/results/post-v1-fnbp1-broad-cell-localization.json" \
  --out-tsv "$ROOT/gravlax-paper-scripts/results/post-v1-fnbp1-broad-cell-localization.tsv"
