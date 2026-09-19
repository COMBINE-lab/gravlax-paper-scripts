#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-$P}
AIE=${AIE:-$PROJECT_ROOT/env/cargo-target-atlas-v2/release/aie}
BENCH_OUT=${BENCH_OUT:-$PROJECT_ROOT/runs/post-v1/collection-index-v2-r5}

if [[ ! -x "$AIE" ]]; then
  echo "missing executable: $AIE" >&2
  exit 2
fi
ids=(A B C D E F G H)
archives=(
  "$PROJECT_ROOT/runs/post-v1/sez-transcript-end-atlas-r1/donor-A/sez.called.aie"
  "$PROJECT_ROOT/runs/post-v1/sez-transcript-end-atlas-r1/donor-B/sez.called.aie"
  "$PROJECT_ROOT/runs/post-v1/sez-transcript-end-atlas-r1/donor-C/sez.called.aie"
  "$PROJECT_ROOT/runs/post-v1/sez-transcript-end-atlas-r1/donor-D/sez.called.aie"
  "$PROJECT_ROOT/runs/post-v1/sez-transcript-end-atlas-r1/donor-E/sez.called.aie"
  "$PROJECT_ROOT/runs/post-v1/sez-transcript-end-atlas-r1/donor-F/sez.called.aie"
  "$PROJECT_ROOT/runs/post-v1/sez-transcript-end-atlas-r1/donor-G/sez.called.aie"
  "$PROJECT_ROOT/runs/post-v1/sez-transcript-end-atlas-r1/donor-H/sez.called.aie"
)
digests=(
  029de28423c59a8b9ca283e73687a9501cd7533cfd508788b61b15e79b751463
  4feced152203593f22ddb5a13207dbcd7dd0be42e08d229c5209a9433105e474
  9067a7d94acc472085a882b139e512a8e28dd1b39f761f9b9d2b4c6c0a697076
  4f69add3177ebee40584aea2b4c1e8e40a81a1b64a32c455c27e502cd1e0b1f9
  3afe043952a1d1c1b4d831365d3c12b4bdbd2c040d09e4c5a99defde964d76e5
  8850c98d1453965fbce69fbb5cad17a71ee8714e5c1d2303dabff1c4e1e8ec55
  19e9672274ef24bb194926ac0fc66e3c3c1b788ece74f1e478e9e04f7a7ddee3
  d26731d109a6c0b74afe95e32f1471a290fa8e274f1e75c4ec7e55d70673375c
)

run_naive_kind() {
  local kind=$1
  local output_dir=$2
  for index in "${!ids[@]}"; do
    local id=${ids[$index]}
    local archive=${archives[$index]}
    case "$kind" in
      dense)
        "$AIE" query "$archive" junction chr9:129925157-129927194 --json --top 0 \
          > "$output_dir/$id-dense.json"
        ;;
      sparse)
        if ! "$AIE" query "$archive" junction chr9:129861033-129861542 --json --top 0 \
          > "$output_dir/$id-sparse.json" 2> "$output_dir/$id-sparse.stderr.txt"; then
          grep -F "not present in the archive" "$output_dir/$id-sparse.stderr.txt" >/dev/null
          printf '{"present":false,"umis":0,"cells":0}\n' \
            > "$output_dir/$id-sparse.json"
        fi
        ;;
      absent)
        if "$AIE" query "$archive" junction chr9:129861034-129861542 --json --top 0 \
          > "$output_dir/$id-absent.stdout.txt" \
          2> "$output_dir/$id-absent.stderr.txt"; then
          echo "expected absent junction unexpectedly resolved in sample $id" >&2
          exit 1
        fi
        grep -F "not present in the archive" "$output_dir/$id-absent.stderr.txt" >/dev/null
        printf '{"present":false,"umis":0,"cells":0}\n' \
          > "$output_dir/$id-absent.json"
        ;;
      region)
        "$AIE" query "$archive" region chr9:84800000-85040000 --json --top 0 \
          > "$output_dir/$id-region.json"
        ;;
      jset)
        "$AIE" query "$archive" jset \
          --include chr9:129925157-129927194 --exclude chr9:129861033-129861542 \
          --json --top 0 > "$output_dir/$id-jset.json"
        ;;
      *)
        echo "unknown naive query kind: $kind" >&2
        exit 2
        ;;
    esac
  done
}

if [[ ${1:-} == "__naive-kind" ]]; then
  [[ $# -eq 3 ]] || { echo "__naive-kind requires KIND OUTPUT_DIR" >&2; exit 2; }
  run_naive_kind "$2" "$3"
  exit 0
fi

if [[ -e "$BENCH_OUT" ]]; then
  echo "refusing to overwrite existing benchmark directory: $BENCH_OUT" >&2
  exit 2
fi
mkdir -p "$BENCH_OUT/naive"

root_args=()
reverse_args=()
base_args=()
extension_args=()
for index in "${!ids[@]}"; do
  pair="${ids[$index]}=${archives[$index]}"
  digest="${ids[$index]}=${digests[$index]}"
  root_args+=(--sample "$pair" --source-digest "$digest")
  if (( index < 4 )); then
    base_args+=(--sample "$pair" --source-digest "$digest")
  else
    extension_args+=(--sample "$pair" --source-digest "$digest")
  fi
done
for ((index=${#ids[@]} - 1; index >= 0; index--)); do
  reverse_args+=(
    --sample "${ids[$index]}=${archives[$index]}"
    --source-digest "${ids[$index]}=${digests[$index]}"
  )
done

/usr/bin/time -v -o "$BENCH_OUT/root-build.time.txt" \
  "$AIE" collection build "${root_args[@]}" --out "$BENCH_OUT/sez8-root.aicollection" \
  > "$BENCH_OUT/root-build.stdout.txt"
/usr/bin/time -v -o "$BENCH_OUT/reverse-build.time.txt" \
  "$AIE" collection build "${reverse_args[@]}" --out "$BENCH_OUT/sez8-reverse.aicollection" \
  > "$BENCH_OUT/reverse-build.stdout.txt"
cmp "$BENCH_OUT/sez8-root.aicollection" "$BENCH_OUT/sez8-reverse.aicollection"

/usr/bin/time -v -o "$BENCH_OUT/base-build.time.txt" \
  "$AIE" collection build "${base_args[@]}" --out "$BENCH_OUT/sez4-base.aicollection" \
  > "$BENCH_OUT/base-build.stdout.txt"
sha256sum "$BENCH_OUT/sez4-base.aicollection" > "$BENCH_OUT/base-before.sha256"
/usr/bin/time -v -o "$BENCH_OUT/extension-build.time.txt" \
  "$AIE" collection build --base "$BENCH_OUT/sez4-base.aicollection" \
  "${extension_args[@]}" --out "$BENCH_OUT/sez8-extension.aicollection" \
  > "$BENCH_OUT/extension-build.stdout.txt"
sha256sum "$BENCH_OUT/sez4-base.aicollection" > "$BENCH_OUT/base-after.sha256"
cmp "$BENCH_OUT/base-before.sha256" "$BENCH_OUT/base-after.sha256"

/usr/bin/time -v -o "$BENCH_OUT/root-inspect.time.txt" \
  "$AIE" collection inspect "$BENCH_OUT/sez8-root.aicollection" \
  > "$BENCH_OUT/root-inspect.json"
/usr/bin/time -v -o "$BENCH_OUT/chain-inspect.time.txt" \
  "$AIE" collection inspect "$BENCH_OUT/sez8-extension.aicollection" \
  > "$BENCH_OUT/chain-inspect.json"
/usr/bin/time -v -o "$BENCH_OUT/root-inspect-verify.time.txt" \
  "$AIE" collection inspect "$BENCH_OUT/sez8-root.aicollection" --verify-content \
  > "$BENCH_OUT/root-inspect-verify.json"

run_collection_queries() {
  local label=$1
  local collection=$2
  /usr/bin/time -v -o "$BENCH_OUT/$label-dense.time.txt" \
    "$AIE" collection junction "$collection" chr9:129925157-129927194 \
    --json --explain --top 0 > "$BENCH_OUT/$label-dense.json"
  /usr/bin/time -v -o "$BENCH_OUT/$label-sparse.time.txt" \
    "$AIE" collection junction "$collection" chr9:129861033-129861542 \
    --json --explain --top 0 > "$BENCH_OUT/$label-sparse.json"
  /usr/bin/time -v -o "$BENCH_OUT/$label-absent.time.txt" \
    "$AIE" collection junction "$collection" chr9:129861034-129861542 \
    --json --explain --top 0 > "$BENCH_OUT/$label-absent.json"
  /usr/bin/time -v -o "$BENCH_OUT/$label-region.time.txt" \
    "$AIE" collection region "$collection" chr9:84800000-85040000 \
    --json --explain --top 0 > "$BENCH_OUT/$label-region.json"
  /usr/bin/time -v -o "$BENCH_OUT/$label-jset.time.txt" \
    "$AIE" collection jset "$collection" \
    --include chr9:129925157-129927194 --exclude chr9:129861033-129861542 \
    --json --explain --top 0 > "$BENCH_OUT/$label-jset.json"
}

run_collection_queries root "$BENCH_OUT/sez8-root.aicollection"
run_collection_queries chain "$BENCH_OUT/sez8-extension.aicollection"
/usr/bin/time -v -o "$BENCH_OUT/root-support-pruned.time.txt" \
  "$AIE" collection junction "$BENCH_OUT/sez8-root.aicollection" \
  chr9:129925157-129927194 --min-support 211 --json --explain --top 0 \
  > "$BENCH_OUT/root-support-pruned.json"
/usr/bin/time -v -o "$BENCH_OUT/root-verify-dense.time.txt" \
  "$AIE" collection junction "$BENCH_OUT/sez8-root.aicollection" \
  chr9:129925157-129927194 --verify-content --json --explain --top 0 \
  > "$BENCH_OUT/root-verify-dense.json"

export PROJECT_ROOT AIE
for kind in dense sparse absent region jset; do
  /usr/bin/time -v -o "$BENCH_OUT/naive-$kind.time.txt" \
    "$0" __naive-kind "$kind" "$BENCH_OUT/naive"
done

sha256sum "$AIE" "$BENCH_OUT"/*.aicollection > "$BENCH_OUT/artifacts.sha256"
printf 'PASS\n'
