#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

image="${TRANSQTL_MAPPABILITY_TEST_IMAGE:-transqtl-mappability-test}"
if [ -z "${TRANSQTL_MAPPABILITY_TEST_IMAGE:-}" ]; then
  docker build \
    --tag "$image" \
    --file envs/transqtl-mappability.Dockerfile \
    .
fi

docker run --rm "$image" genmap --version
docker run --rm "$image" bedtools --version
docker run --rm --entrypoint /bin/bash "$image" \
  -c 'command -v genmap && command -v bedtools'

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

{
  printf '>chrTiny\n'
  printf 'ACGT%.0s' {1..80}
  printf '\n'
} > "$tmpdir/reference.fa"

docker run --rm \
  --volume "$tmpdir:/data" \
  "$image" \
  /bin/bash -c '
    set -euo pipefail
    mkdir -p /data/index /data/output
    genmap index -F /data/reference.fa -I /data/index
    genmap map -K 20 -E 2 -I /data/index -O /data/output -t -w -bg
    bedgraph="$(find /data/output -type f -name "*.bedGraph" -print -quit)"
    test -n "$bedgraph"
    test -s "$bedgraph"
    awk '\''BEGIN { OFS = "\t" } $4 < 1 { print $1, $2, $3 }'\'' "$bedgraph" \
      | bedtools sort \
      | bedtools merge \
      > /data/low_mappability.bed
    test -s /data/low_mappability.bed
  '

echo "TransQTLMappability container smoke test passed"
