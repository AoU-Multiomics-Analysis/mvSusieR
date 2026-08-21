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

docker run --rm \
  "$image" \
  /bin/bash -c '
    set -euo pipefail
    tmpdir="$(mktemp -d)"
    trap '\''rm -rf "$tmpdir"'\'' EXIT
    {
      printf '\''>chrTinyA\n'\''
      printf '\''ACGT%.0s'\'' {1..80}
      printf '\''\n'\''
      printf '\''>chrTinyB\n'\''
      printf '\''ACGT%.0s'\'' {1..80}
      printf '\''\n'\''
    } > "$tmpdir/reference.fa"
    genmap index -F "$tmpdir/reference.fa" -I "$tmpdir/index"
    genmap map -K 20 -E 2 -I "$tmpdir/index" -O "$tmpdir/mappability" -t -w -bg
    bedgraph="$(find "$tmpdir" -type f \( -name "*.bedGraph" -o -name "*.bg" \) -print -quit)"
    test -n "$bedgraph"
    test -s "$bedgraph"
    printf '\''GenMap BEDGraph: %s\n'\'' "$bedgraph"
    head -n 5 "$bedgraph"
    low_records="$(awk '\''$4 < 1 { n++ } END { print n + 0 }'\'' "$bedgraph")"
    printf '\''low_records=%s\n'\'' "$low_records"
    test "$low_records" -gt 0
    awk '\''BEGIN { OFS = "\t" } $4 < 1 { print $1, $2, $3 }'\'' "$bedgraph" \
      | bedtools sort \
      | bedtools merge \
      > "$tmpdir/low_mappability.bed"
    test -s "$tmpdir/low_mappability.bed"
  '

echo "TransQTLMappability container smoke test passed"
