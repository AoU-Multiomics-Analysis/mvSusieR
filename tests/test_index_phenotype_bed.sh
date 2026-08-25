#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

test -x scripts/index_phenotype_bed.sh

fixture_dir="$(mktemp -d)"
trap 'rm -rf "$fixture_dir"' EXIT

Rscript tests/fixtures/trans_window/generate_index_fixture.R "$fixture_dir"
bash scripts/index_phenotype_bed.sh \
  "$fixture_dir/unsorted_phenotypes.tsv" \
  expression \
  "$fixture_dir/indexed" \
  2

indexed="$fixture_dir/indexed/phenotypes.bed.gz"
lookup="$fixture_dir/indexed/phenotype_lookup.tsv.gz"
qc="$fixture_dir/indexed/index_qc.tsv"

test -s "$indexed"
test -s "$indexed.tbi"
test -s "$lookup"
test -s "$qc"

test "$(tabix "$indexed" chr1:101-300 | wc -l | tr -d ' ')" -eq 3
test "$(tabix "$indexed" chr1:101-300 | cut -f4 | sort -u | wc -l | tr -d ' ')" -eq 3

gzip -cd "$indexed" | head -n 1 | \
  grep -Fx $'#chrom\tstart\tend\tphenotype_id\tsample_1\tsample_2\tsample_3'
gzip -cd "$lookup" | head -n 1 | \
  grep -Fx $'phenotype_id\tchrom\tstart\tend'

actual_order="$(gzip -cd "$indexed" | tail -n +2 | cut -f4 | paste -sd, -)"
test "$actual_order" = "feature_a,feature_b,feature_d,feature_c,feature_e"

expected_qc_header=$'modality\tphenotype_rows\tsample_columns\tinput_bytes\tbgzf_bytes\tlookup_bytes\telapsed_seconds'
test "$(head -n 1 "$qc")" = "$expected_qc_header"
test "$(tail -n 1 "$qc" | cut -f1-3)" = $'expression\t5\t3'

echo "Phenotype index artifact tests passed"
