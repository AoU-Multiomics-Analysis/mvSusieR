#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/prepare-window-genotype-header.XXXXXX")"
trap 'rm -rf "$tmp_dir"' EXIT

printf 'window_id\tchrom\tstart\tend\tmodality\tmolecular_trait_id\tp_value\n' \
  > "$tmp_dir/trans_window_associations.tsv"
printf 'w1\tchr1\t99\t200\texpression\texpr_1\t0.000001\n' \
  >> "$tmp_dir/trans_window_associations.tsv"
gzip -c "$tmp_dir/trans_window_associations.tsv" \
  > "$tmp_dir/trans_window_associations.tsv.gz"

run_header_case() {
  local case_name="$1"
  local source_header="$2"
  local skip_lines="$3"
  local case_dir="$tmp_dir/$case_name"
  local input_tsv="$case_dir/genome_dosage.tsv"
  local input_bgz="$input_tsv.gz"
  local output_json="$case_dir/outputs.json"
  mkdir -p "$case_dir/run"

  printf '%s\n' "$source_header" > "$input_tsv"
  printf 'chr1\t100\tA\tG\t0\t1\n' >> "$input_tsv"
  printf 'chr1\t150\tC\tT\t1\t2\n' >> "$input_tsv"
  bgzip -c "$input_tsv" > "$input_bgz"
  tabix -f -s 1 -b 2 -e 2 -S "$skip_lines" "$input_bgz"

  printf '[%s] Running genotype-header case %s.\n' \
    "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$case_name" >&2
  miniwdl run \
    -d "$case_dir/run/." \
    -o "$output_json" \
    --verbose \
    --no-color \
    --no-cache \
    --task PrepareWindowGenotypes \
    workflows/prepare_trans_window.wdl \
    window_id=w1 \
    trans_window_associations="$tmp_dir/trans_window_associations.tsv.gz" \
    genome_dosage="$input_bgz" \
    genome_dosage_tbi="$input_bgz.tbi"

  output_dosage="$(
    jq -er '.outputs["PrepareWindowGenotypes.window_dosage"]' "$output_json"
  )"
  expected_header=$'CHROM\tPOS\tREF\tALT\t1001\t1002'
  actual_header="$(head -n 1 "$output_dosage")"
  if [[ "$actual_header" != "$expected_header" ]]; then
    printf 'Header case %s returned:\n%s\nExpected:\n%s\n' \
      "$case_name" "$actual_header" "$expected_header" >&2
    exit 1
  fi
  if [[ "$(wc -l < "$output_dosage")" -ne 3 ]]; then
    printf 'Header case %s did not return one header and two variants.\n' \
      "$case_name" >&2
    exit 1
  fi
}

run_header_case \
  "compressed_first_line" \
  $'CHROM\tPOS\tREF\tALT\t1001\t1002' \
  1
run_header_case \
  "tabix_header" \
  $'#CHROM\tPOS\tREF\tALT\t1001\t1002' \
  0

printf '[%s] Genotype dosage header WDL smoke tests passed.\n' \
  "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" >&2
