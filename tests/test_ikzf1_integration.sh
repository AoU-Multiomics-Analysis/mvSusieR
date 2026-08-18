#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
data_dir="${IKZF1_DATA_DIR:-/Users/evinmpadhi/Documents/trans sqtls}"
dosage_file="${IKZF1_DOSAGE_FILE:-}"

if [ -z "$dosage_file" ] || [ ! -s "$dosage_file" ]; then
  echo "Set IKZF1_DOSAGE_FILE to the raw IKZF1 regional dosage TSV/TSV.GZ before running this integration test." >&2
  exit 2
fi

test -s "$data_dir/IKZF1_rsem_transcripts_isopct.tsv"
test -s "$data_dir/IKZF1_transcript_structure_annotations.tsv"
test -s "$data_dir/IKZF1_each_isoform_susie_variant_pip.tsv.gz"
test -s "$data_dir/mvsusie_top25_ikzf1_splicing_bundle.rds"

work_dir="$(mktemp -d "${TMPDIR:-/tmp}/ikzf1-trans-window.XXXXXX")"
trap 'rm -rf "$work_dir"' EXIT

cd "$repo_root"
Rscript tests/test_ikzf1_reference.R \
  --prepare-dir "$work_dir" \
  --dosage "$dosage_file"

Rscript scripts/prepare_window.R \
  --windows "$work_dir/windows.tsv" \
  --window-phenotypes "$work_dir/window_phenotypes.tsv" \
  --window-id ikzf1 \
  --dosage "$work_dir/dosage.tsv" \
  --phenotype-files "$work_dir/phenotypes.tsv" \
  --covariate-files "$work_dir/covariates.tsv" \
  --output "$work_dir/prepared.rds"

Rscript scripts/fit_window.R \
  --prepared "$work_dir/prepared.rds" \
  --output "$work_dir/fit.rds"

mkdir -p "$work_dir/window"
Rscript scripts/summarize_window.R \
  --prepared "$work_dir/prepared.rds" \
  --fit "$work_dir/fit.rds" \
  --output-dir "$work_dir/window"

mkdir -p "$work_dir/merged"
Rscript scripts/merge_window_outputs.R \
  --variant-pips "$work_dir/window/variant_pip.tsv.gz" \
  --credible-sets "$work_dir/window/credible_sets.tsv.gz" \
  --component-effects "$work_dir/window/component_effects.tsv.gz" \
  --window-qc "$work_dir/window/window_qc.tsv" \
  --output-dir "$work_dir/merged"

IKZF1_OUTPUT_DIR="$work_dir" Rscript tests/test_ikzf1_reference.R
echo "IKZF1 direct-R integration test passed"
