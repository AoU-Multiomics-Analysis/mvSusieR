#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

image="${1:-trans-window-mvsusie-ci:latest}"
tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/trans-window-wdl.XXXXXX")"
trap 'rm -rf "$tmp_dir"' EXIT
input_dir="$tmp_dir/input"
mkdir -p "$input_dir" "$tmp_dir/raw" "$tmp_dir/prepared"

printf '[%s] Generating the WDL smoke-test inputs.\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" >&2
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$repo_root:/workspace" \
  -v "$input_dir:/inputs" \
  -w /workspace \
  "$image" \
  Rscript tests/fixtures/trans_window/generate_model_fixture.R /inputs

printf '[%s] Running the raw-input WDL path.\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" >&2
miniwdl run workflows/trans_window_mvsusie.wdl \
  -d "$tmp_dir/raw/." \
  -o "$tmp_dir/raw_outputs.json" \
  --no-color \
  --no-cache \
  TransWindowMvSusie.window_id=w1 \
  TransWindowMvSusie.window_manifest="$input_dir/windows.tsv" \
  TransWindowMvSusie.window_phenotypes_tsv="$input_dir/window_phenotypes.tsv" \
  TransWindowMvSusie.dosage="$input_dir/model_dosage.tsv" \
  TransWindowMvSusie.phenotype_data="$input_dir/model_phenotypes.tsv" \
  TransWindowMvSusie.expression_covariates="$input_dir/model_expression_covariates.tsv" \
  TransWindowMvSusie.splicing_covariates="$input_dir/model_splicing_covariates.tsv" \
  TransWindowMvSusie.protein_covariates="$input_dir/model_protein_covariates.tsv" \
  TransWindowMvSusie.start_L=10 \
  TransWindowMvSusie.max_L=10 \
  TransWindowMvSusie.greedy_lbf_cutoff=1000000 \
  TransWindowMvSusie.mashr_seed=1 \
  TransWindowMvSusie.docker_image="$image"

prepared_window="$(jq -r '.["TransWindowMvSusie.prepared_window_output"]' "$tmp_dir/raw_outputs.json")"
test -s "$prepared_window"

printf '[%s] Running the prepared-window WDL path.\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" >&2
miniwdl run workflows/trans_window_mvsusie.wdl \
  -d "$tmp_dir/prepared/." \
  -o "$tmp_dir/prepared_outputs.json" \
  --no-color \
  --no-cache \
  TransWindowMvSusie.window_id=w1 \
  TransWindowMvSusie.prepared_window="$prepared_window" \
  TransWindowMvSusie.start_L=10 \
  TransWindowMvSusie.max_L=10 \
  TransWindowMvSusie.greedy_lbf_cutoff=1000000 \
  TransWindowMvSusie.mashr_seed=1 \
  TransWindowMvSusie.docker_image="$image"

test -s "$(jq -r '.["TransWindowMvSusie.mvsusie_fit"]' "$tmp_dir/prepared_outputs.json")"
if find "$tmp_dir/prepared" -type d -name '*PrepareMvSusieInput*' -print -quit | grep -q .; then
  echo "The prepared-window path ran the preparation task." >&2
  exit 1
fi

printf '[%s] Both single-window WDL paths passed.\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" >&2
