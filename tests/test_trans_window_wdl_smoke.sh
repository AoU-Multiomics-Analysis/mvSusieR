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
miniwdl run \
  -d "$tmp_dir/raw/." \
  -o "$tmp_dir/raw_outputs.json" \
  --verbose \
  --no-color \
  --no-cache \
  workflows/trans_window_mvsusie.wdl \
  window_id=w1 \
  window_manifest="$input_dir/windows.tsv" \
  window_phenotypes_tsv="$input_dir/window_phenotypes.tsv" \
  dosage="$input_dir/model_dosage.tsv" \
  phenotype_data="$input_dir/model_phenotypes.tsv" \
  expression_covariates="$input_dir/model_expression_covariates.tsv" \
  splicing_covariates="$input_dir/model_splicing_covariates.tsv" \
  protein_covariates="$input_dir/model_protein_covariates.tsv" \
  start_L=10 \
  max_L=10 \
  greedy_lbf_cutoff=1000000 \
  mashr_seed=1 \
  docker_image="$image"

prepared_window="$(
  jq -er '.outputs["TransWindowMvSusie.prepared_window_output"]' \
    "$tmp_dir/raw_outputs.json"
)"
assert_output_basenames() {
  local outputs_json="$1"
  while IFS=$'\t' read -r output_name expected_basename; do
    actual_path="$(jq -er ".outputs[\"TransWindowMvSusie.${output_name}\"]" "$outputs_json")"
    if [[ "$(basename "$actual_path")" != "$expected_basename" ]]; then
      printf 'Output %s has basename %s; expected %s.\n' \
        "$output_name" "$(basename "$actual_path")" "$expected_basename" >&2
      exit 1
    fi
  done <<'OUTPUT_NAMES'
prepared_window_output	w1.prepared_window.rds
mvsusie_fit	w1.mvsusie_fit.rds
mashr_training	w1.mashr_training.rds
greedy_L_history	w1.greedy_L_history.tsv
covariate_provenance	w1.covariate_provenance.tsv.gz
run_stdout	w1.run.stdout.log
run_stderr	w1.run.stderr.log
session_info	w1.session_info.txt
variant_pip	w1.variant_pip.tsv.gz
credible_sets	w1.credible_sets.tsv.gz
credible_set_members	w1.credible_set_members.tsv.gz
component_feature_support	w1.component_feature_support.tsv.gz
window_qc	w1.window_qc.tsv
effect_plot_png	w1.effect_plot.png
effect_plot_pdf	w1.effect_plot.pdf
effect_plot_rds	w1.effect_plot.rds
OUTPUT_NAMES
}
assert_output_basenames "$tmp_dir/raw_outputs.json"
printf '[%s] Raw path prepared window: %s\n' \
  "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" \
  "$prepared_window" >&2
if [[ ! -s "$prepared_window" ]]; then
  echo "The raw-input WDL path did not return a prepared window." >&2
  exit 1
fi

printf '[%s] Running the prepared-window WDL path.\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" >&2
miniwdl run \
  -d "$tmp_dir/prepared/." \
  -o "$tmp_dir/prepared_outputs.json" \
  --verbose \
  --no-color \
  --no-cache \
  workflows/trans_window_mvsusie.wdl \
  window_id=w1 \
  prepared_window="$prepared_window" \
  start_L=10 \
  max_L=10 \
  greedy_lbf_cutoff=1000000 \
  mashr_seed=1 \
  docker_image="$image"

prepared_fit="$(
  jq -er '.outputs["TransWindowMvSusie.mvsusie_fit"]' \
    "$tmp_dir/prepared_outputs.json"
)"
if [[ ! -s "$prepared_fit" ]]; then
  echo "The prepared-window WDL path did not return an mvSuSiE fit." >&2
  exit 1
fi
if find "$tmp_dir/prepared" -type d -name '*PrepareMvSusieInput*' -print -quit | grep -q .; then
  echo "The prepared-window path ran the preparation task." >&2
  exit 1
fi
assert_output_basenames "$tmp_dir/prepared_outputs.json"

printf '[%s] Both single-window WDL paths passed.\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" >&2
