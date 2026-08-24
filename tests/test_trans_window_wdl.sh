#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

workflow="workflows/trans_window_mvsusie.wdl"
miniwdl check "$workflow"
miniwdl check workflows/prepare_trans_window.wdl

for input in \
  'File phenotype_data' \
  'File expression_covariates' \
  'File splicing_covariates' \
  'File protein_covariates' \
  'Int start_L = 10' \
  'Int step_L = 5' \
  'Int max_L = 40' \
  'Float greedy_lbf_cutoff = 1.0' \
  'Int mashr_n_pca = 5' \
  'Float mashr_strong_lfsr = 0.05'; do
  rg -Fq "$input" "$workflow"
done

for argument in \
  '--expression-covariates' \
  '--splicing-covariates' \
  '--protein-covariates' \
  '--start-L' \
  '--step-L' \
  '--max-L' \
  '--greedy-lbf-cutoff' \
  '--mashr-output' \
  '--greedy-history-output' \
  '--covariate-provenance-output'; do
  rg -Fq -- "$argument" "$workflow"
done

for output in \
  mashr_training \
  greedy_L_history \
  covariate_provenance \
  credible_set_members \
  component_feature_support \
  effect_plot_png \
  effect_plot_pdf \
  effect_plot_rds \
  run_stdout \
  run_stderr \
  session_info; do
  rg -q "$output" "$workflow"
done

test "$(rg -c 'log[(][)]' "$workflow")" -eq 4
test "$(rg -c 'disks: "local-disk 500 SSD"' "$workflow")" -eq 4
test "$(rg -c 'memory: "16 GiB"' "$workflow")" -eq 4
rg -q 'ghcr.io/aou-multiomics-analysis/mvsusier-trans-window-mvsusie:latest' "$workflow"

if rg -q 'canonical|extreme.deconvolution|mashr_use_ed|prior_method|L_greedy|component_effects|covariate_modalities|estimate_residual_variance|isoform' "$workflow"; then
  echo "The joint workflow contains a removed model or input mode." >&2
  exit 1
fi

if rg -q 'Array\[File\] phenotype_files|Array\[File\] covariate_files' "$workflow"; then
  echo "The joint workflow must use one phenotype file and three explicit covariate files." >&2
  exit 1
fi

echo "Trans-window WDL validation passed"
bash tests/test_prepare_trans_window_wdl.sh
bash tests/test_prepare_trans_window_containers.sh
bash tests/test_trans_window_mvsusie_container.sh
