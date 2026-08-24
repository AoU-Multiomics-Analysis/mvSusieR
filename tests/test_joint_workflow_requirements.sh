#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

prepare_wdl="workflows/prepare_trans_window.wdl"
model_wdl="workflows/trans_window_mvsusie.wdl"
model_script="scripts/trans_window_model.R"
prior_script="scripts/trans_window_prior.R"
preprocess_script="scripts/trans_window_preprocess.R"
image="envs/trans-window-mvsusie.Dockerfile"
image_ci=".github/workflows/trans-window-mvsusie-image.yml"

for contract in \
  'File expression_phenotypes' \
  'File splicing_phenotypes' \
  'File protein_phenotypes' \
  'File target_phenotypes' \
  'Int top_n_expression = 25' \
  'Int top_n_splicing = 25' \
  'Int top_n_protein = 15'; do
  rg -Fq "$contract" "$prepare_wdl"
done
rg -q 'target_phenotypes' scripts/prepare_trans_window.R

for contract in \
  'File expression_covariates' \
  'File splicing_covariates' \
  'File protein_covariates'; do
  rg -Fq "$contract" "$model_wdl"
done
rg -q 'make_genotype_covariates' "$preprocess_script"
rg -q 'collapse_aligned_covariates' "$preprocess_script"
rg -q 'finite_by_row[(]X_raw[)] & finite_by_row[(]Y_raw[)]' "$preprocess_script"
rg -q 'rank_int' "$preprocess_script"
rg -q 'residualize_matrix' "$preprocess_script"

rg -q 'compute_marginal_bhat_shat_matrix' "$model_script"
rg -q 'cov_pca_fun' "$prior_script"
rg -q 'npc = n_pca' "$prior_script"
rg -q 'mash_model_training_scope = "all_snps_in_window"' "$prior_script"
rg -q 'covariance_input_method = "pca_only"' "$prior_script"
rg -q 'prior[$]xUlist <- lapply' "$prior_script"
rg -q 'U / automatic_scale' "$prior_script"
rg -q 'estimate_residual_variance = FALSE' "$model_script"
rg -q 'estimate_prior_variance = FALSE' "$model_script"
rg -q 'estimate_prior_mixture_weights = FALSE' "$model_script"
rg -q 'verbose = TRUE' "$model_script"

for default in \
  'start_L = 10L' \
  'step_L = 5L' \
  'max_L = 40L' \
  'greedy_lbf_cutoff = 1'; do
  rg -Fq "$default" "$model_script"
done
rg -q 'model_init = previous_fit' "$model_script"
rg -q 'greedy_history = scheduled[$]history' "$model_script"
rg -q 'fit = previous_fit' "$model_script"
rg -q 'mvsusieR::mvsusie_plot' scripts/plot_window_mvsusie.R
rg -q 'conditional_effect = TRUE' scripts/plot_window_mvsusie.R

rg -q '65f3586a865fb6748cb4f9df50510ac577706348' "$image"
rg -q 'ebd1133953005fa70c6b338727b5fe9222e2a1c2' "$image"
test "$(rg -c 'log[(][)]' "$model_wdl")" -eq 4
rg -q 'docker run --rm' "$image_ci"
rg -q 'tests/test_trans_window_r[.]sh' "$image_ci"
rg -q 'tests/test_plot_window_mvsusie[.]R' "$image_ci"

Rscript - <<'RS'
source("scripts/trans_window_io.R")
source("scripts/trans_window_model.R")
stopifnot(identical(required_joint_modalities(), c("expression", "splicing", "protein")))
config <- make_model_config()
stopifnot(
  identical(config$start_L, 10L),
  identical(config$step_L, 5L),
  identical(config$max_L, 40L),
  identical(config$greedy_lbf_cutoff, 1)
)
RS

Rscript - <<'RS'
dockstore <- readLines(".dockstore.yml", warn = FALSE)
paths <- sub("^.*primaryDescriptorPath: /", "", grep(
  "primaryDescriptorPath:", dockstore, value = TRUE
))
stopifnot(length(paths) > 0L, all(file.exists(paths)))
RS

echo "Joint workflow requirements passed"
