#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

rg -q 'verbose = TRUE' scripts/trans_window_model.R

Rscript tests/test_build_trans_window_tensorqtl.R scripts/build_trans_window_tensorqtl.R

tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/trans-window-r.XXXXXX")"
trap 'rm -rf "$tmp_dir"' EXIT

Rscript tests/test_trans_window_cli.R
Rscript tests/test_prepare_trans_window.R

reader_dir="$tmp_dir/reader"
Rscript tests/fixtures/trans_window/generate_reader_fixture.R "$reader_dir"
Rscript tests/test_trans_window_r.R "$reader_dir"
Rscript tests/test_mashr_prior.R

input_dir="$tmp_dir/input"
Rscript tests/fixtures/trans_window/generate_model_fixture.R "$input_dir"

Rscript scripts/prepare_window.R \
  --windows "$input_dir/windows.tsv" \
  --window-phenotypes "$input_dir/window_phenotypes.tsv" \
  --window-id w1 \
  --dosage "$input_dir/model_dosage.tsv" \
  --phenotype-files "$input_dir/model_phenotypes.tsv" \
  --expression-covariates "$input_dir/model_expression_covariates.tsv" \
  --splicing-covariates "$input_dir/model_splicing_covariates.tsv" \
  --protein-covariates "$input_dir/model_protein_covariates.tsv" \
  --covariate-provenance-output "$tmp_dir/standalone_covariate_provenance.tsv.gz" \
  --output "$tmp_dir/standalone_prepared_window.rds" \
  2>&1 | tee "$tmp_dir/prepare_window.log"

grep -q 'Reading genotype data' "$tmp_dir/prepare_window.log"
grep -q 'Residualizing genotype and phenotype matrices' "$tmp_dir/prepare_window.log"
grep -q 'Prepared window data saved' "$tmp_dir/prepare_window.log"

if Rscript scripts/fit_window.R \
  --prepared "$tmp_dir/standalone_prepared_window.rds" \
  --window-id w1 \
  --start-L 10 \
  --step-L 0 \
  --max-L 40 \
  --mashr-output "$tmp_dir/invalid_mashr.rds" \
  --greedy-history-output "$tmp_dir/invalid_history.tsv" \
  --covariate-provenance-output "$tmp_dir/invalid_provenance.tsv.gz" \
  --output "$tmp_dir/invalid_fit.rds" \
  >"$tmp_dir/invalid_fit.log" 2>&1; then
  echo "fit_window.R accepted a zero greedy L step." >&2
  exit 1
fi
grep -q 'step_L must be a positive integer' "$tmp_dir/invalid_fit.log"

Rscript scripts/fit_window.R \
  --prepared "$tmp_dir/standalone_prepared_window.rds" \
  --window-id w1 \
  --covariate-provenance-output "$tmp_dir/covariate_provenance.tsv.gz" \
  --start-L 10 \
  --step-L 5 \
  --max-L 10 \
  --greedy-lbf-cutoff 1000000 \
  --mashr-n-pca 5 \
  --mashr-seed 1 \
  --mashr-output "$tmp_dir/mashr_training_bundle.rds" \
  --greedy-history-output "$tmp_dir/greedy_L_history.tsv" \
  --output "$tmp_dir/mvsusie_fit.rds" \
  2>&1 | tee "$tmp_dir/run_window.log"

grep -q 'Reading prepared window data' "$tmp_dir/run_window.log"
grep -q 'Computing the all-SNP cross-product' "$tmp_dir/run_window.log"
grep -q 'Starting 5 PCA covariance inputs' "$tmp_dir/run_window.log"
grep -q 'Starting greedy mvSuSiE round 1 at L = 10' "$tmp_dir/run_window.log"
if grep -q 'extreme deconvolution' "$tmp_dir/run_window.log"; then
  echo "The joint model attempted extreme deconvolution." >&2
  exit 1
fi

Rscript scripts/fit_window.R \
  --prepared "$tmp_dir/standalone_prepared_window.rds" \
  --window-id w1 \
  --start-L 10 \
  --step-L 5 \
  --max-L 10 \
  --greedy-lbf-cutoff 1000000 \
  --mashr-n-pca 5 \
  --mashr-seed 1 \
  --mashr-output "$tmp_dir/resumed_mashr_training_bundle.rds" \
  --greedy-history-output "$tmp_dir/resumed_greedy_L_history.tsv" \
  --covariate-provenance-output "$tmp_dir/resumed_covariate_provenance.tsv.gz" \
  --output "$tmp_dir/resumed_mvsusie_fit.rds" \
  2>&1 | tee "$tmp_dir/fit_window.log"

grep -q 'Reading prepared window data' "$tmp_dir/fit_window.log"
grep -q 'Computing the all-SNP cross-product' "$tmp_dir/fit_window.log"
grep -q 'Using fixed mashr weights' "$tmp_dir/fit_window.log"

Rscript scripts/summarize_window.R \
  --prepared "$tmp_dir/standalone_prepared_window.rds" \
  --fit "$tmp_dir/mvsusie_fit.rds" \
  --output-dir "$tmp_dir/window"

Rscript scripts/merge_window_outputs.R \
  --variant-pips "$tmp_dir/window/variant_pip.tsv.gz" \
  --credible-sets "$tmp_dir/window/credible_sets.tsv.gz" \
  --credible-set-members "$tmp_dir/window/credible_set_members.tsv.gz" \
  --component-feature-support "$tmp_dir/window/component_feature_support.tsv.gz" \
  --window-qc "$tmp_dir/window/window_qc.tsv" \
  --output-dir "$tmp_dir/merged"

for output in \
  "$tmp_dir/standalone_prepared_window.rds" \
  "$tmp_dir/standalone_covariate_provenance.tsv.gz" \
  "$tmp_dir/covariate_provenance.tsv.gz" \
  "$tmp_dir/resumed_covariate_provenance.tsv.gz" \
  "$tmp_dir/mvsusie_fit.rds" \
  "$tmp_dir/resumed_mvsusie_fit.rds" \
  "$tmp_dir/mashr_training_bundle.rds" \
  "$tmp_dir/resumed_mashr_training_bundle.rds" \
  "$tmp_dir/greedy_L_history.tsv" \
  "$tmp_dir/resumed_greedy_L_history.tsv" \
  "$tmp_dir/window/variant_pip.tsv.gz" \
  "$tmp_dir/window/credible_sets.tsv.gz" \
  "$tmp_dir/window/credible_set_members.tsv.gz" \
  "$tmp_dir/window/component_feature_support.tsv.gz" \
  "$tmp_dir/window/window_qc.tsv" \
  "$tmp_dir/merged/variant_pip.tsv.gz" \
  "$tmp_dir/merged/credible_sets.tsv.gz" \
  "$tmp_dir/merged/credible_set_members.tsv.gz" \
  "$tmp_dir/merged/component_feature_support.tsv.gz" \
  "$tmp_dir/merged/window_qc.tsv"; do
  test -s "$output"
done

Rscript - "$tmp_dir/standalone_prepared_window.rds" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
prepared <- readRDS(args[[1L]])
stopifnot(inherits(prepared$input_checksums, "data.frame"))
stopifnot(nrow(prepared$input_checksums) == 7L)
stopifnot(all(nzchar(prepared$input_checksums$md5)))
stopifnot(nrow(prepared$covariate_provenance) == 6L)
RS

Rscript - \
  "$tmp_dir/greedy_L_history.tsv" \
  "$tmp_dir/window/credible_sets.tsv.gz" \
  "$tmp_dir/window/credible_set_members.tsv.gz" \
  "$tmp_dir/window/component_feature_support.tsv.gz" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
history <- data.table::fread(args[[1L]], check.names = FALSE)
stopifnot(identical(
  names(history),
  c(
    "round", "requested_L", "fitted_L", "niter", "minimum_lbf",
    "credible_set_count", "supported_component_count", "maximum_alpha", "action"
  )
))
credible_sets <- data.table::fread(args[[2L]], check.names = FALSE)
stopifnot(all(c(
  "window_id", "component", "credible_set_size", "sentinel_variant_id",
  "sentinel_alpha", "coverage", "purity_min", "purity_mean"
) %in% names(credible_sets)))
members <- data.table::fread(args[[3L]], check.names = FALSE)
stopifnot(all(c(
  "window_id", "component", "variant_id", "alpha", "pip", "is_sentinel"
) %in% names(members)))
support <- data.table::fread(args[[4L]], check.names = FALSE)
stopifnot(all(c(
  "window_id", "component", "outcome_key", "modality", "phenotype_id",
  "single_effect_lfsr", "outcome_lbf"
) %in% names(support)))
stopifnot(nrow(support) == 10L * 6L)
stopifnot(!file.exists(file.path(dirname(args[[2L]]), "component_effects.tsv.gz")))
RS

Rscript - "$tmp_dir/window/window_qc.tsv" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
actual <- data.table::fread(args[[1L]], check.names = FALSE)
expected <- list(
  window_id = "w1", input_samples = 50L, shared_samples = 50L,
  input_variants = 12L, retained_variants = 12L, excluded_variants = 0L,
  input_phenotypes = 6L, retained_phenotypes = 6L,
  excluded_phenotypes = 0L, excluded_samples = 0L, covariate_rank = 5L,
  start_L = 10L, step_L = 5L, max_L = 10L,
  greedy_lbf_cutoff = 1e6, L_final = 10L
)
for (column in names(expected)) stopifnot(identical(actual[[column]][[1L]], expected[[column]]))
RS

Rscript - "$tmp_dir/mvsusie_fit.rds" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
fit <- readRDS(args[[1L]])
stopifnot(identical(fit$metadata$prior, "mashr_pca_only"))
stopifnot(identical(fit$metadata$mash_model_training_scope, "all_snps_in_window"))
stopifnot(identical(fit$metadata$covariance_training_scope, "strong_snps_in_window"))
stopifnot(identical(fit$metadata$prior_mixture_weights_mode, "fixed_from_mashr"))
stopifnot(identical(fit$metadata$covariance_input_method, "pca_only"))
stopifnot(nrow(fit$fit$alpha) == 10L)
stopifnot(identical(fit$metadata$L_final, 10L))
RS

Rscript - \
  "$tmp_dir/resumed_mvsusie_fit.rds" \
  "$tmp_dir/standalone_prepared_window.rds" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
fit <- readRDS(args[[1L]])
prepared <- readRDS(args[[2L]])
stopifnot(identical(fit$metadata$prior, "mashr_pca_only"))
stopifnot(identical(fit$metadata$mash_model_training_scope, "all_snps_in_window"))
stopifnot(identical(fit$metadata$covariance_training_scope, "strong_snps_in_window"))
stopifnot(identical(fit$metadata$prior_mixture_weights_mode, "fixed_from_mashr"))
stopifnot(identical(fit$metadata$covariance_input_method, "pca_only"))
stopifnot(identical(fit$metadata$residual_variance_mode, "fixed_initial_covariance"))
stopifnot(nrow(fit$fit$alpha) == 10L)
stopifnot(identical(fit$metadata$config$start_L, 10L))
stopifnot(identical(fit$metadata$config$step_L, 5L))
stopifnot(identical(fit$metadata$config$max_L, 10L))
stopifnot(identical(fit$metadata$config$greedy_lbf_cutoff, 1e6))
stopifnot(isTRUE(all.equal(
  fit$fit$sigma2,
  stats::cov(prepared$Y),
  tolerance = 1e-10
)))
RS

Rscript - "$tmp_dir/mashr_training_bundle.rds" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
mashr_training <- readRDS(args[[1L]])
stopifnot(identical(dim(mashr_training$Bhat), c(12L, 6L)))
stopifnot(identical(dim(mashr_training$Shat), c(12L, 6L)))
stopifnot(identical(mashr_training$pca_requested, 5L))
stopifnot(identical(mashr_training$covariance_input_method, "pca_only"))
RS

echo "Task 4 entrypoint tests passed"
