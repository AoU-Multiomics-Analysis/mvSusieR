#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

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

Rscript scripts/run_window_mvsusie.R \
  --windows "$input_dir/windows.tsv" \
  --window-phenotypes "$input_dir/window_phenotypes.tsv" \
  --window-id w1 \
  --dosage "$input_dir/model_dosage.tsv" \
  --phenotype-files "$input_dir/model_expression.tsv,$input_dir/model_splicing.tsv,$input_dir/model_isoform.tsv" \
  --covariate-files "$input_dir/model_covariates.tsv" \
  --prior-method mashr \
  --mashr-n-pca 2 \
  --mashr-seed 1 \
  --prepared-output "$tmp_dir/prepared_window.rds" \
  --fit-output "$tmp_dir/mvsusie_fit.rds"

Rscript scripts/summarize_window.R \
  --prepared "$tmp_dir/prepared_window.rds" \
  --fit "$tmp_dir/mvsusie_fit.rds" \
  --output-dir "$tmp_dir/window"

Rscript scripts/merge_window_outputs.R \
  --variant-pips "$tmp_dir/window/variant_pip.tsv.gz" \
  --credible-sets "$tmp_dir/window/credible_sets.tsv.gz" \
  --component-effects "$tmp_dir/window/component_effects.tsv.gz" \
  --window-qc "$tmp_dir/window/window_qc.tsv" \
  --output-dir "$tmp_dir/merged"

for output in \
  "$tmp_dir/prepared_window.rds" \
  "$tmp_dir/mvsusie_fit.rds" \
  "$tmp_dir/window/variant_pip.tsv.gz" \
  "$tmp_dir/window/credible_sets.tsv.gz" \
  "$tmp_dir/window/component_effects.tsv.gz" \
  "$tmp_dir/window/window_qc.tsv" \
  "$tmp_dir/merged/variant_pip.tsv.gz" \
  "$tmp_dir/merged/credible_sets.tsv.gz" \
  "$tmp_dir/merged/component_effects.tsv.gz" \
  "$tmp_dir/merged/window_qc.tsv"; do
  test -s "$output"
done

Rscript - "$tmp_dir/window/window_qc.tsv" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
actual <- data.table::fread(args[[1L]], check.names = FALSE)
expected <- list(
  window_id = "w1", input_samples = 50L, shared_samples = 50L,
  input_variants = 6L, retained_variants = 6L, excluded_variants = 0L,
  input_phenotypes = 3L, retained_phenotypes = 3L,
  excluded_phenotypes = 0L, excluded_samples = 0L, covariate_rank = 4L
)
for (column in names(expected)) stopifnot(identical(actual[[column]][[1L]], expected[[column]]))
RS

Rscript - "$tmp_dir/mvsusie_fit.rds" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
fit <- readRDS(args[[1L]])
stopifnot(identical(fit$metadata$prior, "mashr"))
stopifnot(identical(fit$metadata$covariance_training_scope, "all_snps_in_window"))
stopifnot(isTRUE(fit$metadata$extreme_deconvolution_used))
RS

echo "Task 4 entrypoint tests passed"
