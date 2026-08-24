#!/usr/bin/env Rscript

source("scripts/trans_window_io.R")
source("scripts/trans_window_preprocess.R")
source("scripts/trans_window_model.R")
source("scripts/trans_window_prior.R")
source("scripts/trans_window_cli.R")

args <- parse_cli_args(
  option_list = list(
    optparse::make_option("--windows", type = "character"),
    optparse::make_option("--window-phenotypes", type = "character"),
    optparse::make_option("--window-id", type = "character"),
    optparse::make_option("--dosage", type = "character"),
    optparse::make_option("--phenotype-files", type = "character"),
    optparse::make_option("--covariate-files", type = "character"),
    optparse::make_option("--covariate-modalities", type = "character", default = "shared"),
    optparse::make_option("--keep-samples", type = "character", default = NULL),
    optparse::make_option("--min-nonzero-fraction", type = "double", default = NULL),
    optparse::make_option("--min-genotype-variance", type = "double", default = 1e-8),
    optparse::make_option("--min-phenotype-variance", type = "double", default = 1e-8),
    optparse::make_option("--L", type = "integer", default = 10L),
    optparse::make_option("--L-greedy", type = "double", default = NULL),
    optparse::make_option(
      "--greedy-lbf-cutoff",
      type = "double",
      default = 0.1
    ),
    optparse::make_option("--max-iter", type = "integer", default = 100L),
    optparse::make_option("--tol", type = "double", default = 1e-4),
    optparse::make_option("--coverage", type = "double", default = 0.95),
    optparse::make_option("--min-abs-corr", type = "double", default = 0.5),
    optparse::make_option("--n-thread", type = "integer", default = 1L),
    optparse::make_option("--prior-method", type = "character", default = "canonical"),
    optparse::make_option("--mashr-n-pca", type = "integer", default = 5L),
    optparse::make_option("--mashr-seed", type = "integer", default = NULL),
    optparse::make_option("--mashr-strong-lfsr", type = "double", default = 0.05),
    optparse::make_option("--mashr-skip-ed", action = "store_true", default = FALSE),
    optparse::make_option(
      "--fix-residual-variance",
      action = "store_true",
      default = FALSE
    ),
    optparse::make_option("--marginal-output", type = "character", default = NULL),
    optparse::make_option("--prepared-output", type = "character"),
    optparse::make_option("--fit-output", type = "character")
  ),
  description = "Prepare and fit mvSusieR for one trans window."
)

pipeline_log("Reading window and phenotype manifests.")
windows <- read_windows_manifest(require_cli_arg(args, "windows"))
phenotype_manifest <- read_window_phenotypes_manifest(
  require_cli_arg(args, "window_phenotypes")
)
window_id <- require_cli_arg(args, "window_id")
window <- windows[windows[["window_id"]] == window_id]
if (nrow(window) != 1L) {
  stop("Expected exactly one window row for window_id: ", window_id, call. = FALSE)
}

phenotype_files <- split_cli_paths(require_cli_arg(args, "phenotype_files"))
covariate_files <- split_cli_paths(require_cli_arg(args, "covariate_files"))
covariate_modalities <- split_cli_paths(
  optional_cli_arg(args, "covariate_modalities", "shared")
)
pipeline_log("Reading genotype data.")
dosage <- read_wide_dosage(require_cli_arg(args, "dosage"))
pipeline_log(sprintf(
  "Genotype data loaded: %d samples and %d variants.",
  nrow(dosage$X), ncol(dosage$X)
))
pipeline_log("Reading selected phenotype data.")
phenotype_data <- read_window_phenotypes(window_id, phenotype_manifest, phenotype_files)
pipeline_log(sprintf(
  "Phenotype data loaded: %d samples and %d outcomes.",
  nrow(phenotype_data$Y), ncol(phenotype_data$Y)
))
pipeline_log("Reading modality-specific covariates.")
covariates_by_modality <- read_covariate_matrices(
  paths = covariate_files,
  modalities = covariate_modalities
)

min_nonzero_fraction <- optional_cli_arg(args, "min_nonzero_fraction")
if (!is.null(min_nonzero_fraction)) min_nonzero_fraction <- as.numeric(min_nonzero_fraction)
mashr_seed <- optional_cli_arg(args, "mashr_seed")
if (!is.null(mashr_seed)) mashr_seed <- as_cli_integer(args, "mashr_seed", 0L)
L_greedy <- optional_cli_arg(args, "L_greedy")
if (!is.null(L_greedy)) L_greedy <- as_cli_numeric(args, "L_greedy", 0)
pipeline_log("Residualizing genotype and phenotype matrices.")
prepared <- prepare_window_data(
  window = window,
  phenotype_data = phenotype_data,
  dosage = dosage,
  covariates_by_modality = covariates_by_modality,
  keep_samples = optional_cli_arg(args, "keep_samples"),
  min_genotype_variance = as_cli_numeric(args, "min_genotype_variance", 1e-8),
  min_phenotype_variance = as_cli_numeric(args, "min_phenotype_variance", 1e-8),
  min_nonzero_fraction = min_nonzero_fraction
)
pipeline_log(sprintf(
  "Preprocessing complete: %d samples, %d variants, and %d outcomes retained.",
  nrow(prepared$X), ncol(prepared$X), ncol(prepared$Y)
))
pipeline_log("Saving prepared window data.")
save_rds_checked(prepared, require_cli_arg(args, "prepared_output"))
pipeline_log("Prepared window data saved.")

config <- make_model_config(
  L = as_cli_integer(args, "L", 10L),
  L_greedy = L_greedy,
  greedy_lbf_cutoff = as_cli_numeric(args, "greedy_lbf_cutoff", 0.1),
  max_iter = as_cli_integer(args, "max_iter", 100L),
  tol = as_cli_numeric(args, "tol", 1e-4),
  coverage = as_cli_numeric(args, "coverage", 0.95),
  min_abs_corr = as_cli_numeric(args, "min_abs_corr", 0.5),
  n_thread = as_cli_integer(args, "n_thread", 1L),
  prior_method = optional_cli_arg(args, "prior_method", "canonical"),
  mashr_n_pca = as_cli_integer(args, "mashr_n_pca", 5L),
  mashr_seed = mashr_seed,
  mashr_strong_lfsr = as_cli_numeric(args, "mashr_strong_lfsr", 0.05),
  mashr_use_ed = !isTRUE(args$mashr_skip_ed),
  estimate_residual_variance = !isTRUE(args$fix_residual_variance),
  marginal_output = optional_cli_arg(args, "marginal_output")
)
result <- fit_window_mvsusie(prepared, config)
bundle <- list(
  fit = result$fit,
  metadata = result$metadata,
  window = prepared$window,
  variant_metadata = prepared$variant_metadata,
  phenotype_metadata = prepared$phenotype_metadata,
  samples = prepared$samples,
  covariate_rank = prepared$covariate_rank,
  phenotype_covariate_rank = prepared$phenotype_covariate_rank,
  qc = prepared$qc
)
pipeline_log("Saving the mvSuSiE fit bundle.")
save_rds_checked(bundle, require_cli_arg(args, "fit_output"))
pipeline_log("The mvSuSiE fit bundle was saved.")
