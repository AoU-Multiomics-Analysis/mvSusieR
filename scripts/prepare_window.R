#!/usr/bin/env Rscript

source("scripts/trans_window_io.R")
source("scripts/trans_window_preprocess.R")
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
    optparse::make_option("--output", type = "character")
  ),
  description = "Prepare one trans-window mvSusie input bundle."
)
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
dosage <- read_wide_dosage(require_cli_arg(args, "dosage"))
phenotype_data <- read_window_phenotypes(window_id, phenotype_manifest, phenotype_files)
covariates_by_modality <- read_covariate_matrices(
  paths = covariate_files,
  modalities = covariate_modalities
)

min_nonzero_fraction <- optional_cli_arg(args, "min_nonzero_fraction")
if (!is.null(min_nonzero_fraction)) min_nonzero_fraction <- as.numeric(min_nonzero_fraction)
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

save_rds_checked(prepared, require_cli_arg(args, "output"))
