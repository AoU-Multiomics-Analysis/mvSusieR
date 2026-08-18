#!/usr/bin/env Rscript

source("scripts/trans_window_io.R")
source("scripts/trans_window_preprocess.R")
source("scripts/trans_window_cli.R")

args <- parse_cli_args()
windows <- read_windows_manifest(require_cli_arg(args, "windows"))
phenotype_manifest <- read_window_phenotypes_manifest(
  require_cli_arg(args, "window-phenotypes")
)
window_id <- require_cli_arg(args, "window-id")
window <- windows[windows[["window_id"]] == window_id]
if (nrow(window) != 1L) {
  stop("Expected exactly one window row for window_id: ", window_id, call. = FALSE)
}

phenotype_files <- split_cli_paths(require_cli_arg(args, "phenotype-files"))
covariate_files <- split_cli_paths(require_cli_arg(args, "covariate-files"))
dosage <- read_wide_dosage(require_cli_arg(args, "dosage"))
phenotype_data <- read_window_phenotypes(window_id, phenotype_manifest, phenotype_files)
covariates <- read_covariate_matrix(covariate_files)

min_nonzero_fraction <- optional_cli_arg(args, "min-nonzero-fraction")
if (!is.null(min_nonzero_fraction)) min_nonzero_fraction <- as.numeric(min_nonzero_fraction)
prepared <- prepare_window_data(
  window = window,
  phenotype_data = phenotype_data,
  dosage = dosage,
  covariates = covariates,
  keep_samples = optional_cli_arg(args, "keep-samples"),
  min_genotype_variance = as_cli_numeric(args, "min-genotype-variance", 1e-8),
  min_phenotype_variance = as_cli_numeric(args, "min-phenotype-variance", 1e-8),
  min_nonzero_fraction = min_nonzero_fraction
)

save_rds_checked(prepared, require_cli_arg(args, "output"))
