#!/usr/bin/env Rscript

source("scripts/trans_window_model.R")
source("scripts/trans_window_cli.R")

args <- parse_cli_args()
prepared <- readRDS(require_cli_arg(args, "prepared"))
config <- make_model_config(
  L = as_cli_integer(args, "L", 10L),
  max_iter = as_cli_integer(args, "max-iter", 100L),
  tol = as_cli_numeric(args, "tol", 1e-4),
  coverage = as_cli_numeric(args, "coverage", 0.95),
  min_abs_corr = as_cli_numeric(args, "min-abs-corr", 0.5),
  n_thread = as_cli_integer(args, "n-thread", 1L)
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
  qc = prepared$qc
)
save_rds_checked(bundle, require_cli_arg(args, "output"))
