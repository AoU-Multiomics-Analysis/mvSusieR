#!/usr/bin/env Rscript

source("scripts/trans_window_model.R")
source("scripts/trans_window_cli.R")

args <- parse_cli_args(
  option_list = list(
    optparse::make_option("--prepared", type = "character"),
    optparse::make_option("--L", type = "integer", default = 10L),
    optparse::make_option("--max-iter", type = "integer", default = 100L),
    optparse::make_option("--tol", type = "double", default = 1e-4),
    optparse::make_option("--coverage", type = "double", default = 0.95),
    optparse::make_option("--min-abs-corr", type = "double", default = 0.5),
    optparse::make_option("--n-thread", type = "integer", default = 1L),
    optparse::make_option("--output", type = "character")
  ),
  description = "Fit mvSusieR for one prepared trans window."
)
prepared <- readRDS(require_cli_arg(args, "prepared"))
config <- make_model_config(
  L = as_cli_integer(args, "L", 10L),
  max_iter = as_cli_integer(args, "max_iter", 100L),
  tol = as_cli_numeric(args, "tol", 1e-4),
  coverage = as_cli_numeric(args, "coverage", 0.95),
  min_abs_corr = as_cli_numeric(args, "min_abs_corr", 0.5),
  n_thread = as_cli_integer(args, "n_thread", 1L)
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
