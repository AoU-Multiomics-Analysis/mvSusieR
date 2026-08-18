#!/usr/bin/env Rscript

source("scripts/trans_window_model.R")
source("scripts/trans_window_cli.R")

args <- parse_cli_args()
prepared <- readRDS(require_cli_arg(args, "prepared"))
bundle <- readRDS(require_cli_arg(args, "fit"))
output_dir <- require_cli_arg(args, "output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

window_id <- as.character(prepared$qc$window_id)
config <- bundle$metadata$config

pip <- extract_variant_pips(bundle$fit, prepared)
pip[, window_id := window_id]
pip <- merge(
  pip,
  data.table::copy(prepared$variant_metadata),
  by = "variant_id",
  all.x = TRUE,
  sort = FALSE
)
data.table::setcolorder(
  pip,
  c("window_id", "variant_id", "CHROM", "POS", "REF", "ALT", "pip")
)

credible_sets <- extract_credible_sets(bundle$fit, prepared, config)
credible_sets[, window_id := window_id]
data.table::setcolorder(
  credible_sets,
  c("window_id", "component", "variant_id", "alpha", "pip", "coverage", "purity_min", "purity_mean")
)

effects <- extract_component_effects(bundle$fit, prepared)
effects[, window_id := window_id]
data.table::setcolorder(
  effects,
  c("window_id", "component", "variant_id", "phenotype_id", "posterior_mean", "posterior_sd")
)

qc <- data.table::as.data.table(prepared$qc)
qc[, `:=`(
  converged = isTRUE(bundle$metadata$converged),
  niter = bundle$metadata$niter,
  n_credible_sets = length(bundle$fit$sets$cs),
  mvsusieR_version = bundle$metadata$mvsusieR_version,
  prior = bundle$metadata$prior,
  residual_variance_mode = bundle$metadata$residual_variance_mode
)]

data.table::fwrite(
  pip,
  file.path(output_dir, "variant_pip.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)
data.table::fwrite(
  credible_sets,
  file.path(output_dir, "credible_sets.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)
data.table::fwrite(
  effects,
  file.path(output_dir, "component_effects.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)
data.table::fwrite(qc, file.path(output_dir, "window_qc.tsv"), sep = "\t")
