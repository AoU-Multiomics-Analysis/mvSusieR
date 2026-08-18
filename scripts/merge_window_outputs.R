#!/usr/bin/env Rscript

source("scripts/trans_window_cli.R")
suppressPackageStartupMessages(library(data.table))

args <- parse_cli_args(
  option_list = list(
    optparse::make_option("--variant-pips", type = "character"),
    optparse::make_option("--credible-sets", type = "character"),
    optparse::make_option("--component-effects", type = "character"),
    optparse::make_option("--window-qc", type = "character"),
    optparse::make_option("--output-dir", type = "character")
  ),
  description = "Merge trans-window mvSusie summary tables."
)
output_dir <- require_cli_arg(args, "output_dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_table_list <- function(name) {
  paths <- split_cli_paths(require_cli_arg(args, name))
  if (!length(paths)) stop("No files supplied for --", name, call. = FALSE)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing --", name, " files: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  lapply(paths, fread, check.names = FALSE)
}

merge_gz_table <- function(name, output_name) {
  merged <- rbindlist(read_table_list(name), fill = TRUE, use.names = TRUE)
  fwrite(merged, file.path(output_dir, output_name), sep = "\t", compress = "gzip")
}

merge_gz_table("variant_pips", "variant_pip.tsv.gz")
merge_gz_table("credible_sets", "credible_sets.tsv.gz")
merge_gz_table("component_effects", "component_effects.tsv.gz")

qc_tables <- read_table_list("window_qc")
fwrite(
  rbindlist(qc_tables, fill = TRUE, use.names = TRUE),
  file.path(output_dir, "window_qc.tsv"),
  sep = "\t"
)
