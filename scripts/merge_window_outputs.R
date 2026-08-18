#!/usr/bin/env Rscript

source("scripts/trans_window_cli.R")
suppressPackageStartupMessages(library(data.table))

args <- parse_cli_args()
output_dir <- require_cli_arg(args, "output-dir")
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

merge_gz_table("variant-pips", "variant_pip.tsv.gz")
merge_gz_table("credible-sets", "credible_sets.tsv.gz")
merge_gz_table("component-effects", "component_effects.tsv.gz")

qc_tables <- read_table_list("window-qc")
fwrite(
  rbindlist(qc_tables, fill = TRUE, use.names = TRUE),
  file.path(output_dir, "window_qc.tsv"),
  sep = "\t"
)
