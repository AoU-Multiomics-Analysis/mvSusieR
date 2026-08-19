#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: generate_reader_fixture.R OUTPUT_DIR", call. = FALSE)
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_ids <- paste0("X", 1:6)

fwrite(
  data.table(
    window_id = "w1", chrom = "chr1", start = 100L, end = 300L,
    dosage_file = "window_1_dosage.tsv"
  ),
  file.path(output_dir, "windows.tsv"), sep = "\t"
)
fwrite(
  data.table(
    window_id = "w1",
    phenotype_id = c("ENSG000001.1", "splice_1", "tx_1"),
    modality = c("expression", "splicing", "isoform_usage"),
    phenotype_file = c("expression.tsv", "splicing.tsv", "isoform_usage.tsv")
  ),
  file.path(output_dir, "window_phenotypes.tsv"), sep = "\t"
)

dosage <- data.table(
  CHROM = c("chr1", "chr1"), POS = c(101L, 202L),
  REF = c("A", "C"), ALT = c("G", "T")
)
dosage[, (sample_ids) := list(c(0, 1), c(1, 0), c(2, 1), c(0, 1), c(1, 2), c(2, 0))]
fwrite(dosage, file.path(output_dir, "window_1_dosage.tsv"), sep = "\t")

write_feature_file <- function(path, metadata, values) {
  table <- cbind(metadata, as.data.frame(values))
  setnames(table, (ncol(metadata) + 1L):ncol(table), sample_ids)
  fwrite(table, path, sep = "\t")
}

write_feature_file(
  file.path(output_dir, "expression.tsv"),
  data.table(chr = "chr2", start = 1001L, end = 1100L, gene_id = "ENSG000001.1"),
  matrix(c(1, 2, 3, 1.5, 2.5, 3.5), nrow = 1L)
)
write_feature_file(
  file.path(output_dir, "splicing.tsv"),
  data.table(chr = "chr3", start = 2001L, end = 2100L, phenotype_id = "splice_1"),
  matrix(c(0.6, 0.1, 0.4, 0.2, 0.5, 0.3), nrow = 1L)
)
write_feature_file(
  file.path(output_dir, "isoform_usage.tsv"),
  data.table(transcript_id = "tx_1", transcript_name = "TX1"),
  matrix(c(20, 50, 10, 60, 30, 40), nrow = 1L)
)

covariates <- data.table(
  covariate = c("COV1", "COV2"),
  X1 = c(0.2, -0.4), X2 = c(-0.3, 0.2), X3 = c(0.1, 0.5),
  X4 = c(0.4, -0.1), X5 = c(-0.1, 0.3), X6 = c(0.3, -0.2)
)
fwrite(covariates, file.path(output_dir, "covariates.tsv"), sep = "\t")

fwrite(
  data.table(
    covariate = "EXPR_COV",
    X1 = 0.7, X2 = -0.2, X3 = 0.4,
    X4 = -0.5, X5 = 0.1, X6 = 0.3
  ),
  file.path(output_dir, "expression_covariates.tsv"),
  sep = "\t"
)
fwrite(
  data.table(
    covariate = "SPLICE_COV",
    X1 = -0.1, X2 = 0.6, X3 = -0.3,
    X4 = 0.2, X5 = 0.5, X6 = -0.4
  ),
  file.path(output_dir, "splicing_covariates.tsv"),
  sep = "\t"
)
