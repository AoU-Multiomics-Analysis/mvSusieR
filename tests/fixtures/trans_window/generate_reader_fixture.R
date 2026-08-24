#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: generate_reader_fixture.R OUTPUT_DIR", call. = FALSE)
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

canonical_samples <- paste0("X", 1:6)

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
    outcome_key = c(
      "expression::ENSG000001.1",
      "splicing::splice_1",
      "protein::prot_1"
    ),
    phenotype_id = c("ENSG000001.1", "splice_1", "prot_1"),
    modality = c("expression", "splicing", "protein"),
    phenotype_file = c("expression.tsv", "splicing.tsv", "protein.tsv")
  ),
  file.path(output_dir, "window_phenotypes.tsv"), sep = "\t"
)

dosage <- data.table(
  CHROM = c("chr1", "chr1"), POS = c(101L, 202L),
  REF = c("A", "C"), ALT = c("G", "T")
)
dosage[, (canonical_samples) := list(
  c(0, 1), c(1, 0), c(2, 1), c(0, 1), c(1, 2), c(2, 0)
)]
fwrite(dosage, file.path(output_dir, "window_1_dosage.tsv"), sep = "\t")

write_feature_file <- function(path, outcome_key, values, sample_order) {
  aligned <- stats::setNames(values, canonical_samples)
  table <- data.table(
    chr = "chr2", start = 1001L, end = 1100L,
    phenotype_id = outcome_key
  )
  table[, (sample_order) := as.list(aligned[sample_order])]
  fwrite(table, path, sep = "\t")
}

write_feature_file(
  file.path(output_dir, "expression.tsv"),
  "expression::ENSG000001.1",
  c(1, 2, 3, 1.5, 2.5, 3.5),
  c("X3", "X1", "X6", "X2", "X5", "X4")
)
write_feature_file(
  file.path(output_dir, "splicing.tsv"),
  "splicing::splice_1",
  c(0.6, 0.1, 0.4, 0.2, 0.5, 0.3),
  c("X6", "X4", "X2", "X5", "X1", "X3")
)
write_feature_file(
  file.path(output_dir, "protein.tsv"),
  "protein::prot_1",
  c(20, 50, 10, 60, 30, 40),
  c("X2", "X5", "X1", "X4", "X6", "X3")
)

write_covariates <- function(path, pc1, sample_order, covariate_order) {
  aligned_pc1 <- stats::setNames(pc1, canonical_samples)
  aligned_shared <- stats::setNames(c(2, 4, 6, 8, 10, 12), canonical_samples)
  values <- rbind(PC1 = aligned_pc1, SHARED = aligned_shared)
  values <- values[covariate_order, sample_order, drop = FALSE]
  table <- data.table(covariate = rownames(values))
  table[, (sample_order) := as.data.table(values)]
  fwrite(table, path, sep = "\t")
}

write_covariates(
  file.path(output_dir, "expression_covariates.tsv"),
  c(1, 2, 3, 4, 5, 6),
  c("X2", "X1", "X4", "X3", "X6", "X5"),
  c("PC1", "SHARED")
)
write_covariates(
  file.path(output_dir, "splicing_covariates.tsv"),
  c(-1, 1, -1, 1, -1, 1),
  c("X6", "X5", "X4", "X3", "X2", "X1"),
  c("SHARED", "PC1")
)
write_covariates(
  file.path(output_dir, "protein_covariates.tsv"),
  c(0.5, 0.2, -0.1, -0.4, 0.7, 0.9),
  c("X3", "X4", "X1", "X6", "X2", "X5"),
  c("PC1", "SHARED")
)
