#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: generate_model_fixture.R OUTPUT_DIR", call. = FALSE)
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(2001)
n_samples <- 50L
n_variants <- 12L
sample_ids <- paste0("X", seq_len(n_samples))
variant_metadata <- data.table::data.table(
  CHROM = rep("chr7", n_variants),
  POS = seq(50300000L, by = 1000L, length.out = n_variants),
  REF = rep(c("A", "C", "G"), length.out = n_variants),
  ALT = rep(c("G", "T", "A"), length.out = n_variants)
)
dosage <- matrix(
  sample(0:2, n_samples * n_variants, replace = TRUE),
  nrow = n_samples
)
dosage[, 1L] <- rep(c(0, 1, 2, 1, 0), length.out = n_samples)
dosage_table <- cbind(variant_metadata, as.data.frame(t(dosage)))
names(dosage_table)[-(1:4)] <- sample_ids
data.table::fwrite(
  dosage_table,
  file.path(output_dir, "model_dosage.tsv"),
  sep = "\t",
  quote = FALSE
)

write_feature_file <- function(path, modality, phenotype_ids) {
  outcome_keys <- paste(modality, phenotype_ids, sep = "::")
  metadata <- data.table::data.table(
    chr = rep("chr2", length(phenotype_ids)),
    start = seq(1001L, by = 100L, length.out = length(phenotype_ids)),
    end = seq(1100L, by = 100L, length.out = length(phenotype_ids)),
    phenotype_id = outcome_keys
  )
  values <- matrix(rnorm(n_samples * length(phenotype_ids)), nrow = length(phenotype_ids))
  table <- cbind(metadata, as.data.frame(values))
  names(table)[-(seq_len(ncol(metadata)))] <- sample_ids
  data.table::fwrite(table, path, sep = "\t", quote = FALSE)
}

phenotype_ids <- list(
  expression = c("ENSG_MODEL_EXPR.1", "ENSG_MODEL_EXPR.2"),
  splicing = c("splice_model_1", "splice_model_2"),
  protein = c("protein_model_1", "protein_model_2")
)
write_feature_file(
  file.path(output_dir, "model_expression.tsv"),
  "expression",
  phenotype_ids$expression
)
write_feature_file(
  file.path(output_dir, "model_splicing.tsv"),
  "splicing",
  phenotype_ids$splicing
)
write_feature_file(
  file.path(output_dir, "model_protein.tsv"),
  "protein",
  phenotype_ids$protein
)

shared <- rnorm(n_samples)
write_covariates <- function(path, pc1) {
  covariate_table <- rbind(
    c("PC1", pc1),
    c("SHARED", shared)
  )
  colnames(covariate_table) <- c("covariate", sample_ids)
  data.table::fwrite(
    as.data.frame(covariate_table),
    path,
    sep = "\t",
    quote = FALSE,
    col.names = TRUE
  )
}
write_covariates(
  file.path(output_dir, "model_expression_covariates.tsv"),
  rnorm(n_samples)
)
write_covariates(
  file.path(output_dir, "model_splicing_covariates.tsv"),
  rnorm(n_samples)
)
write_covariates(
  file.path(output_dir, "model_protein_covariates.tsv"),
  rnorm(n_samples)
)

data.table::fwrite(
  data.table::data.table(
    window_id = "w1",
    chrom = "chr7",
    start = 50299999L,
    end = 50320000L,
    dosage_file = "model_dosage.tsv"
  ),
  file.path(output_dir, "windows.tsv"),
  sep = "\t",
  quote = FALSE
)

manifest <- data.table::rbindlist(lapply(names(phenotype_ids), function(modality) {
  ids <- phenotype_ids[[modality]]
  data.table::data.table(
    window_id = "w1",
    outcome_key = paste(modality, ids, sep = "::"),
    phenotype_id = ids,
    modality = modality,
    phenotype_file = paste0("model_", modality, ".tsv")
  )
}))
data.table::fwrite(
  manifest,
  file.path(output_dir, "window_phenotypes.tsv"),
  sep = "\t",
  quote = FALSE
)
