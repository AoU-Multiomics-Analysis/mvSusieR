args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: generate_model_fixture.R OUTPUT_DIR", call. = FALSE)
output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(2001)
n_samples <- 50L
n_variants <- 6L
sample_ids <- paste0("X", seq_len(n_samples))
variant_metadata <- data.table::data.table(
  CHROM = rep("chr7", n_variants),
  POS = seq(50300000L, by = 1000L, length.out = n_variants),
  REF = rep(c("A", "C", "G"), length.out = n_variants),
  ALT = rep(c("G", "T", "A"), length.out = n_variants)
)
dosage <- matrix(sample(0:2, n_samples * n_variants, replace = TRUE), nrow = n_samples)
dosage[, 1L] <- rep(c(0, 1, 2, 1, 0), length.out = n_samples)
dosage_table <- cbind(variant_metadata, as.data.frame(t(dosage)))
names(dosage_table)[-(1:4)] <- sample_ids
data.table::fwrite(
  dosage_table,
  file.path(output_dir, "model_dosage.tsv"),
  sep = "\t",
  quote = FALSE
)

write_feature_file <- function(path, metadata, values, sample_ids) {
  table <- cbind(metadata, as.data.frame(values))
  names(table)[-(seq_len(ncol(metadata)))] <- sample_ids
  data.table::fwrite(table, path, sep = "\t", quote = FALSE)
}

write_feature_file(
  file.path(output_dir, "model_expression.tsv"),
  data.table::data.table(
    chr = "chr1", start = 1001L, end = 1100L, gene_id = "ENSG_MODEL_EXPR.1"
  ),
  matrix(rnorm(n_samples), nrow = 1L),
  sample_ids
)
write_feature_file(
  file.path(output_dir, "model_splicing.tsv"),
  data.table::data.table(
    chr = "chr2", start = 2001L, end = 2100L, phenotype_id = "splice_model"
  ),
  matrix(rnorm(n_samples), nrow = 1L),
  sample_ids
)
write_feature_file(
  file.path(output_dir, "model_isoform.tsv"),
  data.table::data.table(
    transcript_id = "tx_model", transcript_name = "TX_MODEL"
  ),
  matrix(rnorm(n_samples), nrow = 1L),
  sample_ids
)

covariates <- matrix(rnorm(n_samples * 3L), nrow = 3L)
covariate_table <- rbind(
  c("COV1", covariates[1L, ]),
  c("COV2", covariates[2L, ]),
  c("COV3", covariates[3L, ])
)
colnames(covariate_table) <- c("covariate", sample_ids)
data.table::fwrite(
  as.data.frame(covariate_table),
  file.path(output_dir, "model_covariates.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE
)

data.table::fwrite(
  data.table::data.table(
    window_id = "w1",
    chrom = "chr7",
    start = 50299999L,
    end = 50310000L,
    dosage_file = "model_dosage.tsv"
  ),
  file.path(output_dir, "windows.tsv"),
  sep = "\t",
  quote = FALSE
)
data.table::fwrite(
  data.table::data.table(
    window_id = "w1",
    phenotype_id = c("ENSG_MODEL_EXPR.1", "splice_model", "tx_model"),
    modality = c("expression", "splicing", "isoform_usage"),
    phenotype_file = c("model_expression.tsv", "model_splicing.tsv", "model_isoform.tsv")
  ),
  file.path(output_dir, "window_phenotypes.tsv"),
  sep = "\t",
  quote = FALSE
)
