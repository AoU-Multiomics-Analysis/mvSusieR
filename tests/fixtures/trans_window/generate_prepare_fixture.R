#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: generate_prepare_fixture.R OUTPUT_DIR", call. = FALSE)
}

output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv(
  tibble(
    window_id = c("w1", "w2"),
    chrom = c("chr1", "chr1"),
    start = c(100L, 300L),
    end = c(200L, 400L)
  ),
  file.path(output_dir, "windows.tsv")
)

write_tsv(
  tibble(
    variant_id = c("chr2:500_A_G", "chr3:600_C_T", "chr2:501_A_G"),
    phenotype_id = c("ENSG_TRANS", "splice_trans", "ENSG_TRANS_2"),
    pval = c(1e-10, 2e-10, 5e-9),
    b = c(0.4, -0.3, 0.2),
    b_se = c(0.05, 0.04, 0.03),
    af = c(0.2, 0.3, 0.25),
    modality = c("expression", "splicing", "expression")
  ),
  file.path(output_dir, "w1.tensorqtl.tsv.gz")
)

write_tsv(
  tibble(
    chr = c("chr1", "chr1", "chr2", "chr1", "chr2"),
    start = c(150L, 250L, 500L, 350L, 520L),
    end = c(160L, 260L, 510L, 360L, 530L),
    phenotype_id = c("ENSG_CIS", "ENSG_OUT", "ENSG_TRANS", "ENSG_W2", "ENSG_TRANS_2"),
    sample_1 = c(1, 2, 3, 4, 5),
    sample_2 = c(2, 3, 4, 5, 6)
  ),
  file.path(output_dir, "expression.bed.gz")
)

write_tsv(
  tibble(
    chr = c("chr1", "chr3", "chr1"),
    start = c(199L, 600L, 350L),
    end = c(250L, 610L, 360L),
    phenotype_id = c("splice_cis", "splice_trans", "splice_w2"),
    sample_1 = c(0.1, 0.2, 0.3),
    sample_2 = c(0.2, 0.3, 0.4)
  ),
  file.path(output_dir, "splicing.bed.gz")
)

write_tsv(
  tibble(
    CHROM = "chr1",
    POS = 150L,
    REF = "A",
    ALT = "G",
    sample_1 = 0,
    sample_2 = 1
  ),
  file.path(output_dir, "w1.dosage.tsv")
)
