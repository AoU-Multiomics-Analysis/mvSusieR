#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: generate_prepare_fixture.R OUTPUT_DIR", call. = FALSE)
}

output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

make_associations <- function(modality, prefix, count) {
  tibble(
    window_id = "w1",
    chrom = "chr1",
    start = 100L,
    end = 200L,
    modality = modality,
    molecular_trait_id = sprintf("%s_%02d", prefix, seq_len(count)),
    p_value = 10^-(seq_len(count) + 2L)
  )
}

associations <- bind_rows(
  make_associations("expression", "expr", 27L),
  make_associations("splicing", "splice", 27L),
  make_associations("protein", "protein", 17L),
  tibble(
    window_id = "w2", chrom = "chr1", start = 300L, end = 400L,
    modality = "expression", molecular_trait_id = "expr_w2", p_value = 1e-20
  )
)
write_tsv(
  associations,
  file.path(output_dir, "trans_window_associations.tsv.gz")
)

make_phenotypes <- function(prefix, count, samples, targets = character()) {
  ids <- c(sprintf("%s_%02d", prefix, seq_len(count)), targets)
  index <- seq_along(ids)
  metadata <- tibble(
    chr = ifelse(index %% 2L == 0L, "chr2", "chr3"),
    start = 500L + index * 10L,
    end = 505L + index * 10L,
    phenotype_id = ids
  )
  values <- vapply(
    seq_along(samples),
    function(sample_index) index + sample_index / 10,
    numeric(length(index))
  )
  colnames(values) <- samples
  bind_cols(metadata, as_tibble(values, .name_repair = "minimal"))
}

write_tsv(
  make_phenotypes(
    "expr", 27L,
    samples = c("X1001", "1002", "1003"),
    targets = "expr_target"
  ),
  file.path(output_dir, "expression.bed.gz")
)
write_tsv(
  make_phenotypes(
    "splice", 27L,
    samples = c("1003", "X1001", "1004"),
    targets = c("splice_target_1", "splice_target_2")
  ),
  file.path(output_dir, "splicing.bed.gz")
)
write_tsv(
  make_phenotypes(
    "protein", 17L,
    samples = c("1005", "1003", "1001")
  ),
  file.path(output_dir, "protein.bed.gz")
)

write_tsv(
  tribble(
    ~window_id, ~modality, ~phenotype_id,
    "w1", "expression", "expr_target",
    "w1", "splicing", "splice_target_1",
    "w1", "splicing", "splice_target_2"
  ),
  file.path(output_dir, "target_phenotypes.tsv")
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
