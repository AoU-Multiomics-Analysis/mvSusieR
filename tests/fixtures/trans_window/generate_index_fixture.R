#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: generate_index_fixture.R OUTPUT_DIR", call. = FALSE)
}

output_dir <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fixture <- data.frame(
  contig = c("chr2", "chr1", "chr1", "chr1", "chr2"),
  chromStart = c(400L, 220L, 100L, 120L, 50L),
  chromEnd = c(450L, 260L, 160L, 180L, 90L),
  molecular_trait_id = c("feature_e", "feature_d", "feature_a", "feature_b", "feature_c"),
  sample_1 = c(5.1, 4.1, 1.1, 2.1, 3.1),
  sample_2 = c(5.2, 4.2, 1.2, 2.2, 3.2),
  sample_3 = c(5.3, 4.3, 1.3, 2.3, 3.3),
  check.names = FALSE
)

write.table(
  fixture,
  file = file.path(output_dir, "unsorted_phenotypes.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
