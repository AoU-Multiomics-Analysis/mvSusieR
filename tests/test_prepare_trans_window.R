#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

source("scripts/prepare_trans_window.R")

fixture_dir <- tempfile("prepare-trans-window-fixture-")
dir.create(fixture_dir, recursive = TRUE)
on.exit(unlink(fixture_dir, recursive = TRUE), add = TRUE)

system2(
  command = "Rscript",
  args = c(
    "tests/fixtures/trans_window/generate_prepare_fixture.R",
    fixture_dir
  )
)

windows <- read_prepare_window_manifest(file.path(fixture_dir, "windows.tsv"))
trans_associations <- read_tsv(
  file.path(fixture_dir, "w1.tensorqtl.tsv.gz"),
  show_col_types = FALSE
)
phenotype_inputs <- tibble(
  modality = c("expression", "splicing"),
  phenotype_file = file.path(fixture_dir, c("expression.bed.gz", "splicing.bed.gz"))
)

with_cis <- prepare_trans_window_data(
  windows = windows,
  window_id = "w1",
  trans_associations = trans_associations,
  phenotype_inputs = phenotype_inputs,
  output_dir = file.path(fixture_dir, "with_cis"),
  extract_cis_window_phenotypes = TRUE
)

stopifnot(!"window_dosage" %in% names(with_cis))
stopifnot(!"window_manifest" %in% names(with_cis))
stopifnot(file.exists(with_cis$window_phenotypes))
stopifnot(file.exists(with_cis$window_qc))
stopifnot(length(with_cis$phenotype_subsets) == 2L)

manifest <- read_tsv(with_cis$window_phenotypes, show_col_types = FALSE)
stopifnot(
  identical(
    manifest %>% arrange(modality, phenotype_id) %>% pull(phenotype_id),
    c("ENSG_CIS", "ENSG_TRANS", "splice_cis", "splice_trans")
  )
)
stopifnot(all(manifest$window_id == "w1"))
stopifnot(
  all(file.exists(file.path(dirname(with_cis$window_phenotypes), "phenotype_subsets", manifest$phenotype_file)))
)

expression_subset <- read_tsv(
  file.path(
    dirname(with_cis$window_phenotypes),
    "phenotype_subsets",
    unique(manifest$phenotype_file[manifest$modality == "expression"])
  ),
  show_col_types = FALSE
)
stopifnot(identical(expression_subset$phenotype_id, c("ENSG_CIS", "ENSG_TRANS")))

without_cis <- prepare_trans_window_data(
  windows = windows,
  window_id = "w1",
  trans_associations = trans_associations,
  phenotype_inputs = phenotype_inputs,
  output_dir = file.path(fixture_dir, "without_cis"),
  extract_cis_window_phenotypes = FALSE
)

manifest_without_cis <- read_tsv(
  without_cis$window_phenotypes,
  show_col_types = FALSE
)
stopifnot(
  identical(
    manifest_without_cis %>% arrange(modality, phenotype_id) %>% pull(phenotype_id),
    c("ENSG_TRANS", "splice_trans")
  )
)

message("Single-window preparation tests passed")
