#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

source("scripts/prepare_trans_window.R")

expect_error_matching <- function(expression, pattern) {
  observed <- tryCatch(
    {
      force(expression)
      NA_character_
    },
    error = function(condition) conditionMessage(condition)
  )
  stopifnot(!is.na(observed), grepl(pattern, observed, ignore.case = TRUE))
}

fixture_dir <- tempfile("prepare-trans-window-fixture-")
dir.create(fixture_dir, recursive = TRUE)
on.exit(unlink(fixture_dir, recursive = TRUE), add = TRUE)

status <- system2(
  command = "Rscript",
  args = c("tests/fixtures/trans_window/generate_prepare_fixture.R", fixture_dir)
)
stopifnot(identical(status, 0L))

fixture <- function(name) file.path(fixture_dir, name)
trans_associations <- read_tsv(
  fixture("trans_window_associations.tsv.gz"),
  show_col_types = FALSE
)

# Workflow engines start the container in a task directory. The CLI must load
# its helper relative to the installed script, not relative to that directory.
cli_output_dir <- fixture("cli_from_external_workdir")
external_workdir <- tempfile("prepare-trans-window-workdir-")
dir.create(external_workdir)
old_workdir <- getwd()
setwd(external_workdir)
cli_log <- system2(
  command = "Rscript",
  args = c(
    shQuote(file.path(old_workdir, "scripts", "prepare_trans_window.R")),
    "--window-id", "w1",
    "--trans-associations", shQuote(fixture("trans_window_associations.tsv.gz")),
    "--expression-phenotypes", shQuote(fixture("expression.bed.gz")),
    "--splicing-phenotypes", shQuote(fixture("splicing.bed.gz")),
    "--protein-phenotypes", shQuote(fixture("protein.bed.gz")),
    "--target-phenotypes", shQuote(fixture("target_phenotypes.tsv")),
    "--output-dir", shQuote(cli_output_dir)
  ),
  stdout = TRUE,
  stderr = TRUE
)
setwd(old_workdir)
stopifnot(
  is.null(attr(cli_log, "status")),
  file.exists(file.path(cli_output_dir, "window_phenotypes.tsv")),
  file.exists(file.path(cli_output_dir, "window_phenotypes.bed.gz")),
  file.exists(file.path(cli_output_dir, "window_qc.tsv"))
)

prepare_messages <- capture.output(
  result <- prepare_trans_window_data(
    window_id = "w1",
    trans_associations = trans_associations,
    expression_phenotypes = fixture("expression.bed.gz"),
    splicing_phenotypes = fixture("splicing.bed.gz"),
    protein_phenotypes = fixture("protein.bed.gz"),
    target_phenotypes = fixture("target_phenotypes.tsv"),
    output_dir = fixture("prepared")
  ),
  type = "message"
)

manifest <- read_tsv(result$window_phenotypes, show_col_types = FALSE)
expected_columns <- c(
  "window_id", "outcome_key", "phenotype_id", "modality", "phenotype_file"
)
stopifnot(identical(names(manifest), expected_columns))
expected_counts <- c(expression = 26L, splicing = 27L, protein = 15L)
actual_counts <- table(factor(manifest$modality, levels = names(expected_counts)))
stopifnot(identical(as.integer(actual_counts), unname(expected_counts)))
stopifnot(identical(
  manifest$outcome_key,
  paste(manifest$modality, manifest$phenotype_id, sep = "::")
))
stopifnot(!anyDuplicated(manifest$outcome_key))
stopifnot(all(c(
  "expression::expr_target",
  "splicing::splice_target_1",
  "splicing::splice_target_2"
) %in% manifest$outcome_key))

combined <- read_tsv(result$phenotype_data, show_col_types = FALSE)
stopifnot(identical(combined$phenotype_id, manifest$outcome_key))
stopifnot(nrow(combined) == nrow(manifest))
stopifnot(identical(names(combined)[-(1:4)], c("1001", "1003")))
stopifnot(
  identical(
    as.numeric(combined[combined$phenotype_id == "expression::expr_27", -(1:4)]),
    c(27.1, 27.3)
  ),
  identical(
    as.numeric(combined[combined$phenotype_id == "splicing::splice_27", -(1:4)]),
    c(27.2, 27.1)
  ),
  identical(
    as.numeric(combined[combined$phenotype_id == "protein::protein_17", -(1:4)]),
    c(17.3, 17.2)
  )
)

qc <- read_tsv(result$window_qc, show_col_types = FALSE)
stopifnot(identical(as.integer(qc$top_n), c(25L, 25L, 15L)))
stopifnot(identical(as.integer(qc$n_targets), c(1L, 2L, 0L)))
stopifnot(
  identical(as.integer(qc$n_input_samples), c(3L, 3L, 3L)),
  identical(as.integer(qc$n_shared_samples), c(2L, 2L, 2L)),
  identical(as.integer(qc$n_samples_removed), c(1L, 1L, 1L))
)
stopifnot(all(vapply(c("expression", "splicing", "protein"), function(modality) {
  any(grepl(
    paste0(modality, " phenotype samples: input=3, shared=2, removed=1"),
    prepare_messages,
    fixed = TRUE
  ))
}, logical(1L))))

protein_with_distinct_metadata <- read_tsv(
  fixture("protein.bed.gz"),
  show_col_types = FALSE
)
names(protein_with_distinct_metadata)[seq_len(4L)] <- c(
  "#chrom", "chromStart", "chromEnd", "protein_id"
)
protein_with_distinct_metadata_path <- fixture(
  "protein_distinct_metadata.bed.gz"
)
write_tsv(
  protein_with_distinct_metadata,
  protein_with_distinct_metadata_path
)
distinct_metadata_result <- prepare_trans_window_data(
  "w1", trans_associations,
  fixture("expression.bed.gz"), fixture("splicing.bed.gz"),
  protein_with_distinct_metadata_path, fixture("target_phenotypes.tsv"),
  fixture("distinct_metadata")
)
distinct_metadata_output <- read_tsv(
  distinct_metadata_result$phenotype_data,
  show_col_types = FALSE
)
stopifnot(identical(
  names(distinct_metadata_output)[seq_len(4L)],
  c("chrom", "start", "end", "phenotype_id")
))

disjoint_protein <- read_tsv(fixture("protein.bed.gz"), show_col_types = FALSE)
names(disjoint_protein)[-(1:4)] <- c("2001", "2002", "2003")
disjoint_protein_path <- fixture("disjoint_protein.bed.gz")
write_tsv(disjoint_protein, disjoint_protein_path)
expect_error_matching(
  prepare_trans_window_data(
    "w1", trans_associations,
    fixture("expression.bed.gz"), fixture("splicing.bed.gz"),
    disjoint_protein_path, fixture("target_phenotypes.tsv"),
    fixture("disjoint_samples")
  ),
  "no shared phenotype samples"
)

overlap_targets <- tribble(
  ~window_id, ~modality, ~phenotype_id,
  "w1", "expression", "expr_27",
  "w1", "splicing", "splice_target_1",
  "w1", "splicing", "splice_target_2"
)
overlap_target_path <- fixture("overlap_targets.tsv")
write_tsv(overlap_targets, overlap_target_path)
overlap_result <- prepare_trans_window_data(
  window_id = "w1",
  trans_associations = trans_associations,
  expression_phenotypes = fixture("expression.bed.gz"),
  splicing_phenotypes = fixture("splicing.bed.gz"),
  protein_phenotypes = fixture("protein.bed.gz"),
  target_phenotypes = overlap_target_path,
  output_dir = fixture("overlap")
)
overlap_manifest <- read_tsv(
  overlap_result$window_phenotypes,
  show_col_types = FALSE
)
stopifnot(sum(overlap_manifest$outcome_key == "expression::expr_27") == 1L)

missing_modality <- trans_associations |>
  filter(.data$modality != "protein")
missing_protein_result <- prepare_trans_window_data(
  "w1", missing_modality,
  fixture("expression.bed.gz"), fixture("splicing.bed.gz"),
  fixture("protein.bed.gz"), fixture("target_phenotypes.tsv"),
  fixture("missing_protein")
)
missing_protein_manifest <- read_tsv(
  missing_protein_result$window_phenotypes,
  show_col_types = FALSE
)
stopifnot(!"protein" %in% missing_protein_manifest$modality)
missing_protein_qc <- read_tsv(
  missing_protein_result$window_qc,
  show_col_types = FALSE
)
stopifnot(
  missing_protein_qc$n_trans_selected[missing_protein_qc$modality == "protein"] == 0L
)

single_target_path <- fixture("single_expression_target.tsv")
write_tsv(
  tibble(window_id = "w1", modality = "expression", phenotype_id = "expr_target"),
  single_target_path
)
missing_splicing_result <- prepare_trans_window_data(
  "w1", filter(trans_associations, .data$modality != "splicing"),
  fixture("expression.bed.gz"), fixture("splicing.bed.gz"),
  fixture("protein.bed.gz"), single_target_path,
  fixture("missing_splicing")
)
missing_splicing_manifest <- read_tsv(
  missing_splicing_result$window_phenotypes,
  show_col_types = FALSE
)
stopifnot(!"splicing" %in% missing_splicing_manifest$modality)

single_splicing_target_path <- fixture("single_splicing_target.tsv")
write_tsv(
  tibble(
    window_id = "w1", modality = "splicing",
    phenotype_id = "splice_target_1"
  ),
  single_splicing_target_path
)
missing_expression_result <- prepare_trans_window_data(
  "w1", filter(trans_associations, .data$modality != "expression"),
  fixture("expression.bed.gz"), fixture("splicing.bed.gz"),
  fixture("protein.bed.gz"), single_splicing_target_path,
  fixture("missing_expression")
)
missing_expression_manifest <- read_tsv(
  missing_expression_result$window_phenotypes,
  show_col_types = FALSE
)
stopifnot(!"expression" %in% missing_expression_manifest$modality)

isoform_associations <- trans_associations
isoform_associations$modality[[1L]] <- "isoform_usage"
expect_error_matching(
  prepare_trans_window_data(
    "w1", isoform_associations,
    fixture("expression.bed.gz"), fixture("splicing.bed.gz"),
    fixture("protein.bed.gz"), fixture("target_phenotypes.tsv"),
    fixture("isoform")
  ),
  "unsupported.*modality"
)

missing_target_path <- fixture("missing_target.tsv")
write_tsv(
  tibble(window_id = "w1", modality = "expression", phenotype_id = "absent"),
  missing_target_path
)
expect_error_matching(
  prepare_trans_window_data(
    "w1", trans_associations,
    fixture("expression.bed.gz"), fixture("splicing.bed.gz"),
    fixture("protein.bed.gz"), missing_target_path,
    fixture("missing_target")
  ),
  "target.*absent"
)

duplicate_target_path <- fixture("duplicate_target.tsv")
duplicate_target <- read_tsv(fixture("target_phenotypes.tsv"), show_col_types = FALSE)
write_tsv(bind_rows(duplicate_target, duplicate_target[1L, ]), duplicate_target_path)
expect_error_matching(
  prepare_trans_window_data(
    "w1", trans_associations,
    fixture("expression.bed.gz"), fixture("splicing.bed.gz"),
    fixture("protein.bed.gz"), duplicate_target_path,
    fixture("duplicate_target")
  ),
  "target.*duplicate"
)

message("Joint preparation tests passed")
