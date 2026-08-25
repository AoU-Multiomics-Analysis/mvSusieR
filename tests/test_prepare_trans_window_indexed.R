#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("scripts/prepare_trans_window.R")

indexed_arguments <- c(
  "expression_phenotypes_tbi", "expression_phenotype_lookup",
  "splicing_phenotypes_tbi", "splicing_phenotype_lookup",
  "protein_phenotypes_tbi", "protein_phenotype_lookup"
)
stopifnot(all(indexed_arguments %in% names(formals(prepare_trans_window_data))))

fixture_dir <- tempfile("prepare-trans-window-indexed-")
dir.create(fixture_dir, recursive = TRUE)
on.exit(unlink(fixture_dir, recursive = TRUE), add = TRUE)

status <- system2(
  "Rscript",
  c("tests/fixtures/trans_window/generate_prepare_fixture.R", fixture_dir)
)
stopifnot(identical(status, 0L))

fixture <- function(name) file.path(fixture_dir, name)
modalities <- c("expression", "splicing", "protein")
index_dirs <- setNames(file.path(fixture_dir, paste0(modalities, "_index")), modalities)

for (modality in modalities) {
  status <- system2(
    "bash",
    c(
      "scripts/index_phenotype_bed.sh",
      fixture(paste0(modality, ".bed.gz")),
      modality,
      index_dirs[[modality]],
      "2"
    )
  )
  stopifnot(identical(status, 0L))
}

trans_associations <- read_tsv(
  fixture("trans_window_associations.tsv.gz"),
  show_col_types = FALSE
)

run_prepare <- function(output_name, indexed_modalities) {
  phenotype_path <- function(modality) {
    if (modality %in% indexed_modalities) {
      return(file.path(index_dirs[[modality]], "phenotypes.bed.gz"))
    }
    fixture(paste0(modality, ".bed.gz"))
  }
  index_path <- function(modality) {
    if (!modality %in% indexed_modalities) return(NULL)
    file.path(index_dirs[[modality]], "phenotypes.bed.gz.tbi")
  }
  lookup_path <- function(modality) {
    if (!modality %in% indexed_modalities) return(NULL)
    file.path(index_dirs[[modality]], "phenotype_lookup.tsv.gz")
  }

  prepare_trans_window_data(
    window_id = "w1",
    trans_associations = trans_associations,
    expression_phenotypes = phenotype_path("expression"),
    splicing_phenotypes = phenotype_path("splicing"),
    protein_phenotypes = phenotype_path("protein"),
    target_phenotypes = fixture("target_phenotypes.tsv"),
    output_dir = fixture(output_name),
    expression_phenotypes_tbi = index_path("expression"),
    expression_phenotype_lookup = lookup_path("expression"),
    splicing_phenotypes_tbi = index_path("splicing"),
    splicing_phenotype_lookup = lookup_path("splicing"),
    protein_phenotypes_tbi = index_path("protein"),
    protein_phenotype_lookup = lookup_path("protein")
  )
}

full_result <- run_prepare("full_scan", character())
indexed_messages <- capture.output(
  indexed_result <- run_prepare("indexed", modalities),
  type = "message"
)
mixed_result <- run_prepare("mixed", c("expression", "protein"))

read_result <- function(result) {
  list(
    manifest = read_tsv(result$window_phenotypes, show_col_types = FALSE),
    phenotypes = read_tsv(result$phenotype_data, show_col_types = FALSE),
    qc = read_tsv(result$window_qc, show_col_types = FALSE)
  )
}

full <- read_result(full_result)
indexed <- read_result(indexed_result)
mixed <- read_result(mixed_result)

stopifnot(
  identical(full$manifest, indexed$manifest),
  identical(full$phenotypes, indexed$phenotypes),
  identical(full$manifest, mixed$manifest),
  identical(full$phenotypes, mixed$phenotypes),
  !anyDuplicated(indexed$manifest$outcome_key),
  !anyDuplicated(mixed$manifest$outcome_key)
)

preexisting_qc_columns <- c(
  "window_id", "modality", "n_input", "n_trans_eligible",
  "n_trans_selected", "n_targets", "n_retained", "top_n",
  "n_input_samples", "n_shared_samples", "n_samples_removed"
)
stopifnot(
  identical(
    full$qc[preexisting_qc_columns],
    indexed$qc[preexisting_qc_columns]
  ),
  identical(
    full$qc[preexisting_qc_columns],
    mixed$qc[preexisting_qc_columns]
  ),
  identical(indexed$qc$access_method, rep("tabix", 3L)),
  identical(mixed$qc$access_method, c("tabix", "full_scan", "tabix")),
  all(indexed$qc$query_rows >= indexed$qc$n_retained),
  all(indexed$qc$lookup_seconds >= 0),
  all(indexed$qc$query_seconds >= 0),
  all(indexed$qc$parse_seconds >= 0),
  all(indexed$qc$modality_seconds >= 0)
)

stopifnot(all(vapply(modalities, function(modality) {
  any(grepl(
    paste0(modality, " phenotype access: method=tabix"),
    indexed_messages,
    fixed = TRUE
  ))
}, logical(1L))))

message("Indexed phenotype preparation tests passed")
