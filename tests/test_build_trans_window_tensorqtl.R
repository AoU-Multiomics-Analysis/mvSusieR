#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_path <- if (length(args)) args[[1L]] else "scripts/build_trans_window_tensorqtl.R"
source(script_path)

trans <- tibble::tibble(
  variant_id = c(
    "chr1:100_A_G",
    "chr1:200_A_G",
    "chr1:300_A_G",
    "chr1:2500000_A_G",
    "chr1:3000000_A_G"
  ),
  phenotype_id = c("trans_1", "trans_2", "trans_1", "trans_3", "not_significant"),
  modality = c("expression", "expression", "expression", "splicing", "splicing"),
  pval = c(1e-10, 2e-10, 1e-12, 3e-10, 1e-4),
  b = c(0.1, 0.2, 0.15, 0.3, 0.4),
  b_se = rep(0.01, 5L),
  af = rep(0.25, 5L)
)

input_dir <- tempfile("tensorqtl_inputs_")
dir.create(input_dir)
expression_path <- file.path(input_dir, "expression.tsv")
splicing_path <- file.path(input_dir, "splicing.tsv")
readr::write_tsv(trans[1:2, ], expression_path)
readr::write_tsv(trans[3:4, ], splicing_path)

combined_inputs <- read_trans_input_files(
  paths = c(expression_path, splicing_path),
  labels = c("expression", "splicing")
)
stopifnot(identical(
  sort(unique(combined_inputs$modality)),
  c("expression", "splicing")
))
stopifnot(nrow(combined_inputs) == 4L)

output_dir <- tempfile("tensorqtl_windows_")
result <- build_trans_window_tensorqtl_outputs(
  trans_associations = trans,
  output_dir = output_dir,
  window_size_bp = 2e6,
  trans_p_threshold = 1e-8
)

stopifnot(file.exists(result$associations_path))
stopifnot(basename(result$associations_path) == "trans_window_associations.tsv.gz")
associations <- readr::read_tsv(
  result$associations_path,
  show_col_types = FALSE
)
expected_columns <- c(
  "window_id", "chrom", "start", "end", "modality",
  "molecular_trait_id", "p_value"
)
stopifnot(identical(names(associations), expected_columns))
stopifnot(nrow(associations) == 3L)
stopifnot(!"not_significant" %in% associations$molecular_trait_id)
stopifnot(nrow(dplyr::distinct(associations, window_id, modality, molecular_trait_id)) == 3L)
stopifnot(all(associations$p_value < 1e-8))
stopifnot(identical(sort(unique(associations$window_id)), c("chr1_0_2000000", "chr1_2000000_4000000")))
stopifnot(associations$p_value[associations$molecular_trait_id == "trans_1"] == 1e-12)
stopifnot(isTRUE(all.equal(
  as.data.frame(result$associations),
  as.data.frame(associations),
  check.attributes = FALSE
)))

boundary <- tibble::tibble(
  variant_id = "chr1:2000000_A_G",
  phenotype_id = "boundary_trait",
  modality = "expression",
  pval = 1e-10,
  b = 0.1,
  b_se = 0.01,
  af = 0.25
)
boundary_result <- build_trans_window_tensorqtl_outputs(
  trans_associations = boundary,
  output_dir = tempfile("tensorqtl_boundary_"),
  window_size_bp = 2e6,
  trans_p_threshold = 1e-8
)
stopifnot(boundary_result$associations$window_id == "chr1_0_2000000")

message("Trans TensorQTL window output tests passed")
