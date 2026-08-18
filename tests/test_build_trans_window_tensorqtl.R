#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_path <- if (length(args)) args[[1L]] else "scripts/build_trans_window_tensorqtl.R"
source(script_path)

trans <- tibble::tibble(
  variant_id = c(
    "chr1:100_A_G",
    "chr1:200_A_G",
    "chr1:2500000_A_G",
    "chr1:3000000_A_G"
  ),
  phenotype_id = c("trans_1", "trans_2", "trans_3", "not_significant"),
  modality = c("expression", "expression", "splicing", "splicing"),
  pval = c(1e-10, 2e-10, 3e-10, 1e-4),
  b = c(0.1, 0.2, 0.3, 0.4),
  b_se = rep(0.01, 4L),
  af = rep(0.25, 4L)
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

stopifnot(nrow(result$windows) == 2L)
stopifnot(all(file.exists(result$windows$association_file)))
stopifnot(file.exists(result$windows_path))

expected_columns <- c(
  "variant_id", "phenotype_id", "pval", "b", "b_se", "af",
  "modality", "__index_level_0__"
)

window_files <- result$windows$association_file
window_tables <- purrr::map(window_files, ~ readr::read_tsv(
  .x,
  show_col_types = FALSE
))

stopifnot(all(purrr::map_lgl(window_tables, ~ identical(names(.x), expected_columns))))
stopifnot(all(!purrr::map_lgl(window_tables, ~ "not_significant" %in% .x$phenotype_id)))
stopifnot(sum(purrr::map_int(window_tables, nrow)) == 3L)
stopifnot(all(purrr::map_lgl(
  window_tables,
  ~ identical(as.integer(.x$`__index_level_0__`), seq_len(nrow(.x)) - 1L)
)))

stopifnot(result$windows$n_assoc %>% sum() == 3L)
stopifnot(result$windows$n_variants %>% sum() == 3L)
stopifnot(result$windows$n_phenotypes %>% sum() == 3L)

message("Trans TensorQTL window output tests passed")
