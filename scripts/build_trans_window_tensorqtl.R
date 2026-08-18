#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

tensorqtl_columns <- c(
  "variant_id", "phenotype_id", "pval", "b", "b_se", "af"
)

output_columns <- c(tensorqtl_columns, "modality")

require_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))

  if (length(missing) > 0L) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

format_coordinate <- function(value) {
  format(value, scientific = FALSE, trim = TRUE)
}

parse_variant_locations <- function(associations, window_size_bp) {
  associations %>%
    mutate(
      chrom = str_extract(variant_id, "^[^:]+"),
      pos = as.integer(str_extract(variant_id, "(?<=:)[0-9]+"))
    ) %>%
    filter(!is.na(chrom), !is.na(pos)) %>%
    mutate(
      window_index = floor(pos / window_size_bp),
      start = window_index * window_size_bp,
      end = (window_index + 1) * window_size_bp,
      window_id = paste(
        chrom,
        format_coordinate(start),
        format_coordinate(end),
        sep = "_"
      )
    )
}

build_trans_window_tensorqtl_outputs <- function(
    trans_associations,
    output_dir,
    window_size_bp = 2e6,
    trans_p_threshold = 1e-8
) {
  require_columns(
    trans_associations,
    c(tensorqtl_columns, "modality"),
    "Trans TensorQTL associations"
  )

  if (length(window_size_bp) != 1L ||
      !is.finite(window_size_bp) || window_size_bp <= 0) {
    stop("window_size_bp must be a positive finite number.", call. = FALSE)
  }
  if (length(trans_p_threshold) != 1L ||
      !is.finite(trans_p_threshold) || trans_p_threshold <= 0) {
    stop("trans_p_threshold must be a positive finite number.", call. = FALSE)
  }

  trans <- trans_associations %>%
    arrange(pval) %>%
    distinct(variant_id, phenotype_id, modality, .keep_all = TRUE) %>%
    filter(!is.na(pval), pval < trans_p_threshold) %>%
    parse_variant_locations(window_size_bp)

  if (nrow(trans) == 0L) {
    stop("No trans associations pass the p-value threshold.", call. = FALSE)
  }

  windows <- trans %>%
    distinct(window_id, chrom, start, end) %>%
    arrange(chrom, start)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  association_manifest <- map_dfr(seq_len(nrow(windows)), function(i) {
    current_window <- windows[i, ]
    current_window_id <- current_window$window_id[[1L]]
    safe_id <- str_replace_all(current_window_id, "[^A-Za-z0-9_.-]", "_")
    association_file <- file.path(
      output_dir,
      paste0(safe_id, ".tensorqtl.tsv.gz")
    )

    window_rows <- trans %>%
      filter(window_id == current_window_id) %>%
      select(all_of(output_columns)) %>%
      arrange(variant_id, phenotype_id) %>%
      mutate(`__index_level_0__` = row_number() - 1L)

    write_tsv(window_rows, association_file)

    tibble(
      window_id = current_window_id,
      association_file = normalizePath(association_file, mustWork = TRUE),
      n_assoc = nrow(window_rows),
      n_variants = n_distinct(window_rows$variant_id),
      n_phenotypes = n_distinct(
        window_rows$modality,
        window_rows$phenotype_id
      )
    )
  })

  windows_output <- windows %>%
    left_join(association_manifest, by = "window_id")

  windows_path <- file.path(output_dir, "windows.tsv")
  write_tsv(windows_output, windows_path)

  list(
    windows = windows_output,
    windows_path = normalizePath(windows_path, mustWork = TRUE)
  )
}

read_trans_input_files <- function(paths, labels) {
  paths <- as.character(paths)
  labels <- as.character(labels)

  if (length(paths) == 0L || length(paths) != length(labels)) {
    stop(
      "The number of trans input files must match the number of labels.",
      call. = FALSE
    )
  }
  if (anyNA(paths) || any(!nzchar(paths))) {
    stop("Trans input file paths cannot be empty.", call. = FALSE)
  }
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("Trans input labels cannot be empty.", call. = FALSE)
  }

  map2_dfr(paths, labels, function(path, modality) {
    read_tsv(path, show_col_types = FALSE) %>%
      mutate(modality = modality)
  })
}

main <- function() {
  source("scripts/trans_window_cli.R")
  args <- parse_cli_args(
    option_list = list(
      optparse::make_option("--trans-files", type = "character"),
      optparse::make_option("--trans-labels", type = "character"),
      optparse::make_option("--output-dir", type = "character"),
      optparse::make_option("--window-size-bp", type = "double", default = 2e6),
      optparse::make_option("--trans-p-threshold", type = "double", default = 1e-8)
    ),
    description = "Build trans TensorQTL files by genomic window."
  )

  trans_files <- split_cli_paths(require_cli_arg(args, "trans_files"))
  trans_labels <- split_cli_paths(require_cli_arg(args, "trans_labels"))

  trans_associations <- read_trans_input_files(
    paths = trans_files,
    labels = trans_labels
  )

  result <- build_trans_window_tensorqtl_outputs(
    trans_associations = trans_associations,
    output_dir = require_cli_arg(args, "output_dir"),
    window_size_bp = as.numeric(args$window_size_bp),
    trans_p_threshold = as.numeric(args$trans_p_threshold)
  )

  message("Wrote ", nrow(result$windows), " trans windows.")
  message("Windows: ", result$windows_path)
}

if (sys.nframe() == 0L) {
  main()
}
