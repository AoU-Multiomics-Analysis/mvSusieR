#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

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

read_prepare_window_manifest <- function(path) {
  windows <- read_tsv(
    path,
    col_types = cols(.default = col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  require_columns(
    windows,
    c("window_id", "chrom", "start", "end"),
    "Window manifest"
  )
  windows <- windows %>%
    transmute(
      window_id,
      chrom,
      start = suppressWarnings(as.integer(start)),
      end = suppressWarnings(as.integer(end))
    )
  if (anyDuplicated(windows$window_id)) {
    stop("Window manifest contains duplicate window_id values.", call. = FALSE)
  }
  if (
    any(is.na(windows$start)) || any(is.na(windows$end)) ||
    any(windows$start < 0L) || any(windows$end <= windows$start)
  ) {
    stop(
      "Window coordinates must be valid 0-based half-open intervals.",
      call. = FALSE
    )
  }
  windows
}

select_prepare_window <- function(windows, window_id) {
  selected <- windows %>% filter(.data$window_id == !!window_id)
  if (nrow(selected) != 1L) {
    stop("Expected exactly one window for window_id: ", window_id, call. = FALSE)
  }
  selected
}

select_top_trans_phenotypes <- function(trans_associations, top_n) {
  require_columns(
    trans_associations,
    c("phenotype_id", "modality", "pval"),
    "Trans associations"
  )
  if (
    length(top_n) != 1L || is.na(top_n) ||
    top_n < 1L || top_n != as.integer(top_n)
  ) {
    stop("top_n must be a positive integer.", call. = FALSE)
  }

  associations <- trans_associations %>%
    mutate(.pval = suppressWarnings(as.numeric(.data$pval)))
  if (anyNA(associations$.pval) || any(!is.finite(associations$.pval))) {
    stop("Trans association p-values must be finite numeric values.", call. = FALSE)
  }

  associations %>%
    group_by(.data$modality, .data$phenotype_id) %>%
    summarise(min_pval = min(.data$.pval), .groups = "drop") %>%
    group_by(.data$modality) %>%
    arrange(.data$min_pval, .data$phenotype_id, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    arrange(.data$min_pval, .data$modality, .data$phenotype_id)
}

read_prepare_phenotype_table <- function(path, modality) {
  phenotype_table <- read_tsv(
    path,
    col_types = cols(.default = col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  if (ncol(phenotype_table) < 4L) {
    stop(
      paste0(
        "Phenotype file must contain chromosome, start, end, phenotype ID, ",
        "and at least one sample column: "
      ),
      path,
      call. = FALSE
    )
  }

  id_column <- names(phenotype_table)[[4L]]
  phenotype_table <- phenotype_table %>%
    mutate(
      .phenotype_id = as.character(.data[[id_column]]),
      .chrom = as.character(.data[[names(phenotype_table)[[1L]]]]),
      .start = suppressWarnings(as.integer(.data[[names(phenotype_table)[[2L]]]])),
      .end = suppressWarnings(as.integer(.data[[names(phenotype_table)[[3L]]]])),
      .modality = modality,
      .row_id = row_number()
    )
  if (
    anyNA(phenotype_table$.start) || anyNA(phenotype_table$.end) ||
    any(phenotype_table$.end <= phenotype_table$.start)
  ) {
    stop("Phenotype coordinates are invalid in: ", path, call. = FALSE)
  }
  if (anyDuplicated(phenotype_table$.phenotype_id)) {
    stop("Phenotype file contains duplicate IDs: ", path, call. = FALSE)
  }
  phenotype_table
}

select_prepare_phenotypes <- function(
    phenotype_table,
    window,
    trans_ids,
    extract_cis_window_phenotypes
) {
  cis_rows <- phenotype_table %>%
    filter(
      .data$.chrom == window$chrom[[1L]],
      .data$.start < window$end[[1L]],
      .data$.end > window$start[[1L]]
    ) %>%
    mutate(.selection = "cis")

  trans_rows <- phenotype_table %>%
    filter(.data$.phenotype_id %in% trans_ids) %>%
    mutate(.selection = "trans")

  selected <- if (isTRUE(extract_cis_window_phenotypes)) {
    bind_rows(trans_rows, cis_rows) %>%
      distinct(.data$.phenotype_id, .keep_all = TRUE) %>%
      arrange(.data$.row_id)
  } else {
    trans_rows %>%
      distinct(.data$.phenotype_id, .keep_all = TRUE) %>%
      arrange(.data$.row_id)
  }

  list(
    table = selected,
    n_input = nrow(phenotype_table),
    n_trans = nrow(trans_rows),
    n_cis = if (isTRUE(extract_cis_window_phenotypes)) nrow(cis_rows) else 0L,
    n_retained = nrow(selected)
  )
}

write_prepare_phenotype_subset <- function(selected, modality, output_dir, window_id) {
  subset_dir <- file.path(output_dir, "phenotype_subsets")
  dir.create(subset_dir, recursive = TRUE, showWarnings = FALSE)
  safe_modality <- str_replace_all(modality, "[^A-Za-z0-9_.-]", "_")
  safe_window_id <- str_replace_all(window_id, "[^A-Za-z0-9_.-]", "_")
  output_path <- file.path(
    subset_dir,
    paste0(safe_window_id, ".", safe_modality, ".bed.gz")
  )

  output_table <- selected %>%
    select(-starts_with("."))
  write_tsv(output_table, output_path)
  normalizePath(output_path, mustWork = TRUE)
}

prepare_trans_window_data <- function(
    windows,
    window_id,
    trans_associations,
    phenotype_inputs,
    output_dir,
    extract_cis_window_phenotypes = TRUE,
    top_n_trans_phenotypes = 25L
) {
  window <- select_prepare_window(windows, window_id)
  require_columns(
    trans_associations,
    c("phenotype_id", "modality"),
    "Trans associations"
  )
  require_columns(
    phenotype_inputs,
    c("modality", "phenotype_file"),
    "Phenotype inputs"
  )
  if (anyDuplicated(phenotype_inputs$modality)) {
    stop("Phenotype inputs must contain one file per modality.", call. = FALSE)
  }
  if (any(!file.exists(phenotype_inputs$phenotype_file))) {
    stop("At least one phenotype input file does not exist.", call. = FALSE)
  }

  unknown_modalities <- setdiff(
    unique(trans_associations$modality),
    phenotype_inputs$modality
  )
  if (length(unknown_modalities) > 0L) {
    stop(
      "Trans associations contain modalities without phenotype files: ",
      paste(unknown_modalities, collapse = ", "),
      call. = FALSE
    )
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  trans_associations <- select_top_trans_phenotypes(
    trans_associations,
    top_n = top_n_trans_phenotypes
  )

  per_modality <- map_dfr(seq_len(nrow(phenotype_inputs)), function(index) {
    input <- phenotype_inputs[index, ]
    phenotype_table <- read_prepare_phenotype_table(
      input$phenotype_file[[1L]],
      input$modality[[1L]]
    )
    trans_ids <- trans_associations %>%
      filter(.data$modality == input$modality[[1L]]) %>%
      pull(.data$phenotype_id)
    selected <- select_prepare_phenotypes(
      phenotype_table = phenotype_table,
      window = window,
      trans_ids = trans_ids,
      extract_cis_window_phenotypes = extract_cis_window_phenotypes
    )
    if (selected$n_retained == 0L) return(tibble())

    subset_path <- write_prepare_phenotype_subset(
      selected$table,
      input$modality[[1L]],
      output_dir,
      window_id
    )
    tibble(
      window_id = window_id,
      phenotype_id = selected$table$.phenotype_id,
      modality = input$modality[[1L]],
      phenotype_file = basename(subset_path),
      n_input = selected$n_input,
      n_trans = selected$n_trans,
      n_trans_selected = length(trans_ids),
      n_cis = selected$n_cis,
      n_retained = selected$n_retained
    )
  })

  if (nrow(per_modality) == 0L) {
    stop("No trans or cis phenotypes were selected for window: ", window_id, call. = FALSE)
  }

  manifest_path <- file.path(output_dir, "window_phenotypes.tsv")
  write_tsv(
    per_modality %>% select(all_of(c("window_id", "phenotype_id", "modality", "phenotype_file"))),
    manifest_path
  )

  qc_path <- file.path(output_dir, "window_qc.tsv")
  write_tsv(
    per_modality %>%
      distinct(
        .data$window_id,
        .data$modality,
        .data$n_input,
        .data$n_trans,
        .data$n_trans_selected,
        .data$n_cis,
        .data$n_retained
      ) %>%
      mutate(
        top_n_trans_phenotypes = top_n_trans_phenotypes,
        extract_cis_window_phenotypes = extract_cis_window_phenotypes
      ),
    qc_path
  )

  list(
    window_id = window_id,
    window_phenotypes = normalizePath(manifest_path, mustWork = TRUE),
    phenotype_subsets = unique(file.path(output_dir, "phenotype_subsets", basename(per_modality$phenotype_file))),
    window_qc = normalizePath(qc_path, mustWork = TRUE)
  )
}

main <- function() {
  source("scripts/trans_window_cli.R")
  args <- parse_cli_args(
    option_list = list(
      optparse::make_option("--windows", type = "character"),
      optparse::make_option("--window-id", type = "character"),
      optparse::make_option("--trans-associations", type = "character"),
      optparse::make_option("--phenotype-files", type = "character"),
      optparse::make_option("--phenotype-modalities", type = "character"),
      optparse::make_option(
        "--top-n-trans-phenotypes",
        type = "integer",
        default = 25L
      ),
      optparse::make_option(
        "--extract-cis-window-phenotypes",
        type = "logical",
        default = TRUE
      ),
      optparse::make_option("--output-dir", type = "character")
    ),
    description = "Prepare one trans-window mvSuSiE input bundle."
  )

  phenotype_files <- split_cli_paths(require_cli_arg(args, "phenotype_files"))
  phenotype_modalities <- split_cli_paths(
    require_cli_arg(args, "phenotype_modalities")
  )
  if (length(phenotype_files) != length(phenotype_modalities)) {
    stop(
      "The number of phenotype files must match the number of phenotype modalities.",
      call. = FALSE
    )
  }

  result <- prepare_trans_window_data(
    windows = read_prepare_window_manifest(require_cli_arg(args, "windows")),
    window_id = require_cli_arg(args, "window_id"),
    trans_associations = read_tsv(
      require_cli_arg(args, "trans_associations"),
      show_col_types = FALSE
    ),
    phenotype_inputs = tibble(
      modality = phenotype_modalities,
      phenotype_file = phenotype_files
    ),
    output_dir = require_cli_arg(args, "output_dir"),
    extract_cis_window_phenotypes = args$extract_cis_window_phenotypes,
    top_n_trans_phenotypes = as_cli_integer(args, "top_n_trans_phenotypes", 25L)
  )

  message("Prepared window ", result$window_id, ".")
  message("Phenotype manifest: ", result$window_phenotypes)
}

if (sys.nframe() == 0L) {
  main()
}
