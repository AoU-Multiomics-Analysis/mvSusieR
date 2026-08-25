#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

prepare_log <- function(message_text) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), message_text)
}

require_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      label, " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

required_joint_modalities <- function() {
  c("expression", "splicing", "protein")
}

validate_joint_modalities <- function(modalities, label) {
  actual <- sort(unique(as.character(modalities)))
  unsupported <- setdiff(actual, required_joint_modalities())
  if (!length(actual)) {
    stop(label, " must contain at least one supported modality.", call. = FALSE)
  }
  if (length(unsupported)) {
    stop(
      label, " contains an unsupported modality: ",
      paste(unsupported, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

normalize_trans_window_associations <- function(data) {
  require_columns(
    data,
    c(
      "window_id", "chrom", "start", "end", "modality",
      "molecular_trait_id", "p_value"
    ),
    "Trans-window associations"
  )
  associations <- data |>
    transmute(
      window_id = as.character(.data$window_id),
      chrom = as.character(.data$chrom),
      start = suppressWarnings(as.integer(.data$start)),
      end = suppressWarnings(as.integer(.data$end)),
      modality = as.character(.data$modality),
      molecular_trait_id = as.character(.data$molecular_trait_id),
      p_value = suppressWarnings(as.numeric(.data$p_value))
    )
  if (
    any(!nzchar(associations$window_id)) ||
    any(!nzchar(associations$chrom)) ||
    any(!nzchar(associations$modality)) ||
    any(!nzchar(associations$molecular_trait_id))
  ) {
    stop("Trans-window association identifiers cannot be empty.", call. = FALSE)
  }
  if (
    anyNA(associations$start) || anyNA(associations$end) ||
    any(associations$start < 0L) || any(associations$end <= associations$start)
  ) {
    stop(
      "Trans-window coordinates must be valid 0-based half-open intervals.",
      call. = FALSE
    )
  }
  if (
    anyNA(associations$p_value) ||
    any(!is.finite(associations$p_value)) ||
    any(associations$p_value < 0 | associations$p_value > 1)
  ) {
    stop("Trans-window p-values must be finite values from zero through one.", call. = FALSE)
  }
  associations
}

read_trans_window_associations <- function(path) {
  associations <- read_tsv(
    path,
    col_types = cols(.default = col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  normalize_trans_window_associations(associations)
}

select_prepare_window <- function(trans_associations, window_id) {
  selected <- trans_associations |>
    filter(.data$window_id == !!window_id) |>
    distinct(.data$window_id, .data$chrom, .data$start, .data$end)
  if (nrow(selected) != 1L) {
    stop(
      "Expected exactly one coordinate interval for window_id: ",
      window_id,
      call. = FALSE
    )
  }
  selected
}

validate_top_n_by_modality <- function(top_n_by_modality) {
  expected <- required_joint_modalities()
  values <- suppressWarnings(as.numeric(top_n_by_modality))
  if (
    !identical(sort(names(top_n_by_modality)), sort(expected)) ||
    length(values) != length(expected) ||
    anyNA(values) || any(!is.finite(values)) ||
    any(values < 1L) || any(values != as.integer(values))
  ) {
    stop(
      "top_n_by_modality must contain one positive integer per required modality.",
      call. = FALSE
    )
  }
  stats::setNames(as.integer(values), names(top_n_by_modality))
}

select_top_trans_phenotypes <- function(
    trans_associations,
    top_n_by_modality
) {
  require_columns(
    trans_associations,
    c("molecular_trait_id", "modality", "p_value"),
    "Trans associations"
  )
  validate_joint_modalities(trans_associations$modality, "Trans associations")
  top_n_by_modality <- validate_top_n_by_modality(top_n_by_modality)

  eligible <- trans_associations |>
    transmute(
      modality = as.character(.data$modality),
      molecular_trait_id = as.character(.data$molecular_trait_id),
      p_value = suppressWarnings(as.numeric(.data$p_value))
    )
  if (anyNA(eligible$p_value) || any(!is.finite(eligible$p_value))) {
    stop("Trans association p-values must be finite numeric values.", call. = FALSE)
  }

  eligible <- eligible |>
    group_by(.data$modality, .data$molecular_trait_id) |>
    summarise(min_pval = min(.data$p_value), .groups = "drop") |>
    arrange(.data$modality, .data$min_pval, .data$molecular_trait_id)

  eligible |>
    group_by(.data$modality) |>
    group_modify(function(.x, .y) {
      slice_head(
        .x,
        n = unname(top_n_by_modality[[.y$modality[[1L]]]])
      )
    }) |>
    ungroup()
}

read_target_phenotypes <- function(path) {
  targets <- read_tsv(
    path,
    col_types = cols(.default = col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  require_columns(
    targets,
    c("window_id", "modality", "phenotype_id"),
    "Target phenotypes"
  )
  targets <- targets |>
    transmute(
      window_id = as.character(.data$window_id),
      modality = as.character(.data$modality),
      phenotype_id = as.character(.data$phenotype_id)
    )
  if (
    any(!nzchar(targets$window_id)) ||
    any(!nzchar(targets$modality)) ||
    any(!nzchar(targets$phenotype_id))
  ) {
    stop("Target phenotype identifiers cannot be empty.", call. = FALSE)
  }
  if (any(!targets$modality %in% c("expression", "splicing"))) {
    stop("Target phenotypes may contain expression and splicing only.", call. = FALSE)
  }
  if (anyDuplicated(targets[c("window_id", "modality", "phenotype_id")])) {
    stop("Target phenotype rows contain duplicates.", call. = FALSE)
  }
  targets
}

read_prepare_phenotype_table <- function(path, modality) {
  phenotype_table <- read_tsv(
    path,
    col_types = cols(.default = col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  if (ncol(phenotype_table) < 5L) {
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
  phenotype_table <- phenotype_table |>
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
  if (any(!nzchar(phenotype_table$.phenotype_id))) {
    stop("Phenotype IDs cannot be empty in: ", path, call. = FALSE)
  }
  if (anyDuplicated(phenotype_table$.phenotype_id)) {
    stop("Phenotype file contains duplicate IDs: ", path, call. = FALSE)
  }
  phenotype_table
}

select_joint_phenotype_rows <- function(
    phenotype_table,
    modality,
    trans_ids,
    target_ids
) {
  requested_ids <- unique(c(trans_ids, target_ids))
  missing_ids <- setdiff(requested_ids, phenotype_table$.phenotype_id)
  if (length(missing_ids) > 0L) {
    target_missing <- intersect(missing_ids, target_ids)
    if (length(target_missing) > 0L) {
      stop(
        "Target phenotype is absent from the ", modality, " file: ",
        paste(target_missing, collapse = ", "),
        call. = FALSE
      )
    }
    stop(
      "Selected trans phenotype is absent from the ", modality, " file: ",
      paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }
  selected <- phenotype_table[match(requested_ids, phenotype_table$.phenotype_id), ]
  selected$.outcome_key <- paste(modality, selected$.phenotype_id, sep = "::")
  selected$.selection <- ifelse(
    selected$.phenotype_id %in% target_ids,
    "target",
    "trans"
  )
  selected
}

normalize_prepare_sample_ids <- function(ids) {
  ids <- trimws(as.character(ids))
  sub("^X(?=[0-9])", "", ids, perl = TRUE)
}

align_prepare_phenotype_samples <- function(selected_tables) {
  if (!length(selected_tables) || is.null(names(selected_tables))) {
    stop("Contributing phenotype tables must be a named list.", call. = FALSE)
  }
  sample_columns <- lapply(selected_tables, function(data) {
    output_columns <- names(data)[!startsWith(names(data), ".")]
    output_columns[-seq_len(4L)]
  })
  normalized <- lapply(sample_columns, normalize_prepare_sample_ids)
  duplicated <- vapply(normalized, anyDuplicated, integer(1L)) > 0L
  if (any(duplicated)) {
    stop(
      "Sample ID normalization creates duplicates within modality: ",
      paste(names(selected_tables)[duplicated], collapse = ", "),
      call. = FALSE
    )
  }
  shared <- Reduce(intersect, normalized)
  shared <- normalized[[1L]][normalized[[1L]] %in% shared]
  if (!length(shared)) {
    stop(
      "Contributing modalities have no shared phenotype samples.",
      call. = FALSE
    )
  }
  aligned <- Map(function(data, source_columns, sample_ids) {
    metadata_columns <- names(data)[seq_len(4L)]
    internal_columns <- names(data)[startsWith(names(data), ".")]
    selected_columns <- source_columns[match(shared, sample_ids)]
    output <- data |>
      select(all_of(c(metadata_columns, selected_columns, internal_columns)))
    names(output)[seq.int(5L, 4L + length(shared))] <- shared
    output
  }, selected_tables, sample_columns, normalized)
  names(aligned) <- names(selected_tables)
  qc <- tibble(
    modality = names(selected_tables),
    n_input_samples = lengths(normalized),
    n_shared_samples = length(shared),
    n_samples_removed = lengths(normalized) - length(shared)
  )
  list(tables = aligned, qc = qc, shared_samples = shared)
}

write_prepare_phenotype_subset <- function(selected_tables, output_dir) {
  output_path <- file.path(output_dir, "window_phenotypes.bed.gz")
  output_tables <- map(selected_tables, function(selected) {
    output <- selected |> select(-starts_with("."))
    output[[4L]] <- selected$.outcome_key
    output
  })
  if (length(unique(map(output_tables, names))) != 1L) {
    stop(
      "Phenotype files must have identical metadata and sample columns.",
      call. = FALSE
    )
  }
  output_table <- bind_rows(output_tables)
  names(output_table)[seq_len(4L)] <- c("chrom", "start", "end", "phenotype_id")
  write_tsv(output_table, output_path)
  normalizePath(output_path, mustWork = TRUE)
}

prepare_trans_window_data <- function(
    window_id,
    trans_associations,
    expression_phenotypes,
    splicing_phenotypes,
    protein_phenotypes,
    target_phenotypes,
    output_dir,
    top_n_expression = 25L,
    top_n_splicing = 25L,
    top_n_protein = 15L
) {
  prepare_log(paste0("Starting joint phenotype preparation for window ", window_id, "."))
  trans_associations <- normalize_trans_window_associations(trans_associations)
  window_associations <- trans_associations |>
    filter(.data$window_id == !!window_id)
  if (!nrow(window_associations)) {
    stop("No trans associations found for window: ", window_id, call. = FALSE)
  }
  validate_joint_modalities(
    window_associations$modality,
    paste0("Trans associations for window ", window_id)
  )
  select_prepare_window(window_associations, window_id)

  phenotype_inputs <- c(
    expression = expression_phenotypes,
    splicing = splicing_phenotypes,
    protein = protein_phenotypes
  )
  if (any(!file.exists(phenotype_inputs))) {
    stop("Every joint phenotype input file must exist.", call. = FALSE)
  }
  if (!file.exists(target_phenotypes)) {
    stop("The target phenotype input file does not exist.", call. = FALSE)
  }
  top_n_by_modality <- validate_top_n_by_modality(c(
    expression = top_n_expression,
    splicing = top_n_splicing,
    protein = top_n_protein
  ))

  prepare_log(sprintf(
    "Selecting top trans outcomes: expression=%d, splicing=%d, protein=%d.",
    top_n_by_modality[["expression"]],
    top_n_by_modality[["splicing"]],
    top_n_by_modality[["protein"]]
  ))
  selected_trans <- select_top_trans_phenotypes(
    window_associations,
    top_n_by_modality
  )
  for (modality in required_joint_modalities()) {
    prepare_log(sprintf(
      "%s trans selection retained %d of %d available outcomes.",
      modality,
      sum(selected_trans$modality == modality),
      dplyr::n_distinct(
        window_associations$molecular_trait_id[
          window_associations$modality == modality
        ]
      )
    ))
  }
  targets <- read_target_phenotypes(target_phenotypes) |>
    filter(.data$window_id == !!window_id)
  if (!nrow(targets)) {
    stop("No target phenotypes were provided for window: ", window_id, call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  selected_tables <- map(required_joint_modalities(), function(modality) {
    prepare_log(paste0("Reading and selecting ", modality, " phenotypes."))
    phenotype_table <- read_prepare_phenotype_table(
      phenotype_inputs[[modality]],
      modality
    )
    trans_ids <- selected_trans |>
      filter(.data$modality == !!modality) |>
      pull(.data$molecular_trait_id)
    target_ids <- targets |>
      filter(.data$modality == !!modality) |>
      pull(.data$phenotype_id)
    select_joint_phenotype_rows(
      phenotype_table,
      modality,
      trans_ids,
      target_ids
    )
  })
  names(selected_tables) <- required_joint_modalities()
  contributing_tables <- keep(selected_tables, ~ nrow(.x) > 0L)
  if (!length(contributing_tables)) {
    stop("No outcomes were selected for window: ", window_id, call. = FALSE)
  }
  aligned <- align_prepare_phenotype_samples(contributing_tables)
  for (modality in names(aligned$tables)) {
    sample_qc <- aligned$qc |> filter(.data$modality == !!modality)
    prepare_log(sprintf(
      "%s phenotype samples: input=%d, shared=%d, removed=%d.",
      modality,
      sample_qc$n_input_samples,
      sample_qc$n_shared_samples,
      sample_qc$n_samples_removed
    ))
  }

  phenotype_data_path <- write_prepare_phenotype_subset(
    aligned$tables,
    output_dir
  )
  manifest <- imap_dfr(aligned$tables, function(selected, modality) {
    tibble(
      window_id = window_id,
      outcome_key = selected$.outcome_key,
      phenotype_id = selected$.phenotype_id,
      modality = modality,
      phenotype_file = basename(phenotype_data_path)
    )
  })
  if (anyDuplicated(manifest$outcome_key)) {
    stop("Joint phenotype outcome keys must be unique.", call. = FALSE)
  }

  manifest_path <- file.path(output_dir, "window_phenotypes.tsv")
  write_tsv(manifest, manifest_path)

  qc <- imap_dfr(selected_tables, function(selected, modality) {
    tibble(
      window_id = window_id,
      modality = modality,
      n_input = nrow(read_prepare_phenotype_table(
        phenotype_inputs[[modality]],
        modality
      )),
      n_trans_eligible = dplyr::n_distinct(
        window_associations$molecular_trait_id[
          window_associations$modality == modality
        ]
      ),
      n_trans_selected = sum(selected_trans$modality == modality),
      n_targets = sum(targets$modality == modality),
      n_retained = nrow(selected),
      top_n = top_n_by_modality[[modality]]
    )
  }) |>
    left_join(aligned$qc, by = "modality")
  qc_path <- file.path(output_dir, "window_qc.tsv")
  write_tsv(qc, qc_path)

  required_outputs <- c(phenotype_data_path, manifest_path, qc_path)
  if (any(!file.exists(required_outputs)) || any(file.info(required_outputs)$size == 0)) {
    stop("Joint phenotype preparation did not write every required output.", call. = FALSE)
  }
  prepare_log(sprintf(
    "Joint phenotype preparation complete: %d outcomes retained.",
    nrow(manifest)
  ))

  list(
    window_id = window_id,
    phenotype_data = phenotype_data_path,
    window_phenotypes = normalizePath(manifest_path, mustWork = TRUE),
    window_qc = normalizePath(qc_path, mustWork = TRUE)
  )
}

main <- function() {
  file_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (length(file_argument) != 1L) {
    stop("Cannot determine the prepare_trans_window.R location.", call. = FALSE)
  }
  script_path <- normalizePath(
    gsub(
      "~+~",
      " ",
      sub("^--file=", "", file_argument[[1L]]),
      fixed = TRUE
    ),
    mustWork = TRUE
  )
  source(file.path(dirname(script_path), "trans_window_cli.R"))
  prepare_log(paste0("Loaded command-line helpers from ", dirname(script_path), "."))
  args <- parse_cli_args(
    option_list = list(
      optparse::make_option("--window-id", type = "character"),
      optparse::make_option("--trans-associations", type = "character"),
      optparse::make_option("--expression-phenotypes", type = "character"),
      optparse::make_option("--splicing-phenotypes", type = "character"),
      optparse::make_option("--protein-phenotypes", type = "character"),
      optparse::make_option("--target-phenotypes", type = "character"),
      optparse::make_option("--top-n-expression", type = "integer", default = 25L),
      optparse::make_option("--top-n-splicing", type = "integer", default = 25L),
      optparse::make_option("--top-n-protein", type = "integer", default = 15L),
      optparse::make_option("--output-dir", type = "character")
    ),
    description = "Prepare one joint expression, splicing, and protein window."
  )

  result <- prepare_trans_window_data(
    window_id = require_cli_arg(args, "window_id"),
    trans_associations = read_trans_window_associations(
      require_cli_arg(args, "trans_associations")
    ),
    expression_phenotypes = require_cli_arg(args, "expression_phenotypes"),
    splicing_phenotypes = require_cli_arg(args, "splicing_phenotypes"),
    protein_phenotypes = require_cli_arg(args, "protein_phenotypes"),
    target_phenotypes = require_cli_arg(args, "target_phenotypes"),
    output_dir = require_cli_arg(args, "output_dir"),
    top_n_expression = as_cli_integer(args, "top_n_expression", 25L),
    top_n_splicing = as_cli_integer(args, "top_n_splicing", 25L),
    top_n_protein = as_cli_integer(args, "top_n_protein", 15L)
  )

  prepare_log(paste0("Phenotype manifest saved: ", result$window_phenotypes))
}

if (sys.nframe() == 0L) {
  main()
}
