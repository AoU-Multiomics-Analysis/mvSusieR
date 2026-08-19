suppressPackageStartupMessages(library(data.table))

normalize_sample_ids <- function(ids) {
  ids <- trimws(as.character(ids))
  sub("^X(?=[0-9])", "", ids, perl = TRUE)
}

require_columns <- function(dt, required, label) {
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    stop(
      label, " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

read_windows_manifest <- function(path) {
  dt <- fread(path, check.names = FALSE)
  require_columns(
    dt,
    c("window_id", "chrom", "start", "end", "dosage_file"),
    "Window manifest"
  )
  if (anyDuplicated(dt$window_id)) {
    stop("Window manifest contains duplicate window_id values.", call. = FALSE)
  }
  dt[, start := as.integer(start)]
  dt[, end := as.integer(end)]
  if (anyNA(dt$start) || anyNA(dt$end) || any(dt$start < 0L) || any(dt$end <= dt$start)) {
    stop("Window coordinates must be valid 0-based half-open intervals.", call. = FALSE)
  }
  dt[]
}

read_window_phenotypes_manifest <- function(path) {
  dt <- fread(path, check.names = FALSE)
  require_columns(
    dt,
    c("window_id", "phenotype_id", "modality", "phenotype_file"),
    "Window phenotype manifest"
  )
  if (anyDuplicated(dt[, paste(window_id, phenotype_id, sep = "\r")])) {
    stop(
      "Window phenotype manifest contains duplicate window_id/phenotype_id pairs.",
      call. = FALSE
    )
  }
  allowed <- c("expression", "splicing", "isoform_usage")
  if (any(!dt$modality %in% allowed)) {
    stop(
      "Unsupported phenotype modality. Expected one of: ",
      paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }
  dt[]
}

numeric_matrix <- function(values, label) {
  out <- suppressWarnings(matrix(as.numeric(as.matrix(values)), nrow = nrow(values)))
  if (any(!is.finite(out) & !is.na(out))) {
    stop(label, " contains non-numeric values.", call. = FALSE)
  }
  out
}

read_wide_dosage <- function(path) {
  dt <- fread(path, check.names = FALSE)
  required <- c("CHROM", "POS", "REF", "ALT")
  require_columns(dt, required, "Dosage file")
  sample_columns <- setdiff(names(dt), required)
  if (!length(sample_columns)) {
    stop("Dosage file has no sample columns.", call. = FALSE)
  }
  sample_ids <- normalize_sample_ids(sample_columns)
  if (anyDuplicated(sample_ids)) {
    stop("Dosage file contains duplicate sample IDs.", call. = FALSE)
  }
  X <- numeric_matrix(dt[, ..sample_columns], "Dosage values")
  X <- t(X)
  rownames(X) <- sample_ids
  variant_ids <- paste0(dt$CHROM, ":", dt$POS, "_", dt$REF, "_", dt$ALT)
  if (anyDuplicated(variant_ids)) {
    stop("Dosage file contains duplicate variant IDs.", call. = FALSE)
  }
  colnames(X) <- variant_ids
  list(
    metadata = dt[, ..required],
    X = X,
    sample_ids = sample_ids,
    variant_ids = variant_ids
  )
}

phenotype_layout <- function(modality) {
  if (modality %in% c("expression", "splicing")) {
    return(list(id_column = 4L, metadata_columns = 1:4))
  }
  if (identical(modality, "isoform_usage")) {
    return(list(id_column = 1L, metadata_columns = 1:2))
  }
  stop("Unsupported phenotype modality: ", modality, call. = FALSE)
}

read_phenotype_rows <- function(path, modality, phenotype_ids) {
  dt <- fread(path, check.names = FALSE)
  layout <- phenotype_layout(modality)
  if (ncol(dt) <= max(layout$metadata_columns)) {
    stop("Phenotype file has no sample columns: ", path, call. = FALSE)
  }
  id_values <- as.character(dt[[layout$id_column]])
  if (anyDuplicated(id_values)) {
    stop("Phenotype file contains duplicate feature IDs: ", path, call. = FALSE)
  }
  sample_columns <- setdiff(names(dt), names(dt)[layout$metadata_columns])
  sample_ids <- normalize_sample_ids(sample_columns)
  if (anyDuplicated(sample_ids)) {
    stop("Phenotype file contains duplicate sample IDs: ", path, call. = FALSE)
  }
  phenotype_ids <- as.character(phenotype_ids)
  missing_ids <- setdiff(phenotype_ids, id_values)
  if (length(missing_ids)) {
    stop(
      "Missing phenotype IDs in ", path, ": ",
      paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }
  row_index <- match(phenotype_ids, id_values)
  values <- numeric_matrix(dt[row_index, ..sample_columns], "Phenotype values")
  Y <- t(values)
  rownames(Y) <- sample_ids
  colnames(Y) <- phenotype_ids
  metadata_columns <- names(dt)[layout$metadata_columns]
  list(
    Y = Y,
    metadata = dt[row_index, ..metadata_columns],
    sample_ids = sample_ids,
    phenotype_ids = phenotype_ids,
    modality = modality,
    source_file = path
  )
}

resolve_file_reference <- function(reference, files) {
  files <- as.character(files)
  direct <- files[files == reference]
  if (length(direct) == 1L) return(direct)
  matches <- files[basename(files) == basename(reference)]
  if (length(matches) != 1L) {
    stop(
      "Could not resolve file reference ", reference,
      " among supplied files.",
      call. = FALSE
    )
  }
  matches
}

read_window_phenotypes <- function(window_id, phenotype_manifest, phenotype_files) {
  rows <- phenotype_manifest[phenotype_manifest[["window_id"]] == window_id]
  if (!nrow(rows)) {
    stop("No phenotype rows found for window: ", window_id, call. = FALSE)
  }
  grouped_rows <- split(seq_len(nrow(rows)), rows$phenotype_file)
  parts <- lapply(grouped_rows, function(indices) {
    source_file <- resolve_file_reference(rows$phenotype_file[[indices[[1L]]]], phenotype_files)
    list(
      indices = indices,
      data = read_phenotype_rows(
        source_file,
        rows$modality[[indices[[1L]]]],
        rows$phenotype_id[indices]
      )
    )
  })
  common_samples <- Reduce(intersect, lapply(parts, function(part) part$data$sample_ids))
  if (!length(common_samples)) {
    stop("No shared phenotype samples for window: ", window_id, call. = FALSE)
  }
  Y <- do.call(cbind, lapply(parts, function(part) {
    part$data$Y[common_samples, , drop = FALSE]
  }))
  column_order <- order(unlist(lapply(parts, `[[`, "indices")))
  Y <- Y[, column_order, drop = FALSE]
  metadata <- rbindlist(lapply(parts, function(part) part$data$metadata), fill = TRUE)
  metadata <- metadata[column_order]
  list(
    Y = Y,
    metadata = metadata,
    sample_ids = common_samples,
    phenotype_ids = colnames(Y),
    modalities = rows$modality,
    source_files = rows$phenotype_file
  )
}

read_covariate_matrix <- function(paths) {
  paths <- as.character(paths)
  if (!length(paths)) stop("At least one covariate file is required.", call. = FALSE)
  matrices <- lapply(paths, function(path) {
    dt <- fread(path, check.names = FALSE)
    if (ncol(dt) < 2L) stop("Covariate file has no sample columns: ", path, call. = FALSE)
    covariate_ids <- as.character(dt[[1L]])
    if (anyDuplicated(covariate_ids)) {
      stop("Covariate file contains duplicate covariate IDs: ", path, call. = FALSE)
    }
    sample_columns <- names(dt)[-1L]
    sample_ids <- normalize_sample_ids(sample_columns)
    if (anyDuplicated(sample_ids)) {
      stop("Covariate file contains duplicate sample IDs: ", path, call. = FALSE)
    }
    values <- numeric_matrix(dt[, -1L, with = FALSE], "Covariate values")
    values <- t(values)
    rownames(values) <- sample_ids
    colnames(values) <- covariate_ids
    values
  })
  sample_ids <- rownames(matrices[[1L]])
  if (anyDuplicated(unlist(lapply(matrices, colnames)))) {
    stop("Covariate files contain duplicate covariate IDs.", call. = FALSE)
  }
  for (i in seq_along(matrices)[-1L]) {
    if (!all(sample_ids %in% rownames(matrices[[i]]))) {
      stop("Covariate files do not contain the same sample IDs.", call. = FALSE)
    }
    matrices[[i]] <- matrices[[i]][sample_ids, , drop = FALSE]
  }
  do.call(cbind, matrices)
}
