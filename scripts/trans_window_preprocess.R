rank_int <- function(x) {
  if (length(x) == 0L) return(numeric())
  qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
}

read_keep_samples <- function(keep_samples) {
  if (is.null(keep_samples)) return(NULL)
  if (length(keep_samples) == 1L && file.exists(keep_samples)) {
    keep_samples <- readLines(keep_samples, warn = FALSE)
  }
  keep_samples <- normalize_sample_ids(keep_samples)
  keep_samples[nzchar(keep_samples)]
}

residualize_matrix <- function(M, covariate_model) {
  qr.resid(covariate_model, M)
}

prepare_window_data <- function(
  window,
  phenotype_data,
  dosage,
  covariates = NULL,
  covariates_by_modality = NULL,
  keep_samples = NULL,
  min_genotype_variance = 1e-8,
  min_phenotype_variance = 1e-8,
  min_nonzero_fraction = NULL
) {
  required_window <- c("window_id", "chrom", "start", "end")
  if (!all(required_window %in% names(window))) {
    stop("Window row is missing required columns.", call. = FALSE)
  }
  if (is.null(phenotype_data$Y) || is.null(dosage$X)) {
    stop("Phenotype and dosage matrices are required.", call. = FALSE)
  }
  if (anyDuplicated(colnames(phenotype_data$Y))) {
    stop("Selected phenotype IDs must be unique.", call. = FALSE)
  }
  if (anyDuplicated(colnames(dosage$X))) {
    stop("Selected variant IDs must be unique.", call. = FALSE)
  }

  if (is.null(covariates_by_modality)) {
    if (is.null(covariates)) {
      stop("Covariates are required.", call. = FALSE)
    }
    covariates_by_modality <- list(shared = covariates)
  }
  if (!is.list(covariates_by_modality) || !length(covariates_by_modality)) {
    stop("Covariates by modality must be a non-empty list.", call. = FALSE)
  }
  phenotype_modalities <- as.character(phenotype_data$modalities)
  if (length(phenotype_modalities) != ncol(phenotype_data$Y)) {
    stop("Phenotype modalities must match the phenotype matrix columns.", call. = FALSE)
  }
  if (is.null(names(covariates_by_modality))) {
    stop("Covariates by modality must be a named list.", call. = FALSE)
  }
  modality_covariates <- lapply(unique(phenotype_modalities), function(modality) {
    matrices <- list()
    if ("shared" %in% names(covariates_by_modality)) {
      matrices <- c(matrices, list(covariates_by_modality[["shared"]]))
    }
    if (modality %in% names(covariates_by_modality)) {
      matrices <- c(matrices, list(covariates_by_modality[[modality]]))
    }
    if (!length(matrices)) {
      stop("No covariates supplied for modality: ", modality, call. = FALSE)
    }
    unique_covariate_columns(matrices)
  })
  names(modality_covariates) <- unique(phenotype_modalities)
  genotype_covariates <- unique_covariate_columns(modality_covariates)

  requested_samples <- read_keep_samples(keep_samples)
  covariate_samples <- Reduce(
    intersect,
    lapply(c(modality_covariates, list(genotype_covariates)), rownames)
  )
  samples <- dosage$sample_ids[
    dosage$sample_ids %in% phenotype_data$sample_ids &
      dosage$sample_ids %in% covariate_samples
  ]
  if (!is.null(requested_samples)) {
    samples <- samples[samples %in% requested_samples]
  }
  if (!length(samples)) {
    stop("No shared samples for window: ", window$window_id, call. = FALSE)
  }

  X_raw <- dosage$X[samples, , drop = FALSE]
  Y_raw <- phenotype_data$Y[samples, , drop = FALSE]
  complete <- apply(is.finite(X_raw), 1L, all) &
    apply(is.finite(Y_raw), 1L, all) &
    apply(is.finite(genotype_covariates[samples, , drop = FALSE]), 1L, all) &
    all(vapply(modality_covariates, function(matrix) {
      all(apply(is.finite(matrix[samples, , drop = FALSE]), 1L, all))
    }, logical(1L)))
  if (!all(complete)) {
    X_raw <- X_raw[complete, , drop = FALSE]
    Y_raw <- Y_raw[complete, , drop = FALSE]
    C <- C[complete, , drop = FALSE]
    samples <- samples[complete]
  }
  if (!length(samples)) {
    stop("No complete samples for window: ", window$window_id, call. = FALSE)
  }

  raw_phenotype_sd <- apply(Y_raw, 2L, sd)
  keep_phenotype <- is.finite(raw_phenotype_sd) &
    raw_phenotype_sd > min_phenotype_variance
  if (!is.null(min_nonzero_fraction)) {
    nonzero_fraction <- colMeans(Y_raw != 0)
    keep_phenotype <- keep_phenotype & nonzero_fraction >= min_nonzero_fraction
  }
  if (!any(keep_phenotype)) {
    stop("No usable phenotypes for window: ", window$window_id, call. = FALSE)
  }
  Y_raw <- Y_raw[, keep_phenotype, drop = FALSE]
  phenotype_ids <- colnames(Y_raw)

  Y_int <- apply(Y_raw, 2L, rank_int)
  Y_int <- as.matrix(Y_int)
  if (is.null(dim(Y_int))) Y_int <- matrix(Y_int, ncol = 1L)
  colnames(Y_int) <- phenotype_ids
  rownames(Y_int) <- samples

  genotype_model <- cbind(genotype_covariates[samples, , drop = FALSE], intercept = 1)
  genotype_qr <- qr(genotype_model)
  X_resid <- residualize_matrix(X_raw, genotype_qr)
  Y_resid <- matrix(NA_real_, nrow = nrow(Y_int), ncol = ncol(Y_int))
  phenotype_covariate_rank <- setNames(integer(length(modality_covariates)), names(modality_covariates))
  for (modality in names(modality_covariates)) {
    phenotype_indices <- which(phenotype_modalities == modality)
    if (!length(phenotype_indices)) next
    phenotype_model <- cbind(
      modality_covariates[[modality]][samples, , drop = FALSE],
      intercept = 1
    )
    phenotype_qr <- qr(phenotype_model)
    Y_resid[, phenotype_indices] <- residualize_matrix(
      Y_int[, phenotype_indices, drop = FALSE],
      phenotype_qr
    )
    phenotype_covariate_rank[[modality]] <- phenotype_qr$rank
  }
  colnames(Y_resid) <- phenotype_ids
  rownames(Y_resid) <- samples

  genotype_sd <- apply(X_resid, 2L, sd)
  keep_variant <- is.finite(genotype_sd) & genotype_sd > min_genotype_variance
  if (!any(keep_variant)) {
    stop("No usable variants for window: ", window$window_id, call. = FALSE)
  }
  X_resid <- X_resid[, keep_variant, drop = FALSE]

  phenotype_resid_sd <- apply(Y_resid, 2L, sd)
  keep_phenotype_resid <- is.finite(phenotype_resid_sd) &
    phenotype_resid_sd > min_phenotype_variance
  if (!any(keep_phenotype_resid)) {
    stop("No usable residualized phenotypes for window: ", window$window_id, call. = FALSE)
  }
  Y_resid <- Y_resid[, keep_phenotype_resid, drop = FALSE]
  phenotype_ids <- phenotype_ids[keep_phenotype_resid]
  Y_scaled <- scale(Y_resid, center = TRUE, scale = TRUE)
  colnames(Y_scaled) <- phenotype_ids
  rownames(Y_scaled) <- samples

  variant_metadata <- data.table::copy(dosage$metadata)[keep_variant]
  variant_metadata[, variant_id := dosage$variant_ids[keep_variant]]
  phenotype_metadata <- data.table::copy(phenotype_data$metadata)[keep_phenotype]
  phenotype_metadata <- phenotype_metadata[keep_phenotype_resid]
  if (!"modality" %in% names(phenotype_metadata)) {
    phenotype_metadata[, modality := phenotype_modalities[keep_phenotype][keep_phenotype_resid]]
  }
  phenotype_metadata[, phenotype_id := phenotype_ids]

  list(
    window = as.list(window),
    X = X_resid,
    Y = Y_scaled,
    variant_metadata = variant_metadata,
    phenotype_metadata = phenotype_metadata,
    samples = samples,
    covariate_rank = genotype_qr$rank,
    phenotype_covariate_rank = phenotype_covariate_rank,
    qc = list(
      window_id = as.character(window$window_id),
      input_samples = length(dosage$sample_ids),
      shared_samples = length(samples),
      input_variants = ncol(dosage$X),
      retained_variants = ncol(X_resid),
      excluded_variants = sum(!keep_variant),
      input_phenotypes = ncol(phenotype_data$Y),
      retained_phenotypes = ncol(Y_scaled),
      excluded_phenotypes = ncol(phenotype_data$Y) - ncol(Y_scaled),
      excluded_samples = length(dosage$sample_ids) - length(samples),
      covariate_rank = genotype_qr$rank,
      phenotype_covariate_rank = phenotype_covariate_rank,
      min_genotype_variance = min_genotype_variance,
      min_phenotype_variance = min_phenotype_variance,
      min_nonzero_fraction = min_nonzero_fraction
    )
  )
}
