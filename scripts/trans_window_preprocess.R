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
  covariates,
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

  requested_samples <- read_keep_samples(keep_samples)
  samples <- dosage$sample_ids[
    dosage$sample_ids %in% phenotype_data$sample_ids &
      dosage$sample_ids %in% rownames(covariates)
  ]
  if (!is.null(requested_samples)) {
    samples <- samples[samples %in% requested_samples]
  }
  if (!length(samples)) {
    stop("No shared samples for window: ", window$window_id, call. = FALSE)
  }

  X_raw <- dosage$X[samples, , drop = FALSE]
  Y_raw <- phenotype_data$Y[samples, , drop = FALSE]
  C <- covariates[samples, , drop = FALSE]
  complete <- apply(is.finite(X_raw), 1L, all) &
    apply(is.finite(Y_raw), 1L, all) &
    apply(is.finite(C), 1L, all)
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

  C_model <- cbind(C, intercept = 1)
  qrC <- qr(C_model)
  X_resid <- residualize_matrix(X_raw, qrC)
  Y_resid <- residualize_matrix(Y_int, qrC)

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
  phenotype_metadata[, phenotype_id := phenotype_ids]

  list(
    window = as.list(window),
    X = X_resid,
    Y = Y_scaled,
    variant_metadata = variant_metadata,
    phenotype_metadata = phenotype_metadata,
    samples = samples,
    covariate_rank = qrC$rank,
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
      covariate_rank = qrC$rank,
      min_genotype_variance = min_genotype_variance,
      min_phenotype_variance = min_phenotype_variance,
      min_nonzero_fraction = min_nonzero_fraction
    )
  )
}
