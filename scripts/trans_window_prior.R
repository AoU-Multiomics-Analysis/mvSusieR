source("scripts/trans_window_logging.R")

validate_marginal_summary_statistics <- function(Bhat, Shat) {
  if (!is.matrix(Bhat) || !is.numeric(Bhat)) {
    stop("Bhat must be a numeric matrix.", call. = FALSE)
  }
  if (!is.matrix(Shat) || !is.numeric(Shat)) {
    stop("Shat must be a numeric matrix.", call. = FALSE)
  }
  if (!identical(dim(Bhat), dim(Shat))) {
    stop("Bhat and Shat must have identical dimensions.", call. = FALSE)
  }
  if (any(!is.finite(Bhat)) || any(!is.finite(Shat))) {
    stop("Bhat and Shat must contain only finite values.", call. = FALSE)
  }
  if (any(Shat <= 0)) {
    stop("Shat must contain strictly positive values.", call. = FALSE)
  }
  invisible(TRUE)
}

make_mashr_data <- function(Bhat, Shat) {
  validate_marginal_summary_statistics(Bhat, Shat)
  if (!requireNamespace("mashr", quietly = TRUE)) {
    stop("The mashr package is required to learn a data-driven prior.", call. = FALSE)
  }
  mashr::mash_set_data(Bhat = Bhat, Shat = Shat, alpha = 0)
}

prepare_mashr_prior_for_mvsusie <- function(prior, Y) {
  if (!inherits(prior, "mash_prior") || !length(prior$xUlist)) {
    stop("prior must be a non-empty mash prior.", call. = FALSE)
  }
  if (!is.matrix(Y) || !is.numeric(Y) || ncol(Y) < 1L) {
    stop("Y must be a numeric matrix with at least one outcome.", call. = FALSE)
  }
  outcome_n <- colSums(is.finite(Y))
  outcome_sd <- apply(Y, 2L, stats::sd, na.rm = TRUE)
  outcome_se <- outcome_sd / sqrt(outcome_n)
  if (any(!is.finite(outcome_se)) || any(outcome_se <= 0)) {
    stop("Each outcome must have a finite positive standard error scale.", call. = FALSE)
  }
  if (any(vapply(prior$xUlist, nrow, integer(1L)) != ncol(Y))) {
    stop("Mash prior matrices must match the number of outcomes.", call. = FALSE)
  }

  automatic_scale <- tcrossprod(outcome_se)
  prior$xUlist <- lapply(prior$xUlist, function(U) U / automatic_scale)
  attr(prior, "mvsusie_outcome_se_scale") <- outcome_se
  prior
}

compute_marginal_bhat_shat_matrix <- function(X, Y, block_size = 1000L) {
  if (!is.matrix(X) || !is.numeric(X) || !is.matrix(Y) || !is.numeric(Y)) {
    stop("X and Y must be numeric matrices.", call. = FALSE)
  }
  if (nrow(X) != nrow(Y) || nrow(X) < 2L) {
    stop("X and Y must have the same number of at least two rows.", call. = FALSE)
  }
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("X and Y must contain only finite values.", call. = FALSE)
  }
  if (
    length(block_size) != 1L || is.na(block_size) ||
    block_size < 1L || block_size != as.integer(block_size)
  ) {
    stop("block_size must be a positive integer.", call. = FALSE)
  }

  start_time <- proc.time()[["elapsed"]]
  n <- nrow(X)
  n_variants <- ncol(X)
  n_outcomes <- ncol(Y)
  pipeline_log(sprintf(
    "Starting marginal associations: N=%d, J=%d, R=%d.",
    n, n_variants, n_outcomes
  ))
  pipeline_log("Centering outcomes and computing their sums of squares.")
  Y_centered <- sweep(Y, 2L, colMeans(Y), "-")
  y_ss <- colSums(Y_centered^2)
  predictor_weight <- n - 1
  Bhat <- matrix(NA_real_, nrow = n_variants, ncol = n_outcomes)
  Shat <- matrix(NA_real_, nrow = n_variants, ncol = n_outcomes)
  if (!is.null(colnames(X))) rownames(Bhat) <- rownames(Shat) <- colnames(X)
  if (!is.null(colnames(Y))) colnames(Bhat) <- colnames(Shat) <- colnames(Y)
  block_starts <- seq.int(1L, n_variants, by = as.integer(block_size))

  for (block_index in seq_along(block_starts)) {
    first <- block_starts[[block_index]]
    last <- min(first + as.integer(block_size) - 1L, n_variants)
    indices <- first:last
    pipeline_log(sprintf(
      "Computing the all-SNP cross-product block %d of %d (%d-%d).",
      block_index, length(block_starts), first, last
    ))
    X_block <- X[, indices, drop = FALSE]
    X_block <- sweep(X_block, 2L, colMeans(X_block), "-")
    x_ss <- colSums(X_block^2)
    if (any(!is.finite(x_ss)) || any(x_ss <= 0)) {
      stop("X contains a non-finite or zero-variance column.", call. = FALSE)
    }
    x_scale <- sqrt(x_ss / predictor_weight)
    X_standardized <- sweep(X_block, 2L, x_scale, "/")
    block_Bhat <- crossprod(X_standardized, Y) / predictor_weight
    Bhat[indices, ] <- block_Bhat

    X_standardized_centered <- sweep(
      X_standardized,
      2L,
      colMeans(X_standardized),
      "-"
    )
    residual_crossproduct <- crossprod(X_standardized_centered, Y_centered)
    standardized_x_ss <- colSums(X_standardized_centered^2)
    residual_ss <- matrix(
      y_ss,
      nrow = length(indices),
      ncol = n_outcomes,
      byrow = TRUE
    ) - 2 * block_Bhat * residual_crossproduct +
      block_Bhat^2 * standardized_x_ss

    y_ss_matrix <- matrix(
      y_ss,
      nrow = length(indices),
      ncol = n_outcomes,
      byrow = TRUE
    )
    unstable <- !is.finite(residual_ss) |
      residual_ss <= sqrt(.Machine$double.eps) * y_ss_matrix
    if (any(unstable)) {
      pipeline_log(sprintf(
        paste(
          "Recomputing %d near-perfect association residual sums of squares",
          "directly."
        ),
        sum(unstable)
      ))
      fallback_chunk_size <- 100L
      for (outcome_index in which(colSums(unstable) > 0L)) {
        local_indices <- which(unstable[, outcome_index])
        fallback_starts <- seq.int(
          1L,
          length(local_indices),
          by = fallback_chunk_size
        )
        for (fallback_first in fallback_starts) {
          fallback_last <- min(
            fallback_first + fallback_chunk_size - 1L,
            length(local_indices)
          )
          selected <- local_indices[fallback_first:fallback_last]
          fitted_values <- sweep(
            X_standardized_centered[, selected, drop = FALSE],
            2L,
            block_Bhat[selected, outcome_index],
            "*"
          )
          residuals <- -sweep(
            fitted_values,
            1L,
            Y_centered[, outcome_index],
            "-"
          )
          residual_ss[selected, outcome_index] <- colSums(residuals^2)
        }
      }
    }
    residual_variance <- pmax(residual_ss / predictor_weight, 1e-64)
    Shat[indices, ] <- sqrt(residual_variance) / sqrt(predictor_weight)
  }
  validate_marginal_summary_statistics(Bhat, Shat)

  pipeline_log(sprintf(
    "Marginal associations complete in %.2f seconds.",
    proc.time()[["elapsed"]] - start_time
  ))
  list(Bhat = Bhat, Shat = Shat)
}

write_marginal_association_table <- function(Bhat, Shat, path) {
  validate_marginal_summary_statistics(Bhat, Shat)
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("path must be one non-empty file path.", call. = FALSE)
  }
  variant_ids <- rownames(Bhat)
  if (is.null(variant_ids)) variant_ids <- paste0("variant_", seq_len(nrow(Bhat)))
  feature_ids <- colnames(Bhat)
  if (is.null(feature_ids)) feature_ids <- paste0("feature_", seq_len(ncol(Bhat)))

  pipeline_log(sprintf(
    "Writing the complete marginal association table with %d rows.",
    length(Bhat)
  ))
  stage_time <- proc.time()[["elapsed"]]
  output_parent <- dirname(path)
  if (!dir.exists(output_parent)) {
    dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
  }
  z <- as.vector(Bhat) / as.vector(Shat)
  associations <- data.table::data.table(
    variant_id = rep(variant_ids, times = ncol(Bhat)),
    feature_id = rep(feature_ids, each = nrow(Bhat)),
    bhat = as.vector(Bhat),
    shat = as.vector(Shat),
    z = z,
    p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  )
  data.table::fwrite(associations, path, sep = "\t", compress = "auto")
  pipeline_log(sprintf(
    "The marginal association table was saved in %.2f seconds.",
    proc.time()[["elapsed"]] - stage_time
  ))
  invisible(path)
}

select_mashr_covariance_rows <- function(mash_data, n_pca, lfsr_threshold) {
  stage_time <- proc.time()[["elapsed"]]
  pipeline_log(sprintf(
    "Starting the one-by-one mashr fit for strong-row selection at lfsr <= %.3g.",
    lfsr_threshold
  ))
  one_by_one <- mashr::mash_1by1(mash_data)
  strong_rows <- mashr::get_significant_results(
    one_by_one,
    thresh = lfsr_threshold
  )
  significant_n <- length(strong_rows)
  fallback_used <- FALSE
  if (significant_n < n_pca) {
    fallback_used <- TRUE
    minimum_lfsr <- apply(one_by_one$result$lfsr, 1L, min)
    strong_rows <- head(order(minimum_lfsr), n_pca)
    pipeline_log(sprintf(
      paste(
        "Only %d rows passed the strong-row threshold; using the %d rows",
        "with the smallest lfsr values for covariance training."
      ),
      significant_n, length(strong_rows)
    ))
  }
  pipeline_log(sprintf(
    paste(
      "Selected %d strong SNP rows from %d total rows in %.2f seconds",
      "for PCA and extreme deconvolution."
    ),
    length(strong_rows), nrow(mash_data$Bhat),
    proc.time()[["elapsed"]] - stage_time
  ))
  list(
    rows = as.integer(strong_rows),
    significant_n = significant_n,
    fallback_used = fallback_used
  )
}

learn_mashr_prior <- function(
  Bhat,
  Shat,
  n_pca = 5L,
  seed = NULL,
  strong_lfsr = 0.05,
  use_extreme_deconvolution = TRUE
) {
  validate_marginal_summary_statistics(Bhat, Shat)
  if (length(n_pca) != 1L || is.na(n_pca) || n_pca < 2L) {
    stop("n_pca must be an integer of at least two.", call. = FALSE)
  }
  if (
    length(strong_lfsr) != 1L || is.na(strong_lfsr) ||
    strong_lfsr <= 0 || strong_lfsr >= 1
  ) {
    stop("strong_lfsr must be between zero and one.", call. = FALSE)
  }
  if (
    length(use_extreme_deconvolution) != 1L ||
    is.na(use_extreme_deconvolution)
  ) {
    stop("use_extreme_deconvolution must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  mash_data <- make_mashr_data(Bhat, Shat)
  n_pca <- min(as.integer(n_pca), ncol(Bhat), nrow(Bhat))
  covariance_selection <- select_mashr_covariance_rows(
    mash_data,
    n_pca = n_pca,
    lfsr_threshold = strong_lfsr
  )
  covariance_rows <- covariance_selection$rows
  stage_time <- proc.time()[["elapsed"]]
  pipeline_log(sprintf(
    "Starting %d PCA covariance inputs on %d selected rows.",
    n_pca, length(covariance_rows)
  ))
  pca_covariances <- mashr::cov_pca(
    mash_data,
    npc = n_pca,
    subset = covariance_rows
  )
  pipeline_log(sprintf(
    "PCA covariance inputs complete in %.2f seconds.",
    proc.time()[["elapsed"]] - stage_time
  ))

  covariance_inputs <- pca_covariances
  if (isTRUE(use_extreme_deconvolution)) {
    stage_time <- proc.time()[["elapsed"]]
    pipeline_log(sprintf(
      "Starting extreme deconvolution on %d selected rows.",
      length(covariance_rows)
    ))
    covariance_inputs <- mashr::cov_ed(
      mash_data,
      Ulist_init = pca_covariances,
      subset = covariance_rows
    )
    pipeline_log(sprintf(
      "Extreme deconvolution complete in %.2f seconds.",
      proc.time()[["elapsed"]] - stage_time
    ))
  } else {
    pipeline_log("Skipping extreme deconvolution; using PCA covariance inputs.")
  }

  stage_time <- proc.time()[["elapsed"]]
  pipeline_log(sprintf(
    "Starting the mashr mixture fit on all %d SNP rows.",
    nrow(Bhat)
  ))
  mash_fit <- mashr::mash(
    data = mash_data,
    Ulist = covariance_inputs,
    usepointmass = TRUE,
    outputlevel = 0,
    verbose = TRUE
  )
  pipeline_log(sprintf(
    "The mashr mixture fit completed in %.2f seconds.",
    proc.time()[["elapsed"]] - stage_time
  ))
  prior <- mvsusieR::create_mixture_prior(
    fitted_g = mash_fit$fitted_g,
    null_weight = 0
  )
  fallback_to_input_covariances <- FALSE
  if (!length(prior$xUlist)) {
    prior <- mvsusieR::create_mixture_prior(
      mixture_prior = list(matrices = covariance_inputs),
      null_weight = 0,
      weights_tol = 0
    )
    fallback_to_input_covariances <- TRUE
  }

  list(
    prior = prior,
    Bhat = Bhat,
    Shat = Shat,
    mash_model_training_scope = "all_snps_in_window",
    mash_model_training_n = nrow(Bhat),
    covariance_training_scope = "strong_snps_in_window",
    covariance_training_n = length(covariance_rows),
    covariance_significant_n = covariance_selection$significant_n,
    covariance_selection_lfsr = strong_lfsr,
    covariance_selection_fallback_used = covariance_selection$fallback_used,
    pca_covariance_inputs = length(pca_covariances),
    extreme_deconvolution_used = isTRUE(use_extreme_deconvolution),
    covariance_input_method = if (isTRUE(use_extreme_deconvolution)) {
      "pca_then_extreme_deconvolution"
    } else {
      "pca_only"
    },
    n_covariance_inputs = length(covariance_inputs),
    n_prior_components = length(prior$xUlist),
    fallback_to_input_covariances = fallback_to_input_covariances,
    fitted_g = mash_fit$fitted_g
  )
}
