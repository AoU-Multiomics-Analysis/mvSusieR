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

compute_marginal_bhat_shat_matrix <- function(X, Y) {
  if (!is.matrix(X) || !is.numeric(X) || !is.matrix(Y) || !is.numeric(Y)) {
    stop("X and Y must be numeric matrices.", call. = FALSE)
  }
  if (nrow(X) != nrow(Y) || nrow(X) < 2L) {
    stop("X and Y must have the same number of at least two rows.", call. = FALSE)
  }
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("X and Y must contain only finite values.", call. = FALSE)
  }

  start_time <- proc.time()[["elapsed"]]
  n <- nrow(X)
  pipeline_log(sprintf(
    "Starting marginal associations: N=%d, J=%d, R=%d.",
    n, ncol(X), ncol(Y)
  ))
  pipeline_log("Computing centered sums of squares.")
  x_mean <- colMeans(X)
  y_mean <- colMeans(Y)
  x_ss <- colSums(X^2) - n * x_mean^2
  y_ss <- colSums(Y^2) - n * y_mean^2
  if (any(!is.finite(x_ss)) || any(x_ss <= 0)) {
    stop("X contains a non-finite or zero-variance column.", call. = FALSE)
  }

  pipeline_log("Computing the all-SNP cross-product.")
  xy_centered <- crossprod(X, Y) - n * outer(x_mean, y_mean)
  x_scale <- sqrt(x_ss / (n - 1))
  xy_standardized <- sweep(xy_centered, 1L, x_scale, "/")
  predictor_weight <- n - 1
  Bhat <- xy_standardized / predictor_weight

  explained_ss <- xy_standardized^2 / predictor_weight
  residual_ss <- matrix(
    y_ss,
    nrow = nrow(explained_ss),
    ncol = ncol(explained_ss),
    byrow = TRUE
  ) - explained_ss
  residual_variance <- pmax(residual_ss / predictor_weight, 1e-64)
  Shat <- sqrt(residual_variance) / sqrt(predictor_weight)
  validate_marginal_summary_statistics(Bhat, Shat)

  pipeline_log(sprintf(
    "Marginal associations complete in %.2f seconds.",
    proc.time()[["elapsed"]] - start_time
  ))
  list(Bhat = Bhat, Shat = Shat)
}

learn_mashr_prior <- function(Bhat, Shat, n_pca = 5L, seed = NULL) {
  validate_marginal_summary_statistics(Bhat, Shat)
  if (length(n_pca) != 1L || is.na(n_pca) || n_pca < 1L) {
    stop("n_pca must be a positive integer.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  mash_data <- make_mashr_data(Bhat, Shat)
  n_pca <- min(as.integer(n_pca), ncol(Bhat))
  stage_time <- proc.time()[["elapsed"]]
  pipeline_log(sprintf("Starting %d PCA covariance inputs.", n_pca))
  pca_covariances <- mashr::cov_pca(mash_data, npc = n_pca)
  pipeline_log(sprintf(
    "PCA covariance inputs complete in %.2f seconds.",
    proc.time()[["elapsed"]] - stage_time
  ))

  stage_time <- proc.time()[["elapsed"]]
  pipeline_log("Starting extreme deconvolution.")
  ed_covariances <- mashr::cov_ed(mash_data, Ulist_init = pca_covariances)
  pipeline_log(sprintf(
    "Extreme deconvolution complete in %.2f seconds.",
    proc.time()[["elapsed"]] - stage_time
  ))

  stage_time <- proc.time()[["elapsed"]]
  pipeline_log("Starting the mashr mixture fit.")
  mash_fit <- mashr::mash(
    data = mash_data,
    Ulist = ed_covariances,
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
  fallback_to_ed_covariances <- FALSE
  if (!length(prior$xUlist)) {
    prior <- mvsusieR::create_mixture_prior(
      mixture_prior = list(matrices = ed_covariances),
      null_weight = 0,
      weights_tol = 0
    )
    fallback_to_ed_covariances <- TRUE
  }

  list(
    prior = prior,
    Bhat = Bhat,
    Shat = Shat,
    covariance_training_scope = "all_snps_in_window",
    covariance_training_n = nrow(Bhat),
    pca_covariance_inputs = length(pca_covariances),
    extreme_deconvolution_used = TRUE,
    n_covariance_inputs = length(ed_covariances),
    n_prior_components = length(prior$xUlist),
    fallback_to_ed_covariances = fallback_to_ed_covariances,
    fitted_g = mash_fit$fitted_g
  )
}
