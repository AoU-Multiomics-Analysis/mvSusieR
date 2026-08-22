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

learn_mashr_prior <- function(Bhat, Shat, n_pca = 5L, seed = NULL) {
  validate_marginal_summary_statistics(Bhat, Shat)
  if (length(n_pca) != 1L || is.na(n_pca) || n_pca < 1L) {
    stop("n_pca must be a positive integer.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  mash_data <- make_mashr_data(Bhat, Shat)
  n_pca <- min(as.integer(n_pca), ncol(Bhat))
  pca_covariances <- mashr::cov_pca(mash_data, npc = n_pca)
  ed_covariances <- mashr::cov_ed(mash_data, Ulist_init = pca_covariances)
  mash_fit <- mashr::mash(
    data = mash_data,
    Ulist = ed_covariances,
    usepointmass = TRUE,
    outputlevel = 0,
    verbose = FALSE
  )
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
