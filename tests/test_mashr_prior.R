#!/usr/bin/env Rscript

source("scripts/trans_window_prior.R")

scale_test_prior <- mvsusieR::create_mixture_prior(
  mixture_prior = list(
    matrices = list(matrix(c(0.04, 0.03, 0.03, 0.09), 2L, 2L)),
    weights = 1
  ),
  null_weight = 0
)
scale_test_Y <- cbind(c(-1, 0, 1), c(-2, 0, 2))
scale_test_prepared <- prepare_mashr_prior_for_mvsusie(
  scale_test_prior,
  scale_test_Y
)
scale_test_expected <- matrix(c(0.12, 0.045, 0.045, 0.0675), 2L, 2L)
stopifnot(isTRUE(all.equal(
  scale_test_prepared$xUlist[[1L]],
  scale_test_expected,
  tolerance = 1e-12
)))
scale_test_sigma <- c(1 / sqrt(3), 2 / sqrt(3))
scale_test_round_trip <-
  t(scale_test_prepared$xUlist[[1L]] * scale_test_sigma) * scale_test_sigma
stopifnot(isTRUE(all.equal(
  scale_test_round_trip,
  scale_test_prior$xUlist[[1L]],
  tolerance = 1e-12
)))

set.seed(20260821)
X <- sweep(matrix(rnorm(80L * 12L), nrow = 80L), 2L, seq_len(12L), "+")
Y <- sweep(matrix(rnorm(80L * 6L), nrow = 80L), 2L, seq_len(6L), "+")
reference <- susieR::compute_marginal_bhat_shat(
  X = scale(X, center = TRUE, scale = TRUE),
  Y = Y
)
marginal_messages <- capture.output(
  marginal <- compute_marginal_bhat_shat_matrix(X, Y, block_size = 5L),
  type = "message"
)
stopifnot(isTRUE(all.equal(marginal$Bhat, reference$Bhat, tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(marginal$Shat, reference$Shat, tolerance = 1e-10)))
stopifnot(sum(grepl("cross-product block", marginal_messages, fixed = TRUE)) == 3L)

X_large_offset <- X + 1e12
large_offset_reference <- susieR::compute_marginal_bhat_shat(
  X = scale(X_large_offset, center = TRUE, scale = TRUE),
  Y = Y
)
large_offset_marginal <- suppressMessages(
  compute_marginal_bhat_shat_matrix(X_large_offset, Y)
)
stopifnot(isTRUE(all.equal(
  large_offset_marginal$Bhat,
  large_offset_reference$Bhat,
  tolerance = 1e-8
)))
stopifnot(isTRUE(all.equal(
  large_offset_marginal$Shat,
  large_offset_reference$Shat,
  tolerance = 1e-8
)))

set.seed(20260820)
X_near_perfect <- matrix(rnorm(100L * 3L), nrow = 100L)
X_near_perfect_scaled <- scale(X_near_perfect)
Y_near_perfect <- cbind(
  2 * X_near_perfect_scaled[, 1L] + rnorm(100L, sd = 1e-10),
  -3 * X_near_perfect_scaled[, 2L] + rnorm(100L, sd = 1e-8)
)
near_perfect_reference <- susieR::compute_marginal_bhat_shat(
  X = X_near_perfect_scaled,
  Y = Y_near_perfect
)
near_perfect_marginal <- suppressMessages(
  compute_marginal_bhat_shat_matrix(X_near_perfect, Y_near_perfect)
)
relative_shat_error <- abs(
  near_perfect_marginal$Shat - near_perfect_reference$Shat
) / near_perfect_reference$Shat
stopifnot(max(relative_shat_error) < 1e-6)

set.seed(20260822)
Bhat <- matrix(rnorm(72L * 6L), nrow = 72L, ncol = 6L)
Shat <- matrix(runif(72L * 6L, min = 0.05, max = 0.2), nrow = 72L, ncol = 6L)
rownames(Bhat) <- rownames(Shat) <- paste0("variant_", seq_len(nrow(Bhat)))
colnames(Bhat) <- colnames(Shat) <- paste0("feature_", seq_len(ncol(Bhat)))

observed_pca <- new.env(parent = emptyenv())
recording_cov_pca <- function(data, npc, subset) {
  observed_pca$npc <- npc
  observed_pca$subset <- subset
  mashr::cov_pca(data, npc = npc, subset = subset)
}
prior_messages <- capture.output(
  prior_fit <- learn_joint_mashr_prior(
    Bhat = Bhat,
    Shat = Shat,
    n_pca = 5L,
    seed = 1L,
    strong_lfsr = 0.05,
    cov_pca_fun = recording_cov_pca
  ),
  type = "message"
)

stopifnot(identical(observed_pca$npc, 5L))
stopifnot(length(observed_pca$subset) >= 5L)
stopifnot(inherits(prior_fit$raw_prior, "mash_prior"))
stopifnot(identical(prior_fit$mash_model_training_scope, "all_snps_in_window"))
stopifnot(identical(prior_fit$mash_model_training_n, 72L))
stopifnot(identical(prior_fit$covariance_training_scope, "strong_snps_in_window"))
stopifnot(identical(prior_fit$pca_requested, 5L))
stopifnot(prior_fit$pca_returned >= 1L)
stopifnot(identical(prior_fit$covariance_input_method, "pca_only"))
stopifnot(identical(prior_fit$Bhat, Bhat))
stopifnot(identical(prior_fit$Shat, Shat))
stopifnot(length(prior_fit$fitted_weights) >= 1L)
stopifnot(abs(sum(prior_fit$fitted_weights) - 1) < 1e-8)
stopifnot(any(grepl("one-by-one", prior_messages, fixed = TRUE)))
stopifnot(any(grepl("PCA covariance", prior_messages, fixed = TRUE)))
stopifnot(any(grepl("mashr mixture", prior_messages, fixed = TRUE)))
stopifnot(!any(grepl("extreme deconvolution", prior_messages, fixed = TRUE)))

set.seed(20260824)
small_Bhat <- matrix(rnorm(24L), nrow = 12L, ncol = 2L)
small_Shat <- matrix(runif(24L, min = 0.05, max = 0.2), nrow = 12L, ncol = 2L)
rownames(small_Bhat) <- rownames(small_Shat) <- paste0("small_variant_", 1:12)
colnames(small_Bhat) <- colnames(small_Shat) <- c("expression::target", "protein::p1")
small_observed <- new.env(parent = emptyenv())
small_prior_fit <- suppressMessages(learn_joint_mashr_prior(
  Bhat = small_Bhat,
  Shat = small_Shat,
  n_pca = 5L,
  seed = 1L,
  cov_pca_fun = function(data, npc, subset) {
    small_observed$npc <- npc
    mashr::cov_pca(data, npc = npc, subset = subset)
  }
))
stopifnot(identical(small_observed$npc, 2L))
stopifnot(identical(small_prior_fit$pca_requested, 5L))
stopifnot(identical(small_prior_fit$pca_used, 2L))

univariate_Bhat <- matrix(
  rnorm(12L), nrow = 12L, ncol = 1L,
  dimnames = list(paste0("univariate_variant_", 1:12), "expression::target")
)
univariate_Shat <- matrix(
  runif(12L, min = 0.05, max = 0.2), nrow = 12L, ncol = 1L,
  dimnames = dimnames(univariate_Bhat)
)
univariate_prior_fit <- suppressMessages(learn_joint_mashr_prior(
  Bhat = univariate_Bhat,
  Shat = univariate_Shat,
  n_pca = 5L,
  seed = 1L,
  cov_pca_fun = function(...) stop("cov_pca must not run for one outcome")
))
stopifnot(identical(univariate_prior_fit$pca_requested, 5L))
stopifnot(identical(univariate_prior_fit$pca_used, 1L))
stopifnot(identical(
  univariate_prior_fit$covariance_input_method,
  "univariate_pca_equivalent"
))
stopifnot(all(vapply(
  univariate_prior_fit$raw_prior$xUlist,
  function(U) identical(dim(U), c(1L, 1L)),
  logical(1L)
)))

bad_covariance <- list(diag(6L))
bad_covariance[[1L]][1L, 1L] <- Inf
bad_covariance_error <- tryCatch(
  validate_mashr_covariances(bad_covariance, 6L),
  error = identity
)
stopifnot(inherits(bad_covariance_error, "error"))

wrong_dimension_error <- tryCatch(
  validate_mashr_covariances(list(diag(5L)), 6L),
  error = identity
)
stopifnot(inherits(wrong_dimension_error, "error"))

message("Joint PCA-only mashr prior tests passed")
