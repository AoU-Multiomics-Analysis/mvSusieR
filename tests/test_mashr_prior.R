source("scripts/trans_window_prior.R")

set.seed(20260822)
Bhat <- matrix(rnorm(72L * 3L), nrow = 72L, ncol = 3L)
Shat <- matrix(runif(72L * 3L, min = 0.05, max = 0.2), nrow = 72L, ncol = 3L)

prior_fit <- learn_mashr_prior(
  Bhat = Bhat,
  Shat = Shat,
  n_pca = 2L,
  seed = 1L
)

stopifnot(inherits(prior_fit$prior, "mash_prior"))
stopifnot(identical(prior_fit$covariance_training_scope, "all_snps_in_window"))
stopifnot(identical(prior_fit$covariance_training_n, nrow(Bhat)))
stopifnot(isTRUE(prior_fit$extreme_deconvolution_used))
stopifnot(prior_fit$n_covariance_inputs >= 1L)

message("All-SNP mashr prior tests passed")
