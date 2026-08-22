source("scripts/trans_window_prior.R")

set.seed(20260821)
X <- sweep(matrix(rnorm(80L * 12L), nrow = 80L), 2L, seq_len(12L), "+")
Y <- sweep(matrix(rnorm(80L * 4L), nrow = 80L), 2L, seq_len(4L), "+")
reference <- susieR::compute_marginal_bhat_shat(
  X = scale(X, center = TRUE, scale = TRUE),
  Y = Y
)
marginal_messages <- capture.output(
  marginal <- compute_marginal_bhat_shat_matrix(X, Y),
  type = "message"
)
stopifnot(isTRUE(all.equal(marginal$Bhat, reference$Bhat, tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(marginal$Shat, reference$Shat, tolerance = 1e-10)))
stopifnot(any(grepl("cross-product", marginal_messages, fixed = TRUE)))
stopifnot(any(grepl("Marginal associations complete", marginal_messages, fixed = TRUE)))

set.seed(20260822)
Bhat <- matrix(rnorm(72L * 3L), nrow = 72L, ncol = 3L)
Shat <- matrix(runif(72L * 3L, min = 0.05, max = 0.2), nrow = 72L, ncol = 3L)

prior_messages <- capture.output(
  prior_fit <- learn_mashr_prior(
    Bhat = Bhat,
    Shat = Shat,
    n_pca = 2L,
    seed = 1L
  ),
  type = "message"
)

stopifnot(inherits(prior_fit$prior, "mash_prior"))
stopifnot(identical(prior_fit$covariance_training_scope, "all_snps_in_window"))
stopifnot(identical(prior_fit$covariance_training_n, nrow(Bhat)))
stopifnot(isTRUE(prior_fit$extreme_deconvolution_used))
stopifnot(prior_fit$n_covariance_inputs >= 1L)
stopifnot(any(grepl("PCA covariance", prior_messages, fixed = TRUE)))
stopifnot(any(grepl("extreme deconvolution", prior_messages, fixed = TRUE)))
stopifnot(any(grepl("mashr mixture", prior_messages, fixed = TRUE)))

message("All-SNP mashr prior tests passed")
