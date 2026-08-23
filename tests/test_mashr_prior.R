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
Bhat <- matrix(rnorm(72L * 3L), nrow = 72L, ncol = 3L)
Shat <- matrix(runif(72L * 3L, min = 0.05, max = 0.2), nrow = 72L, ncol = 3L)
rownames(Bhat) <- rownames(Shat) <- paste0("variant_", seq_len(nrow(Bhat)))
colnames(Bhat) <- colnames(Shat) <- paste0("feature_", seq_len(ncol(Bhat)))

association_path <- tempfile(fileext = ".tsv.gz")
association_messages <- capture.output(
  write_marginal_association_table(Bhat, Shat, association_path),
  type = "message"
)
associations <- data.table::fread(association_path, check.names = FALSE)
stopifnot(nrow(associations) == length(Bhat))
stopifnot(identical(
  names(associations),
  c("variant_id", "feature_id", "bhat", "shat", "z", "p_value")
))
stopifnot(any(grepl("association table", association_messages, fixed = TRUE)))

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
stopifnot(identical(prior_fit$mash_model_training_scope, "all_snps_in_window"))
stopifnot(identical(prior_fit$mash_model_training_n, nrow(Bhat)))
stopifnot(identical(prior_fit$covariance_training_scope, "strong_snps_in_window"))
stopifnot(prior_fit$covariance_training_n <= nrow(Bhat))
stopifnot(prior_fit$covariance_training_n >= 2L)
stopifnot(isTRUE(prior_fit$extreme_deconvolution_used))
stopifnot(prior_fit$n_covariance_inputs >= 1L)
stopifnot(any(grepl("one-by-one", prior_messages, fixed = TRUE)))
stopifnot(any(grepl("strong SNP rows", prior_messages, fixed = TRUE)))
stopifnot(any(grepl("PCA covariance", prior_messages, fixed = TRUE)))
stopifnot(any(grepl("extreme deconvolution", prior_messages, fixed = TRUE)))
stopifnot(any(grepl("mashr mixture", prior_messages, fixed = TRUE)))

pca_only_messages <- capture.output(
  pca_only_prior <- learn_mashr_prior(
    Bhat = Bhat,
    Shat = Shat,
    n_pca = 2L,
    seed = 1L,
    use_extreme_deconvolution = FALSE
  ),
  type = "message"
)
stopifnot(!isTRUE(pca_only_prior$extreme_deconvolution_used))
stopifnot(identical(pca_only_prior$covariance_input_method, "pca_only"))
stopifnot(any(grepl("Skipping extreme deconvolution", pca_only_messages, fixed = TRUE)))
stopifnot(!any(grepl("Starting extreme deconvolution", pca_only_messages, fixed = TRUE)))

message("All-SNP mashr prior tests passed")
