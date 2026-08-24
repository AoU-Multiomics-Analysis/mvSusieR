#!/usr/bin/env Rscript

source("scripts/trans_window_model.R")

mu <- array(seq_len(12L) / 10, dim = c(2L, 3L, 2L))
mu2_values <- mu^2 + 0.25
stopifnot(identical(
  get_fit_mu2(list(mu2_diag = mu2_values)),
  mu2_values
))
stopifnot(identical(
  get_fit_mu2(list(mu2 = mu2_values)),
  mu2_values
))
missing_mu2_error <- tryCatch(
  get_fit_mu2(list()),
  error = identity
)
stopifnot(inherits(missing_mu2_error, "error"))

set.seed(20260824)
n <- 100L
p <- 5L
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- cbind(
  IKZF1_expression = 5 * X[, 1L] + rnorm(n, sd = 0.1),
  IKZF1_splicing = -4 * X[, 1L] + rnorm(n, sd = 0.1)
)
colnames(X) <- paste0("chr7:", 50300000L + seq_len(p), "_A_G")
prior <- mvsusieR::create_mixture_prior(R = 2L, null_weight = 0)
fit <- mvsusieR::mvsusie(
  X = X,
  Y = Y,
  L = 2L,
  prior_variance = prior,
  verbose = FALSE
)
stopifnot(length(fit$sets$cs) >= 1L)

prepared <- list(
  X = X,
  Y = Y,
  variant_metadata = data.table::data.table(
    CHROM = "chr7",
    POS = 50300000L + seq_len(p),
    REF = "A",
    ALT = "G",
    variant_id = colnames(X)
  ),
  phenotype_metadata = data.table::data.table(
    outcome_key = colnames(Y),
    phenotype_id = c("IKZF1", "IKZF1:splice_1"),
    modality = c("expression", "splicing")
  )
)
bundle <- list(fit = fit, metadata = list(window_id = "w1"))

output_dir <- tempfile("mvsusie-plot-")
dir.create(output_dir)
prepared_path <- file.path(output_dir, "prepared.rds")
fit_path <- file.path(output_dir, "fit.rds")
png_path <- file.path(output_dir, "credible_set_by_feature.png")
pdf_path <- file.path(output_dir, "credible_set_by_feature.pdf")
plot_rds_path <- file.path(output_dir, "mvsusie_plot_result.rds")
saveRDS(prepared, prepared_path)
saveRDS(bundle, fit_path)

status <- system2(
  "Rscript",
  c(
    "scripts/plot_window_mvsusie.R",
    "--prepared", prepared_path,
    "--fit", fit_path,
    "--png", png_path,
    "--pdf", pdf_path,
    "--plot-rds", plot_rds_path
  )
)
stopifnot(identical(status, 0L))
stopifnot(all(file.exists(c(png_path, pdf_path, plot_rds_path))))
stopifnot(all(file.info(c(png_path, pdf_path, plot_rds_path))$size > 0))
plot_result <- readRDS(plot_rds_path)
stopifnot(!is.null(plot_result$effect_plot))
stopifnot(any(grepl("IKZF1", rownames(plot_result$effects))))

message("mvSuSiE API plotting tests passed")
