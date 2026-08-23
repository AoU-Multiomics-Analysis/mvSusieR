source("scripts/trans_window_logging.R")

make_model_config <- function(
  L = 10L,
  max_iter = 100L,
  tol = 1e-4,
  coverage = 0.95,
  min_abs_corr = 0.5,
  n_thread = 1L,
  prior_method = "canonical",
  mashr_n_pca = 5L,
  mashr_seed = NULL,
  mashr_strong_lfsr = 0.05,
  mashr_use_ed = TRUE,
  estimate_residual_variance = TRUE,
  marginal_output = NULL
) {
  if (!prior_method %in% c("canonical", "mashr")) {
    stop("prior_method must be either canonical or mashr.", call. = FALSE)
  }
  list(
    L = as.integer(L),
    max_iter = as.integer(max_iter),
    tol = as.numeric(tol),
    coverage = as.numeric(coverage),
    min_abs_corr = as.numeric(min_abs_corr),
    n_thread = as.integer(n_thread),
    prior_method = prior_method,
    mashr_n_pca = as.integer(mashr_n_pca),
    mashr_seed = mashr_seed,
    mashr_strong_lfsr = as.numeric(mashr_strong_lfsr),
    mashr_use_ed = isTRUE(mashr_use_ed),
    estimate_residual_variance = isTRUE(estimate_residual_variance),
    marginal_output = marginal_output
  )
}

make_canonical_prior <- function(n_outcomes) {
  if (length(n_outcomes) != 1L || n_outcomes < 1L) {
    stop("n_outcomes must be a positive integer.", call. = FALSE)
  }
  mvsusieR::create_mixture_prior(R = as.integer(n_outcomes), null_weight = 0)
}

fit_window_mvsusie <- function(prepared, config) {
  if (nrow(prepared$X) != nrow(prepared$Y)) {
    stop("Prepared genotype and phenotype matrices have different sample counts.", call. = FALSE)
  }
  prior_details <- if (identical(config$prior_method, "mashr")) {
    pipeline_log("Preparing all-SNP associations for mashr.")
    marginal <- compute_marginal_bhat_shat_matrix(prepared$X, prepared$Y)
    if (!is.null(config$marginal_output)) {
      write_marginal_association_table(
        marginal$Bhat,
        marginal$Shat,
        config$marginal_output
      )
    }
    learn_mashr_prior(
      Bhat = marginal$Bhat,
      Shat = marginal$Shat,
      n_pca = config$mashr_n_pca,
      seed = config$mashr_seed,
      strong_lfsr = config$mashr_strong_lfsr,
      use_extreme_deconvolution = config$mashr_use_ed
    )
  } else {
    canonical_prior <- make_canonical_prior(ncol(prepared$Y))
    list(
      prior = canonical_prior,
      covariance_training_scope = "not_applicable",
      covariance_training_n = 0L,
      extreme_deconvolution_used = FALSE,
      n_prior_components = length(canonical_prior$xUlist),
      n_covariance_inputs = length(canonical_prior$xUlist)
    )
  }
  prior <- prior_details$prior
  fix_mashr_mixture_weights <- identical(config$prior_method, "mashr")
  prior_scale_conversion <- "not_applicable"
  prior_outcome_se_range <- c(NA_real_, NA_real_)
  if (fix_mashr_mixture_weights) {
    raw_prior_diagonal <- unlist(lapply(prior$xUlist, diag), use.names = FALSE)
    prior <- prepare_mashr_prior_for_mvsusie(prior, prepared$Y)
    prior_outcome_se_range <- range(
      attr(prior, "mvsusie_outcome_se_scale")
    )
    prepared_prior_diagonal <- unlist(
      lapply(prior$xUlist, diag),
      use.names = FALSE
    )
    prior_scale_conversion <- "preserve_mashr_effect_covariance"
    pipeline_log(sprintf(
      "Raw mashr prior diagonal range: %.6g to %.6g.",
      min(raw_prior_diagonal), max(raw_prior_diagonal)
    ))
    pipeline_log(sprintf(
      paste(
        "Pre-scaled mashr prior diagonal range: %.6g to %.6g;",
        "mvSuSiE standardization will restore the raw mashr scale."
      ),
      min(prepared_prior_diagonal), max(prepared_prior_diagonal)
    ))
    pipeline_log("Using the fitted mashr mixture weights without re-estimation.")
  }
  if (!isTRUE(config$estimate_residual_variance)) {
    pipeline_log("Using the initial residual covariance without re-estimation.")
  }
  pipeline_log("Starting mvSuSiE with verbose iteration output.")
  fit <- mvsusieR::mvsusie(
    X = prepared$X,
    Y = prepared$Y,
    L = config$L,
    prior_variance = prior,
    residual_variance = NULL,
    standardize = TRUE,
    intercept = FALSE,
    estimate_residual_variance = config$estimate_residual_variance,
    estimate_prior_variance = FALSE,
    estimate_prior_mixture_weights = !fix_mashr_mixture_weights,
    coverage = config$coverage,
    min_abs_corr = config$min_abs_corr,
    precompute_cache = TRUE,
    n_thread = config$n_thread,
    max_iter = config$max_iter,
    tol = config$tol,
    verbose = TRUE
  )
  if (!isTRUE(fit$converged)) {
    stop(
      "mvsusie did not converge for window: ",
      prepared$qc$window_id,
      call. = FALSE
    )
  }
  pipeline_log(sprintf("mvSuSiE converged after %d iterations.", fit$niter))
  list(
    fit = fit,
    metadata = list(
      window_id = prepared$qc$window_id,
      prior = config$prior_method,
      prior_components = prior_details$n_prior_components,
      prior_covariance_inputs = prior_details$n_covariance_inputs,
      pca_covariance_inputs = prior_details$pca_covariance_inputs,
      mash_model_training_scope = prior_details$mash_model_training_scope,
      mash_model_training_n = prior_details$mash_model_training_n,
      covariance_training_scope = prior_details$covariance_training_scope,
      covariance_training_n = prior_details$covariance_training_n,
      covariance_significant_n = prior_details$covariance_significant_n,
      covariance_selection_lfsr = prior_details$covariance_selection_lfsr,
      covariance_selection_fallback_used =
        prior_details$covariance_selection_fallback_used,
      extreme_deconvolution_used = prior_details$extreme_deconvolution_used,
      covariance_input_method = prior_details$covariance_input_method,
      prior_mixture_weights_mode = if (fix_mashr_mixture_weights) {
        "fixed_from_mashr"
      } else {
        "estimated_by_mvsusie"
      },
      prior_scale_conversion = prior_scale_conversion,
      prior_outcome_se_min = prior_outcome_se_range[[1L]],
      prior_outcome_se_max = prior_outcome_se_range[[2L]],
      residual_variance_mode = if (config$estimate_residual_variance) {
        "estimated_by_mvsusie"
      } else {
        "fixed_initial_covariance"
      },
      mvsusieR_version = as.character(utils::packageVersion("mvsusieR")),
      config = config,
      converged = isTRUE(fit$converged),
      niter = fit$niter
    )
  )
}

empty_credible_set_table <- function() {
  data.table::data.table(
    component = integer(),
    variant_id = character(),
    alpha = numeric(),
    pip = numeric(),
    coverage = numeric(),
    purity_min = numeric(),
    purity_mean = numeric()
  )
}

extract_variant_pips <- function(fit, prepared) {
  data.table::data.table(
    variant_id = colnames(prepared$X),
    pip = as.numeric(fit$pip)
  )
}

extract_credible_sets <- function(fit, prepared, config) {
  cs_obj <- susieR::susie_get_cs(
    fit,
    X = prepared$X,
    coverage = config$coverage,
    min_abs_corr = config$min_abs_corr
  )
  if (!length(cs_obj$cs)) return(empty_credible_set_table())
  purity <- cs_obj$purity
  rows <- lapply(seq_along(cs_obj$cs), function(i) {
    members <- cs_obj$cs[[i]]
    component <- as.integer(sub("^L", "", names(cs_obj$cs)[[i]]))
    purity_row <- if (is.null(dim(purity))) purity else purity[i, ]
    purity_values <- suppressWarnings(as.numeric(unlist(purity_row, use.names = FALSE)))
    data.table::data.table(
      component = component,
      variant_id = colnames(prepared$X)[members],
      alpha = as.numeric(fit$alpha[component, members]),
      pip = as.numeric(fit$pip[members]),
      coverage = config$coverage,
      purity_min = if (length(purity_values)) min(purity_values, na.rm = TRUE) else NA_real_,
      purity_mean = if (length(purity_values)) mean(purity_values, na.rm = TRUE) else NA_real_
    )
  })
  data.table::rbindlist(rows, fill = TRUE)
}

extract_component_effects <- function(fit, prepared) {
  mu <- fit$mu
  mu2 <- fit$mu2
  if (length(dim(mu)) != 3L) {
    stop("Expected mvsusie posterior means with three dimensions.", call. = FALSE)
  }
  dims <- dim(mu)
  if (!identical(dims[2L], ncol(prepared$X)) || !identical(dims[3L], ncol(prepared$Y))) {
    stop("Unexpected mvsusie posterior dimension order.", call. = FALSE)
  }
  idx <- expand.grid(
    component = seq_len(dims[1L]),
    variant_index = seq_len(dims[2L]),
    phenotype_index = seq_len(dims[3L])
  )
  posterior_mean <- as.vector(mu)
  posterior_sd <- sqrt(pmax(as.vector(mu2) - posterior_mean^2, 0))
  data.table::data.table(
    component = idx$component,
    variant_id = colnames(prepared$X)[idx$variant_index],
    phenotype_id = colnames(prepared$Y)[idx$phenotype_index],
    posterior_mean = posterior_mean,
    posterior_sd = posterior_sd
  )
}
