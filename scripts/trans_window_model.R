make_model_config <- function(
  L = 10L,
  max_iter = 100L,
  tol = 1e-4,
  coverage = 0.95,
  min_abs_corr = 0.5,
  n_thread = 1L
) {
  list(
    L = as.integer(L),
    max_iter = as.integer(max_iter),
    tol = as.numeric(tol),
    coverage = as.numeric(coverage),
    min_abs_corr = as.numeric(min_abs_corr),
    n_thread = as.integer(n_thread)
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
  prior <- make_canonical_prior(ncol(prepared$Y))
  fit <- mvsusieR::mvsusie(
    X = prepared$X,
    Y = prepared$Y,
    L = config$L,
    prior_variance = prior,
    residual_variance = NULL,
    standardize = TRUE,
    intercept = FALSE,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = FALSE,
    coverage = config$coverage,
    min_abs_corr = config$min_abs_corr,
    precompute_cache = TRUE,
    n_thread = config$n_thread,
    max_iter = config$max_iter,
    tol = config$tol,
    verbose = FALSE
  )
  if (!isTRUE(fit$converged)) {
    stop(
      "mvsusie did not converge for window: ",
      prepared$qc$window_id,
      call. = FALSE
    )
  }
  list(
    fit = fit,
    metadata = list(
      window_id = prepared$qc$window_id,
      prior = "canonical",
      residual_variance_mode = "mvsusieR_default",
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
