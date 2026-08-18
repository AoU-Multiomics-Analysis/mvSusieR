source("scripts/trans_window_io.R")

windows <- read_windows_manifest("tests/fixtures/trans_window/windows.tsv")
stopifnot(nrow(windows) == 1L, windows$window_id == "w1")

phenotype_manifest <- read_window_phenotypes_manifest(
  "tests/fixtures/trans_window/window_phenotypes.tsv"
)
stopifnot(nrow(phenotype_manifest) == 3L)

dosage <- read_wide_dosage("tests/fixtures/trans_window/window_1_dosage.tsv")
stopifnot(identical(dim(dosage$X), c(6L, 2L)))
stopifnot(identical(dosage$sample_ids, as.character(1:6)))
stopifnot(identical(dosage$variant_ids, c("chr1:101_A_G", "chr1:202_C_T")))

phenotypes <- read_phenotype_rows(
  "tests/fixtures/trans_window/expression.tsv",
  "expression",
  "ENSG000001.1"
)
stopifnot(ncol(phenotypes$Y) == 1L, nrow(phenotypes$Y) == 6L)
stopifnot(identical(colnames(phenotypes$Y), "ENSG000001.1"))

covariates <- read_covariate_matrix(
  "tests/fixtures/trans_window/covariates.tsv"
)
stopifnot(identical(dim(covariates), c(6L, 2L)))
stopifnot(identical(rownames(covariates), as.character(1:6)))

source("scripts/trans_window_preprocess.R")

phenotype_data <- read_window_phenotypes(
  window_id = "w1",
  phenotype_manifest = phenotype_manifest,
  phenotype_files = c(
    "tests/fixtures/trans_window/expression.tsv",
    "tests/fixtures/trans_window/splicing.tsv",
    "tests/fixtures/trans_window/isoform_usage.tsv"
  )
)
prepared <- prepare_window_data(
  window = windows[1],
  phenotype_data = phenotype_data,
  dosage = dosage,
  covariates = covariates
)
stopifnot(nrow(prepared$X) == length(prepared$samples))
stopifnot(nrow(prepared$Y) == length(prepared$samples))
stopifnot(all(is.finite(prepared$X)), all(is.finite(prepared$Y)))
stopifnot(prepared$covariate_rank >= 1L)
C_model <- cbind(covariates[prepared$samples, , drop = FALSE], intercept = 1)
stopifnot(abs(max(abs(crossprod(C_model, prepared$X)))) < 1e-6)
stopifnot(abs(max(abs(crossprod(C_model, prepared$Y)))) < 1e-6)

bad_covariates <- covariates
rownames(bad_covariates) <- paste0("missing_", seq_len(nrow(bad_covariates)))
no_shared_samples <- tryCatch(
  prepare_window_data(
    window = windows[1],
    phenotype_data = phenotype_data,
    dosage = dosage,
    covariates = bad_covariates
  ),
  error = identity
)
stopifnot(inherits(no_shared_samples, "error"))
stopifnot(grepl("No shared samples", conditionMessage(no_shared_samples), fixed = TRUE))

constant_dosage <- dosage
constant_dosage$X <- cbind(dosage$X, constant = rep(1, nrow(dosage$X)))
constant_dosage$variant_ids <- c(dosage$variant_ids, "chr1:303_G_A")
constant_dosage$metadata <- data.table::rbindlist(list(
  dosage$metadata,
  data.table::data.table(CHROM = "chr1", POS = 303L, REF = "G", ALT = "A")
))
filtered <- prepare_window_data(
  window = windows[1],
  phenotype_data = phenotype_data,
  dosage = constant_dosage,
  covariates = covariates
)
stopifnot(filtered$qc$excluded_variants == 1L)
stopifnot(ncol(filtered$X) == 2L)

source("scripts/trans_window_model.R")

set.seed(1001)
config <- make_model_config()
model_X <- matrix(rnorm(50L * 6L), nrow = 50L, ncol = 6L)
model_Y <- matrix(rnorm(50L * 3L), nrow = 50L, ncol = 3L)
colnames(model_X) <- paste0("variant_", seq_len(ncol(model_X)))
colnames(model_Y) <- paste0("phenotype_", seq_len(ncol(model_Y)))
model_prepared <- list(
  X = model_X,
  Y = scale(model_Y),
  qc = list(window_id = "synthetic_model")
)
prior <- make_canonical_prior(ncol(model_prepared$Y))
stopifnot(inherits(prior, "mash_prior"))

result <- fit_window_mvsusie(model_prepared, config)
stopifnot(isTRUE(result$fit$converged))
stopifnot(identical(result$metadata$residual_variance_mode, "mvsusieR_default"))

pip <- extract_variant_pips(result$fit, model_prepared)
stopifnot(all(c("variant_id", "pip") %in% names(pip)))
credible_sets <- extract_credible_sets(result$fit, model_prepared, config)
stopifnot(all(c("component", "variant_id", "alpha", "pip") %in% names(credible_sets)))
component_effects <- extract_component_effects(result$fit, model_prepared)
stopifnot(all(c("component", "variant_id", "phenotype_id", "posterior_mean") %in% names(component_effects)))

message("Task 1 reader tests passed")
message("Task 2 preprocessing tests passed")
message("Task 3 model tests passed")
