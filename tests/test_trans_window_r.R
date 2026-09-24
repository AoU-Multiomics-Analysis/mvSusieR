source("scripts/trans_window_io.R")

expect_error_matching <- function(expression, pattern) {
  observed <- tryCatch(
    {
      force(expression)
      NA_character_
    },
    error = function(condition) conditionMessage(condition)
  )
  stopifnot(!is.na(observed), grepl(pattern, observed, ignore.case = TRUE))
}

fixture_dir <- commandArgs(trailingOnly = TRUE)[[1L]]
stopifnot(dir.exists(fixture_dir))
fixture <- function(name) file.path(fixture_dir, name)

windows <- read_windows_manifest(fixture("windows.tsv"))
stopifnot(nrow(windows) == 1L, windows$window_id == "w1")

phenotype_manifest <- read_window_phenotypes_manifest(
  fixture("window_phenotypes.tsv")
)
stopifnot(nrow(phenotype_manifest) == 3L)
stopifnot(identical(
  phenotype_manifest$outcome_key,
  c(
    "expression::ENSG000001.1",
    "splicing::splice_1",
    "protein::prot_1"
  )
))

dosage <- read_wide_dosage(fixture("window_1_dosage.tsv"))
stopifnot(identical(dim(dosage$X), c(6L, 2L)))
stopifnot(identical(dosage$sample_ids, as.character(1:6)))
stopifnot(identical(dosage$variant_ids, c("chr1:101_A_G", "chr1:202_C_T")))

phenotypes <- read_phenotype_rows(
  fixture("expression.tsv"),
  "expression",
  "expression::ENSG000001.1"
)
stopifnot(ncol(phenotypes$Y) == 1L, nrow(phenotypes$Y) == 6L)
stopifnot(identical(colnames(phenotypes$Y), "expression::ENSG000001.1"))

numeric_header_covariates_path <- tempfile(fileext = ".tsv")
writeLines(
  c(
    "ID\t1001\t1002\t1003",
    "PC1\t1\t2\t3",
    "PC2\t4\t5\t6"
  ),
  numeric_header_covariates_path
)
numeric_header_covariates <- read_covariate_file(
  numeric_header_covariates_path
)
stopifnot(identical(rownames(numeric_header_covariates), c("1001", "1002", "1003")))
stopifnot(identical(colnames(numeric_header_covariates), c("PC1", "PC2")))

covariates_by_modality <- read_joint_covariates(
  expression_path = fixture("expression_covariates.tsv"),
  splicing_path = fixture("splicing_covariates.tsv"),
  protein_path = fixture("protein_covariates.tsv")
)
stopifnot(identical(
  names(covariates_by_modality),
  c("expression", "splicing", "protein")
))
stopifnot(all(vapply(covariates_by_modality, ncol, integer(1L)) == 2L))
stopifnot(!identical(
  covariates_by_modality$expression[, "PC1"],
  covariates_by_modality$splicing[, "PC1"]
))
covariates <- covariates_by_modality$expression

phenotype_data <- read_window_phenotypes(
  window_id = "w1",
  phenotype_manifest = phenotype_manifest,
  phenotype_files = c(
    fixture("expression.tsv"),
    fixture("splicing.tsv"),
    fixture("protein.tsv")
  )
)
stopifnot(identical(
  sort(unique(phenotype_data$modalities)),
  c("expression", "protein", "splicing")
))
stopifnot(identical(
  phenotype_data$phenotype_ids,
  phenotype_manifest$outcome_key
))
stopifnot(any(grepl("^protein::", phenotype_data$phenotype_ids)))
stopifnot(identical(
  phenotype_data$metadata$phenotype_id,
  phenotype_manifest$phenotype_id
))

expect_manifest_error <- function(manifest, pattern) {
  path <- tempfile(fileext = ".tsv")
  fwrite(manifest, path, sep = "\t")
  observed <- tryCatch(
    {
      read_window_phenotypes_manifest(path)
      NA_character_
    },
    error = function(condition) conditionMessage(condition)
  )
  stopifnot(!is.na(observed), grepl(pattern, observed, ignore.case = TRUE))
}

protein_missing_manifest <- phenotype_manifest[modality != "protein"]
protein_missing_path <- tempfile(fileext = ".tsv")
fwrite(protein_missing_manifest, protein_missing_path, sep = "\t")
protein_missing_read <- read_window_phenotypes_manifest(protein_missing_path)
stopifnot(identical(
  sort(unique(protein_missing_read$modality)),
  c("expression", "splicing")
))
expect_manifest_error(phenotype_manifest[0L], "at least one outcome")
isoform_manifest <- copy(phenotype_manifest)
isoform_manifest$modality[[3L]] <- "isoform_usage"
isoform_manifest$phenotype_id[[3L]] <- "tx_1"
isoform_manifest$outcome_key[[3L]] <- "isoform_usage::tx_1"
expect_manifest_error(
  isoform_manifest,
  "unsupported.*modality"
)
duplicate_outcome_manifest <- rbindlist(list(
  phenotype_manifest,
  phenotype_manifest[1L]
))
expect_manifest_error(duplicate_outcome_manifest, "duplicate.*outcome")

duplicate_sample_covariates <- tempfile(fileext = ".tsv")
writeLines(
  c("covariate\tX1\t1", "PC1\t1\t2"),
  duplicate_sample_covariates
)
duplicate_sample_error <- tryCatch(
  read_covariate_file(duplicate_sample_covariates),
  error = identity
)
stopifnot(inherits(duplicate_sample_error, "error"))
stopifnot(grepl(
  "duplicate sample IDs",
  conditionMessage(duplicate_sample_error),
  fixed = TRUE
))

no_overlap_covariates <- tempfile(fileext = ".tsv")
writeLines(
  c("covariate\tZ1\tZ2", "PC1\t1\t2"),
  no_overlap_covariates
)
no_overlap_error <- tryCatch(
  read_joint_covariates(
    fixture("expression_covariates.tsv"),
    fixture("splicing_covariates.tsv"),
    no_overlap_covariates
  ),
  error = identity
)
stopifnot(inherits(no_overlap_error, "error"))
stopifnot(grepl(
  "no shared sample IDs",
  conditionMessage(no_overlap_error),
  fixed = TRUE
))

source("scripts/trans_window_logging.R")
source("scripts/trans_window_preprocess.R")

raw_covariates_by_modality <- covariates_by_modality
genotype_union <- make_genotype_covariates(
  raw_covariates_by_modality,
  samples = dosage$sample_ids
)
stopifnot(identical(
  colnames(genotype_union$matrix),
  c("expression::PC1", "SHARED", "splicing::PC1", "protein::PC1")
))
stopifnot(nrow(genotype_union$provenance) == 6L)
stopifnot(sum(genotype_union$provenance$final_name == "SHARED") == 3L)
stopifnot(length(unique(
  genotype_union$provenance$deduplicated_to[
    genotype_union$provenance$original_name == "SHARED"
  ]
)) == 1L)

incomplete_covariates <- lapply(raw_covariates_by_modality, function(matrix) {
  matrix[, "PC1", drop = FALSE]
})
incomplete_covariates$expression["1", "PC1"] <- NA_real_

phenotype_with_constant <- phenotype_data
constant_key <- "expression::constant"
phenotype_with_constant$Y <- cbind(
  phenotype_data$Y[, 1L, drop = FALSE],
  rep(1, nrow(phenotype_data$Y)),
  phenotype_data$Y[, 2:3, drop = FALSE]
)
colnames(phenotype_with_constant$Y) <- c(
  phenotype_data$phenotype_ids[[1L]],
  constant_key,
  phenotype_data$phenotype_ids[2:3]
)
phenotype_with_constant$Y["1", 1L] <- NA_real_
constant_metadata <- data.table::copy(phenotype_data$metadata[1L])
constant_metadata[, `:=`(
  outcome_key = constant_key,
  phenotype_id = "constant",
  modality = "expression"
)]
phenotype_with_constant$metadata <- rbindlist(list(
  phenotype_data$metadata[1L],
  constant_metadata,
  phenotype_data$metadata[2:3]
))
phenotype_with_constant$phenotype_ids <- colnames(phenotype_with_constant$Y)
phenotype_with_constant$modalities <- c(
  "expression", "expression", "splicing", "protein"
)

prepared_joint <- prepare_joint_window_data(
  window = windows[1L],
  phenotype_data = phenotype_with_constant,
  dosage = dosage,
  covariates_by_modality = incomplete_covariates
)
stopifnot(identical(prepared_joint$samples, as.character(2:6)))
stopifnot(identical(
  as.character(prepared_joint$phenotype_metadata$modality),
  c("expression", "splicing", "protein")
))
stopifnot(identical(
  colnames(prepared_joint$Y),
  prepared_joint$phenotype_metadata$outcome_key
))
stopifnot(max(abs(apply(prepared_joint$Y, 2L, sd) - 1)) < 1e-10)
stopifnot(!all(abs(apply(prepared_joint$X, 2L, sd) - 1) < 1e-8))

for (modality in required_joint_modalities()) {
  indices <- which(prepared_joint$phenotype_metadata$modality == modality)
  model <- cbind(
    incomplete_covariates[[modality]][prepared_joint$samples, , drop = FALSE],
    intercept = 1
  )
  stopifnot(max(abs(crossprod(
    model,
    prepared_joint$Y[, indices, drop = FALSE]
  ))) < 1e-6)
}
genotype_model <- cbind(prepared_joint$genotype_covariates, intercept = 1)
stopifnot(max(abs(crossprod(genotype_model, prepared_joint$X))) < 1e-6)
stopifnot(identical(
  colnames(prepared_joint$genotype_covariates),
  c("expression::PC1", "splicing::PC1", "protein::PC1")
))
stopifnot(nrow(prepared_joint$covariate_provenance) == 3L)

prepared <- prepare_joint_window_data(
  window = windows[1],
  phenotype_data = phenotype_data,
  dosage = dosage,
  covariates_by_modality = raw_covariates_by_modality
)
stopifnot(inherits(validate_prepared_window(prepared, "w1"), "list"))
expect_error_matching(
  validate_prepared_window(list(X = matrix(1)), "w1"),
  "missing required fields"
)
expect_error_matching(
  validate_prepared_window(prepared, "wrong_window"),
  "does not match requested window"
)
stopifnot(nrow(prepared$X) == length(prepared$samples))
stopifnot(nrow(prepared$Y) == length(prepared$samples))
stopifnot(all(is.finite(prepared$X)), all(is.finite(prepared$Y)))
stopifnot(prepared$covariate_rank >= 1L)
stopifnot("modality" %in% names(phenotype_data$metadata))
stopifnot("modality" %in% names(prepared$phenotype_metadata))

for (missing_modality in required_joint_modalities()) {
  keep <- phenotype_data$modalities != missing_modality
  subset_data <- phenotype_data
  subset_data$Y <- phenotype_data$Y[, keep, drop = FALSE]
  subset_data$metadata <- data.table::copy(phenotype_data$metadata[keep])
  subset_data$phenotype_ids <- phenotype_data$phenotype_ids[keep]
  subset_data$modalities <- phenotype_data$modalities[keep]
  subset_prepared <- prepare_joint_window_data(
    window = windows[1],
    phenotype_data = subset_data,
    dosage = dosage,
    covariates_by_modality = raw_covariates_by_modality
  )
  stopifnot(!missing_modality %in% subset_prepared$phenotype_metadata$modality)
  stopifnot(identical(
    subset_prepared$phenotype_covariate_rank[[missing_modality]],
    0L
  ))
  stopifnot(identical(
    colnames(subset_prepared$genotype_covariates),
    c("expression::PC1", "SHARED", "splicing::PC1", "protein::PC1")
  ))
}

genotype_model <- cbind(prepared$genotype_covariates, intercept = 1)
stopifnot(abs(max(abs(crossprod(genotype_model, prepared$X)))) < 1e-6)
for (modality in unique(prepared$phenotype_metadata$modality)) {
  phenotype_indices <- which(prepared$phenotype_metadata$modality == modality)
  phenotype_model <- cbind(
    raw_covariates_by_modality[[modality]][prepared$samples, , drop = FALSE],
    intercept = 1
  )
  stopifnot(
    abs(max(abs(crossprod(
      phenotype_model,
      prepared$Y[, phenotype_indices, drop = FALSE]
    )))) < 1e-6
  )
}

expression_pc <- matrix(seq_len(6L), ncol = 1L,
                        dimnames = list(as.character(1:6), "PC1"))
splicing_pc <- matrix(rep(c(-1, 1), 3L), ncol = 1L,
                      dimnames = list(as.character(1:6), "PC1"))
conflicting_modality_covariates <- list(
  expression = expression_pc,
  splicing = splicing_pc,
  protein = expression_pc
)
prepared_conflicting_names <- prepare_joint_window_data(
  window = windows[1],
  phenotype_data = phenotype_data,
  dosage = dosage,
  covariates_by_modality = conflicting_modality_covariates
)
stopifnot(ncol(prepared_conflicting_names$X) == ncol(dosage$X))
stopifnot(all(is.finite(prepared_conflicting_names$X)))
genotype_conflicting_names <- make_genotype_covariates(
  conflicting_modality_covariates,
  samples = dosage$sample_ids
)
stopifnot(all(c(
  "expression::PC1",
  "splicing::PC1"
) %in% colnames(genotype_conflicting_names$matrix)))
stopifnot(!"protein::PC1" %in% colnames(genotype_conflicting_names$matrix))
stopifnot(identical(
  unname(genotype_conflicting_names$matrix[, "expression::PC1"]),
  as.numeric(expression_pc[, "PC1"])
))
stopifnot(identical(
  unname(genotype_conflicting_names$matrix[, "splicing::PC1"]),
  as.numeric(splicing_pc[, "PC1"])
))
conflicting_genotype_model <- cbind(
  genotype_conflicting_names$matrix[
    prepared_conflicting_names$samples,
    ,
    drop = FALSE
  ],
  intercept = 1
)
stopifnot(
  max(abs(crossprod(
    conflicting_genotype_model,
    prepared_conflicting_names$X
  ))) < 1e-6
)

bad_covariates <- lapply(raw_covariates_by_modality, identity)
rownames(bad_covariates$protein) <- paste0(
  "missing_",
  seq_len(nrow(bad_covariates$protein))
)
no_shared_samples <- tryCatch(
  prepare_joint_window_data(
    window = windows[1],
    phenotype_data = phenotype_data,
    dosage = dosage,
    covariates_by_modality = bad_covariates
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
filtered <- prepare_joint_window_data(
  window = windows[1],
  phenotype_data = phenotype_data,
  dosage = constant_dosage,
  covariates_by_modality = raw_covariates_by_modality
)
stopifnot(filtered$qc$excluded_variants == 1L)
stopifnot(ncol(filtered$X) == 2L)

source("scripts/trans_window_model.R")

set.seed(1001)
config <- make_model_config()
stopifnot(
  identical(config$start_L, 10L),
  identical(config$step_L, 5L),
  identical(config$max_L, 40L),
  identical(config$greedy_lbf_cutoff, 1),
  identical(config$mashr_n_pca, 5L),
  identical(config$mashr_strong_lfsr, 0.05),
  identical(config$coverage, 0.95),
  identical(config$min_abs_corr, 0.5)
)
stopifnot(!any(c(
  "prior_method",
  "mashr_use_ed",
  "estimate_residual_variance",
  "marginal_output"
) %in% names(config)))

update_test_Y <- matrix(rnorm(60L), nrow = 20L, ncol = 3L)
raw_update_prior <- mvsusieR::create_mixture_prior(R = 3L, null_weight = 0)
raw_update_prior$xUlist <- lapply(raw_update_prior$xUlist, `*`, 7)
update_spec <- make_mvsusie_update_spec(raw_update_prior, update_test_Y)
stopifnot(
  identical(update_spec$prior_variance, raw_update_prior),
  isTRUE(all.equal(
    update_spec$residual_variance,
    stats::cov(update_test_Y),
    tolerance = 1e-12
  )),
  isTRUE(update_spec$estimate_residual_variance),
  isTRUE(update_spec$estimate_prior_variance),
  isTRUE(update_spec$estimate_prior_mixture_weights)
)
converted_update_prior <- raw_update_prior
attr(converted_update_prior, "mvsusie_outcome_se_scale") <- rep(1, 3L)
converted_prior_error <- tryCatch(
  make_mvsusie_update_spec(converted_update_prior, update_test_Y),
  error = identity
)
stopifnot(
  inherits(converted_prior_error, "error"),
  grepl("raw mashr prior", conditionMessage(converted_prior_error), fixed = TRUE)
)

model_X <- matrix(rnorm(50L * 20L), nrow = 50L, ncol = 20L)
model_Y <- matrix(rnorm(50L * 3L), nrow = 50L, ncol = 3L)
recorded_L <- integer()
recorded_init <- list()
returned_fits <- list()
stub_fit <- function(X, Y, L, model_init, ...) {
  recorded_L <<- c(recorded_L, L)
  recorded_init[length(recorded_init) + 1L] <<- list(model_init)
  lbf <- if (L == 10L) rep(2, 10L) else c(rep(2, L - 1L), 0.8)
  fit <- list(
    alpha = matrix(1 / ncol(X), nrow = L, ncol = ncol(X)),
    lbf = stats::setNames(lbf, paste0("L", seq_len(L))),
    niter = L + 1L,
    converged = TRUE,
    sets = list(cs = list(L1 = 1L)),
    single_effect_lfsr = matrix(0.1, nrow = L, ncol = ncol(Y))
  )
  class(fit) <- "mvsusie"
  returned_fits[[length(returned_fits) + 1L]] <<- fit
  fit
}

scheduled <- fit_mvsusie_greedy_schedule(
  X = model_X,
  Y = model_Y,
  prior = mvsusieR::create_mixture_prior(R = 3L, null_weight = 0),
  start_L = 10L,
  step_L = 5L,
  max_L = 40L,
  greedy_lbf_cutoff = 1,
  fit_fun = stub_fit
)
stopifnot(identical(recorded_L, c(10L, 15L)))
stopifnot(is.null(recorded_init[[1L]]))
stopifnot(identical(recorded_init[[2L]], returned_fits[[1L]]))
stopifnot(identical(scheduled$fit, returned_fits[[2L]]))
stopifnot(identical(names(scheduled), c("fit", "history")))
stopifnot(identical(scheduled$history$action, c("continue", "saturated")))
stopifnot(identical(scheduled$history$requested_L, c(10L, 15L)))
stopifnot(identical(scheduled$history$minimum_lbf, c(2, 0.8)))

for (invalid_config in list(
  list(start_L = 0L),
  list(step_L = 0L),
  list(step_L = 2.5),
  list(start_L = 20L, max_L = 10L),
  list(start_L = 10L, step_L = 7L, max_L = 40L),
  list(greedy_lbf_cutoff = Inf)
)) {
  invalid_result <- tryCatch(
    do.call(make_model_config, invalid_config),
    error = identity
  )
  stopifnot(inherits(invalid_result, "error"))
}

nonconverged_error <- tryCatch(
  fit_mvsusie_greedy_schedule(
    X = model_X,
    Y = model_Y,
    prior = mvsusieR::create_mixture_prior(R = 3L, null_weight = 0),
    fit_fun = function(...) {
      list(
        alpha = matrix(0.05, 10L, 20L),
        lbf = rep(2, 10L),
        niter = 1L,
        converged = FALSE
      )
    }
  ),
  error = identity
)
stopifnot(inherits(nonconverged_error, "error"))
stopifnot(grepl("did not converge", conditionMessage(nonconverged_error)))

nonfinite_error <- tryCatch(
  fit_mvsusie_greedy_schedule(
    X = model_X,
    Y = model_Y,
    prior = mvsusieR::create_mixture_prior(R = 3L, null_weight = 0),
    fit_fun = function(...) {
      list(
        alpha = matrix(NA_real_, 10L, 20L),
        lbf = rep(2, 10L),
        niter = 1L,
        converged = TRUE
      )
    }
  ),
  error = identity
)
stopifnot(inherits(nonfinite_error, "error"))
stopifnot(grepl("non-finite", conditionMessage(nonfinite_error)))

maximum_schedule <- fit_mvsusie_greedy_schedule(
  X = model_X,
  Y = model_Y,
  prior = mvsusieR::create_mixture_prior(R = 3L, null_weight = 0),
  start_L = 10L,
  step_L = 5L,
  max_L = 10L,
  fit_fun = function(X, Y, L, model_init, ...) {
    list(
      alpha = matrix(0.05, L, ncol(X)),
      lbf = rep(2, L),
      niter = 2L,
      converged = TRUE,
      sets = list(cs = NULL),
      single_effect_lfsr = matrix(0.1, L, ncol(Y))
    )
  }
)
stopifnot(identical(maximum_schedule$history$action, "maximum"))

bad_prior <- mvsusieR::create_mixture_prior(R = 2L, null_weight = 0)
prior_dimension_error <- tryCatch(
  validate_prior_for_outcomes(bad_prior, 3L),
    error = identity
  )
stopifnot(inherits(prior_dimension_error, "error"))

message("Task 1 reader tests passed")
message("Task 2 preprocessing tests passed")
message("Task 3 joint model configuration tests passed")
