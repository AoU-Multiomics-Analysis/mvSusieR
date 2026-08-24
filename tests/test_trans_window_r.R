source("scripts/trans_window_io.R")

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

expect_manifest_error(
  phenotype_manifest[modality != "protein"],
  "exactly.*expression.*splicing.*protein"
)
isoform_manifest <- copy(phenotype_manifest)
isoform_manifest$modality[[3L]] <- "isoform_usage"
isoform_manifest$phenotype_id[[3L]] <- "tx_1"
isoform_manifest$outcome_key[[3L]] <- "isoform_usage::tx_1"
expect_manifest_error(
  isoform_manifest,
  "exactly.*expression.*splicing.*protein"
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

# Keep the existing preprocessing tests isolated from the sample-order fix in Task 3.
covariates_by_modality <- lapply(covariates_by_modality, function(matrix) {
  matrix[dosage$sample_ids, , drop = FALSE]
})
covariates <- covariates_by_modality$expression

source("scripts/trans_window_preprocess.R")
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
stopifnot("modality" %in% names(phenotype_data$metadata))
stopifnot("modality" %in% names(prepared$phenotype_metadata))

prepared_modality <- prepare_window_data(
  window = windows[1],
  phenotype_data = phenotype_data,
  dosage = dosage,
  covariates_by_modality = covariates_by_modality
)
genotype_covariates <- make_genotype_covariates(
  lapply(covariates_by_modality, function(matrix) {
    matrix[prepared_modality$samples, , drop = FALSE]
  })
)
genotype_model <- cbind(genotype_covariates, intercept = 1)
stopifnot(abs(max(abs(crossprod(genotype_model, prepared_modality$X)))) < 1e-6)
for (modality in unique(prepared_modality$phenotype_metadata$modality)) {
  phenotype_indices <- which(prepared_modality$phenotype_metadata$modality == modality)
  phenotype_model <- cbind(
    covariates_by_modality[[modality]][prepared_modality$samples, , drop = FALSE],
    intercept = 1
  )
  stopifnot(
    abs(max(abs(crossprod(phenotype_model, prepared_modality$Y[, phenotype_indices, drop = FALSE])))) < 1e-6
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
prepared_conflicting_names <- prepare_window_data(
  window = windows[1],
  phenotype_data = phenotype_data,
  dosage = dosage,
  covariates_by_modality = conflicting_modality_covariates
)
stopifnot(ncol(prepared_conflicting_names$X) == ncol(dosage$X))
stopifnot(all(is.finite(prepared_conflicting_names$X)))
genotype_conflicting_names <- make_genotype_covariates(
  conflicting_modality_covariates
)
stopifnot(all(c(
  "expression::PC1",
  "splicing::PC1",
  "protein::PC1"
) %in% colnames(genotype_conflicting_names)))
stopifnot(identical(
  unname(genotype_conflicting_names[, "expression::PC1"]),
  as.numeric(expression_pc[, "PC1"])
))
stopifnot(identical(
  unname(genotype_conflicting_names[, "splicing::PC1"]),
  as.numeric(splicing_pc[, "PC1"])
))
conflicting_genotype_model <- cbind(
  genotype_conflicting_names[prepared_conflicting_names$samples, , drop = FALSE],
  intercept = 1
)
stopifnot(
  max(abs(crossprod(
    conflicting_genotype_model,
    prepared_conflicting_names$X
  ))) < 1e-6
)

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
stopifnot(identical(result$metadata$residual_variance_mode, "estimated_by_mvsusie"))

greedy_config <- make_model_config(
  L = 4L,
  L_greedy = 2L,
  greedy_lbf_cutoff = 1e6
)
greedy_result <- fit_window_mvsusie(model_prepared, greedy_config)
stopifnot(isTRUE(greedy_result$fit$converged))
stopifnot(nrow(greedy_result$fit$alpha) == 2L)
stopifnot(identical(greedy_result$metadata$config$L_greedy, 2L))
stopifnot(identical(greedy_result$metadata$config$greedy_lbf_cutoff, 1e6))
stopifnot(identical(greedy_result$metadata$L_final, 2L))
stopifnot(isTRUE(greedy_result$metadata$L_greedy_used))

for (invalid_step in list(0, 2.5, 5, Inf, NaN)) {
  invalid_greedy_config <- tryCatch(
    make_model_config(L = 4L, L_greedy = invalid_step),
    error = identity
  )
  stopifnot(inherits(invalid_greedy_config, "error"))
}

for (invalid_cutoff in list(Inf, NaN)) {
  invalid_greedy_config <- tryCatch(
    make_model_config(L = 4L, L_greedy = 2L, greedy_lbf_cutoff = invalid_cutoff),
    error = identity
  )
  stopifnot(inherits(invalid_greedy_config, "error"))
}

pip <- extract_variant_pips(result$fit, model_prepared)
stopifnot(all(c("variant_id", "pip") %in% names(pip)))
credible_sets <- extract_credible_sets(result$fit, model_prepared, config)
stopifnot(all(c("component", "variant_id", "alpha", "pip") %in% names(credible_sets)))
component_effects <- extract_component_effects(result$fit, model_prepared)
stopifnot(all(c("component", "variant_id", "phenotype_id", "posterior_mean") %in% names(component_effects)))

message("Task 1 reader tests passed")
message("Task 2 preprocessing tests passed")
message("Task 3 model tests passed")
