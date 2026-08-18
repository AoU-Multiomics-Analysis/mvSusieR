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

message("Task 1 reader tests passed")
