#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

data_dir <- Sys.getenv("IKZF1_DATA_DIR", "/Users/evinmpadhi/Documents/trans sqtls")
bundle_path <- file.path(data_dir, "mvsusie_top25_ikzf1_splicing_bundle.rds")
pip_path <- file.path(data_dir, "IKZF1_each_isoform_susie_variant_pip.tsv.gz")
phenotype_path <- file.path(data_dir, "IKZF1_rsem_transcripts_isopct.tsv")
annotation_path <- file.path(data_dir, "IKZF1_transcript_structure_annotations.tsv")

stopifnot(file.exists(bundle_path), file.exists(pip_path))
stopifnot(file.exists(phenotype_path), file.exists(annotation_path))

bundle <- readRDS(bundle_path)
stopifnot(length(rownames(bundle$X)) > 100L)
stopifnot(ncol(bundle$X) > 0L, ncol(bundle$Y) > 0L)
reference_pip <- fread(pip_path, nrows = 1L)
stopifnot(all(c("variant", "pip") %in% names(reference_pip)))

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (i == length(args)) stop("Missing value for --", key, call. = FALSE)
    out[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

arg <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) {
    if (required) stop("Missing --", name, call. = FALSE)
    return(default)
  }
  value
}

prepare_inputs <- function(output_dir, dosage_path, start, end) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  dosage <- fread(dosage_path, check.names = FALSE)
  required <- c("CHROM", "POS", "REF", "ALT")
  stopifnot(all(required %in% names(dosage)))
  dosage[, POS := as.integer(POS)]
  dosage <- dosage[POS >= start & POS < end]
  if (!nrow(dosage)) stop("No IKZF1 dosage variants in requested interval.", call. = FALSE)
  dosage_path_local <- file.path(output_dir, "dosage.tsv")
  fwrite(dosage, dosage_path_local, sep = "\t")

  raw_phenotypes <- fread(phenotype_path, header = FALSE, check.names = FALSE)
  phenotype_header <- as.character(raw_phenotypes[1L, ])
  phenotypes <- raw_phenotypes[-1L]
  setnames(phenotypes, phenotype_header)
  phenotype_ids <- as.character(phenotypes[[1L]])
  value_matrix <- as.matrix(phenotypes[, -(1:2), with = FALSE])
  keep <- which(vapply(seq_len(nrow(value_matrix)), function(i) {
    values <- suppressWarnings(as.numeric(value_matrix[i, ]))
    sum(is.finite(values)) > 100L && isTRUE(stats::var(values, na.rm = TRUE) > 0)
  }, logical(1L)))
  if (length(keep) < 3L) stop("Fewer than three usable IKZF1 transcript phenotypes.", call. = FALSE)
  keep <- keep[seq_len(min(5L, length(keep)))]
  selected <- phenotypes[keep]
  setnames(selected, 3:ncol(selected), paste0("X", names(selected)[3:ncol(selected)]))
  phenotype_local <- file.path(output_dir, "phenotypes.tsv")
  fwrite(selected, phenotype_local, sep = "\t")

  covariate_path <- Sys.getenv(
    "IKZF1_COVARIATES_FILE",
    "/Users/evinmpadhi/mvsusie_test/data/mvsusie_test_data/COMB.expression.expression_QTL_covariates.INT.tsv"
  )
  raw_covariates <- fread(covariate_path, header = FALSE, check.names = FALSE)
  covariate_header <- as.character(raw_covariates[1L, ])
  covariates <- raw_covariates[-1L]
  setnames(covariates, covariate_header)
  covariates <- covariates[seq_len(min(5L, nrow(covariates)))]
  setnames(covariates, 2:ncol(covariates), paste0("X", names(covariates)[2:ncol(covariates)]))
  covariate_local <- file.path(output_dir, "covariates.tsv")
  fwrite(covariates, covariate_local, sep = "\t")

  windows_local <- file.path(output_dir, "windows.tsv")
  fwrite(data.table(
    window_id = "ikzf1",
    chrom = as.character(dosage$CHROM[[1L]]),
    start = start,
    end = end,
    dosage_file = basename(dosage_path_local)
  ), windows_local, sep = "\t")
  phenotype_manifest_local <- file.path(output_dir, "window_phenotypes.tsv")
  fwrite(data.table(
    window_id = "ikzf1",
    phenotype_id = as.character(selected[[1L]]),
    modality = "isoform_usage",
    phenotype_file = basename(phenotype_local)
  ), phenotype_manifest_local, sep = "\t")
  invisible(list(
    dosage = dosage_path_local,
    phenotype = phenotype_local,
    covariates = covariate_local,
    windows = windows_local,
    window_phenotypes = phenotype_manifest_local,
    phenotype_ids = as.character(selected[[1L]])
  ))
}

check_outputs <- function(output_dir) {
  pip_path_out <- file.path(output_dir, "window", "variant_pip.tsv.gz")
  qc_path_out <- file.path(output_dir, "window", "window_qc.tsv")
  fit_path_out <- file.path(output_dir, "fit.rds")
  stopifnot(file.exists(pip_path_out), file.exists(qc_path_out), file.exists(fit_path_out))
  pip <- fread(pip_path_out)
  qc <- fread(qc_path_out)
  fit_bundle <- readRDS(fit_path_out)
  stopifnot(nrow(pip) > 0L, nrow(qc) == 1L, isTRUE(qc$converged[[1L]]))
  stopifnot(all(pip$variant_id %in% colnames(bundle$X)))
  stopifnot(all(fit_bundle$samples %in% rownames(bundle$X)))
  stopifnot("phenotype_id" %in% names(fit_bundle$phenotype_metadata))
  stopifnot(nrow(fit_bundle$phenotype_metadata) >= 3L)
  stopifnot(all(nzchar(as.character(fit_bundle$phenotype_metadata$phenotype_id))))
  cat("IKZF1 reference and generated-output checks passed\n")
}

args <- parse_args()
if (!is.null(args[["prepare-dir"]])) {
  prepare_inputs(
    output_dir = args[["prepare-dir"]],
    dosage_path = arg(args, "dosage", required = TRUE),
    start = as.integer(arg(args, "start", "50303540")),
    end = as.integer(arg(args, "end", "50391729"))
  )
}

output_dir <- Sys.getenv("IKZF1_OUTPUT_DIR", "")
if (nzchar(output_dir)) check_outputs(output_dir)
