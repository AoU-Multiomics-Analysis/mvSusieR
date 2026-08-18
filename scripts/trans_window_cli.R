parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args)) return(list())
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected command-line token: ", key, call. = FALSE)
    }
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

require_cli_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) {
    stop("Missing required command-line argument: --", name, call. = FALSE)
  }
  value
}

optional_cli_arg <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (is.null(value) || identical(value, TRUE)) default else value
}

split_cli_paths <- function(value) {
  if (is.null(value) || !nzchar(value)) return(character())
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

as_cli_integer <- function(args, name, default) {
  value <- optional_cli_arg(args, name, as.character(default))
  out <- suppressWarnings(as.integer(value))
  if (is.na(out)) stop("--", name, " must be an integer.", call. = FALSE)
  out
}

as_cli_numeric <- function(args, name, default) {
  value <- optional_cli_arg(args, name, as.character(default))
  out <- suppressWarnings(as.numeric(value))
  if (is.na(out)) stop("--", name, " must be numeric.", call. = FALSE)
  out
}

save_rds_checked <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path)
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Failed to write RDS output: ", path, call. = FALSE)
  }
}
