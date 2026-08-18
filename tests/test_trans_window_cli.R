source("scripts/trans_window_cli.R")

options <- parse_cli_args(
  option_list = list(
    optparse::make_option("--required", type = "character"),
    optparse::make_option("--count", type = "integer", default = 1L),
    optparse::make_option("--max-iter", type = "integer", default = 100L)
  ),
  args = c("--required", "value", "--count", "3", "--max-iter", "25")
)
stopifnot(identical(options$required, "value"))
stopifnot(identical(options$count, 3L))
stopifnot(identical(options$max_iter, 25L))
stopifnot(identical(require_cli_arg(options, "required"), "value"))
stopifnot(identical(optional_cli_arg(options, "missing"), NULL))
message("CLI optparse tests passed")
