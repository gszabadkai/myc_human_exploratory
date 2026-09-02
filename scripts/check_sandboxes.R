# check_sandboxes.R
# =============================================================================
# NOT PART OF THE NUMBERED PIPELINE. A maintenance check, safe to run any time.
#
# WHY IT EXISTS. Every numbered script ends with an `if (FALSE) { ... }` sandbox
# for interactive inspection. source() never executes those blocks, so they are
# the one part of the repo that no run and no dry run ever exercises. When a
# saved object's column is renamed, the pipeline gets updated and the sandbox
# does not, and the rot stays invisible until someone pastes the block into a
# console weeks later. That is exactly how `null_tests$z` survived being renamed
# to `z_mean_abs` in E11 on 2026-09-02.
#
# HOW IT AVOIDS NEEDING MAINTENANCE ITSELF. It does not know which constants a
# script uses. It reads each script's own top-level `PATH_* <- file.path(...)`
# and `<CONST> <- file.path(DIR_*, ...)` assignments, evaluates only those, and
# then runs the sandbox in that environment. Add a script and it is picked up;
# rename a constant and it follows.
#
# TWO KINDS OF SANDBOX, AND ONLY ONE OF THEM CAN GO STALE. Most begin
# `x <- readRDS(PATH_*)` and are self-contained, so they can be checked here. A
# FIGURES script's sandbox instead re-prints plot objects from its own body
# (E04 prints `f3`), which is a legitimate pattern that simply requires the
# script to have been sourced first. When a sandbox fails on a missing object
# that the script itself assigns at top level, that is reported as NEEDS SOURCE
# rather than as staleness.
#
# A script whose results are not on disk is SKIPPED, not failed.
# Run it after ANY rename of a column in a saved object.
# =============================================================================

source(here::here("scripts", "E00_setup_packages.R"))
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

# Top-level path constants: `NAME <- file.path(...)` at column 0, which is how
# every numbered script declares them.
.path_consts <- function(L) {
  keep <- grepl("^[A-Za-z_][A-Za-z0-9_.]*[[:space:]]*<-[[:space:]]*file\\.path\\(",
                L)
  if (!any(keep)) return(character(0))
  L[keep]
}

.sandbox_code <- function(L) {
  i <- grep("^if \\(FALSE\\) \\{", L)
  if (!length(i)) return(NULL)
  j <- length(L)
  while (j > i[1] && !grepl("^\\}", L[j])) j <- j - 1L
  if (j <= i[1] + 1L) return(NULL)
  paste(L[(i[1] + 1L):(j - 1L)], collapse = "\n")
}

files <- sort(list.files(here::here("scripts"), pattern = "^E[0-9]+.*\\.R$",
                         full.names = TRUE))
message("\nchecking sandbox blocks against results/\n", strrep("=", 66))
n_fail <- 0L; n_skip <- 0L
for (f in files) {
  nm <- tools::file_path_sans_ext(basename(f))
  L  <- readLines(f, warn = FALSE)
  code <- .sandbox_code(L)
  if (is.null(code) || !nzchar(trimws(code))) {
    message(sprintf("  %-38s no sandbox", nm)); next
  }
  env <- new.env(parent = globalenv())
  # The script's own constants, and nothing else from it.
  for (ln in .path_consts(L)) {
    try(eval(parse(text = ln), envir = env), silent = TRUE)
  }
  # If the sandbox's own input is missing, that is a skip rather than a failure.
  paths <- unlist(lapply(ls(env), function(v) {
    val <- get(v, envir = env)
    if (is.character(val) && length(val) == 1L) val else NULL
  }))
  rds <- paths[grepl("\\.rds$", paths)]
  if (length(rds) && !any(file.exists(rds))) {
    message(sprintf("  %-38s SKIPPED - %s not on disk", nm,
                    paste(basename(rds), collapse = ", ")))
    n_skip <- n_skip + 1L; next
  }
  res <- tryCatch({ eval(parse(text = code), envir = env); "OK" },
                  error = function(e) conditionMessage(e))
  if (identical(res, "OK")) {
    message(sprintf("  %-38s OK", nm)); next
  }
  # Is the missing object one the script itself creates? Then the sandbox is
  # meant to be run after sourcing, and this is not rot.
  miss <- regmatches(res, regexpr("(?<=object ')[^']+", res, perl = TRUE))
  own <- length(miss) &&
    any(grepl(paste0("^", miss, "[[:space:]]*(<-|=[^=])"), L))
  if (own) {
    message(sprintf("  %-38s NEEDS SOURCE - re-prints `%s` from the script body",
                    nm, miss))
    n_skip <- n_skip + 1L
  } else {
    n_fail <- n_fail + 1L
    message(sprintf("  %-38s FAILED - %s", nm, res))
  }
}
message(strrep("=", 66))
if (n_fail) {
  message(n_fail, " sandbox block(s) are stale. source() never runs them, so ",
          "nothing\nelse in this repo will ever tell you.")
} else {
  message("no stale sandbox blocks",
          if (n_skip) paste0(" (", n_skip,
                             " skipped or needing the script sourced)") else "")
}
