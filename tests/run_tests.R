#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
source(file.path(repo_root, "R", "pwas_time.R"))
source(file.path(repo_root, "examples", "synthetic_data.R"))
source(file.path(repo_root, "specs", "synthetic.R"))
if (!requireNamespace("testthat", quietly = TRUE)) stop("Install testthat to run tests.")
# Let testthat control diagnostic language within this test-only R process.
Sys.unsetenv("LC_ALL")
testthat::test_dir(file.path(repo_root, "tests", "testthat"), reporter = "summary", stop_on_failure = TRUE)
