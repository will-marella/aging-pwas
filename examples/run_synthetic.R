#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this example with Rscript.", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) {
  stop("Usage: Rscript examples/run_synthetic.R [OUTPUT_DIR] [N_CORES]",
       call. = FALSE)
}
output_dir <- if (length(args) >= 1L) args[[1L]] else
  file.path(repo_root, "outputs", "synthetic")
n_cores <- if (length(args) >= 2L) suppressWarnings(as.numeric(args[[2L]])) else 1L
if (length(n_cores) != 1L || !is.finite(n_cores) || n_cores < 1 ||
    n_cores > .Machine$integer.max || n_cores != as.integer(n_cores)) {
  stop("N_CORES must be a positive integer.", call. = FALSE)
}
if (file.exists(output_dir) || dir.exists(output_dir)) {
  stop("Output location already exists; choose a new directory.", call. = FALSE)
}

source(file.path(repo_root, "R", "pwas_time.R"))
source(file.path(repo_root, "examples", "synthetic_data.R"))
config <- new.env(parent = baseenv())
config$pwas_time_spec <- pwas_time_spec
sys.source(file.path(repo_root, "specs", "synthetic.R"), envir = config)
synthetic <- make_synthetic_pwas()
stopifnot(identical(synthetic$preprocessing, config$preprocessing))
result <- run_pwas_time(synthetic$pheno, synthetic$omics, config$spec,
                        synthetic$preprocessing, n_cores = as.integer(n_cores))
print(summarize_pwas_time(result))
write_pwas_time(result, output_dir)
saveRDS(synthetic$pheno, file.path(output_dir, "synthetic_pheno.rds"))
saveRDS(synthetic$omics, file.path(output_dir, "synthetic_omics.rds"))
write.csv(synthetic$truth, file.path(output_dir, "synthetic_truth.csv"),
          row.names = FALSE)
message("Synthetic inputs and results written to ",
        normalizePath(output_dir, mustWork = TRUE))
