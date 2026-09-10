#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this entry point with Rscript.", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: Rscript scripts/run_pwas.R PHENO.rds OMICS.rds CONFIG.R OUTPUT_DIR N_CORES",
       call. = FALSE)
}
pheno_path <- args[[1L]]
omics_path <- args[[2L]]
config_path <- args[[3L]]
output_dir <- args[[4L]]
n_cores <- suppressWarnings(as.numeric(args[[5L]]))
if (length(n_cores) != 1L || !is.finite(n_cores) || n_cores < 1 ||
    n_cores > .Machine$integer.max || n_cores != as.integer(n_cores)) {
  stop("N_CORES must be a positive integer.", call. = FALSE)
}
input_paths <- c(pheno_path, omics_path, config_path)
if (any(!file.exists(input_paths)) || any(dir.exists(input_paths))) {
  stop("PHENO.rds, OMICS.rds, and CONFIG.R must be existing files.", call. = FALSE)
}
if (file.exists(output_dir) || dir.exists(output_dir)) {
  stop("Output location already exists; choose a new directory.", call. = FALSE)
}

source(file.path(repo_root, "R", "pwas_time.R"))
config <- new.env(parent = baseenv())
config$pwas_time_spec <- pwas_time_spec
sys.source(config_path, envir = config)
if (!exists("spec", envir = config, inherits = FALSE) ||
    !exists("preprocessing", envir = config, inherits = FALSE)) {
  stop("CONFIG.R must define spec and preprocessing.", call. = FALSE)
}
pheno <- readRDS(pheno_path)
omics <- readRDS(omics_path)
input_md5 <- tools::md5sum(c(pheno_path, omics_path))
config_md5 <- unname(tools::md5sum(config_path))
if (anyNA(input_md5) || anyNA(config_md5)) {
  stop("Could not compute input or configuration checksums.", call. = FALSE)
}
result <- run_pwas_time(pheno, omics, config$spec, config$preprocessing,
                        n_cores = as.integer(n_cores))
result$metadata$input_md5 <- list(pheno = unname(input_md5[[1L]]),
                                   omics = unname(input_md5[[2L]]))
result$metadata$config_md5 <- config_md5
print(summarize_pwas_time(result))
write_pwas_time(result, output_dir)
message("PWAS results written to ", normalizePath(output_dir, mustWork = TRUE))
