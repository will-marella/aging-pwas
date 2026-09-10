# Source this file from any working directory; no analysis runs on source().
.pwas_entry <- sys.frame(1L)$ofile
if (is.null(.pwas_entry)) stop("Use source('/path/to/R/pwas_time.R').", call. = FALSE)
.pwas_source_dir <- dirname(normalizePath(.pwas_entry, mustWork = TRUE))
for (.pwas_file in c("spec.R", "validate.R", "fit.R", "run.R")) {
  sys.source(file.path(.pwas_source_dir, .pwas_file), envir = environment())
}
rm(.pwas_entry, .pwas_file)
