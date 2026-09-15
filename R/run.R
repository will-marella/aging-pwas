.pwas_object_md5 <- function(object) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 3)
  unname(tools::md5sum(path))
}

.pwas_code_metadata <- function() {
  root <- dirname(.pwas_source_dir)
  files <- sort(list.files(.pwas_source_dir, pattern = "\\.R$", full.names = TRUE))
  hashes <- tools::md5sum(files)
  names(hashes) <- paste0("R/", basename(files))
  git <- if (file.exists(file.path(root, ".git")) && nzchar(Sys.which("git"))) {
    read_git <- function(args) {
      result <- suppressWarnings(system2("git", c("-C", shQuote(root), args), stdout = TRUE, stderr = FALSE))
      if (is.null(attr(result, "status"))) result else NA_character_
    }
    list(commit = read_git(c("rev-parse", "HEAD")), status = read_git(c("status", "--porcelain")))
  } else list(commit = NA_character_, status = "not_a_git_checkout")
  list(source_md5 = hashes, git = git)
}

.pwas_adjust <- function(table, groups, component) {
  table$BH_P_VALUE <- rep(NA_real_, nrow(table))
  families <- data.frame(COMPONENT = character(), FAMILY = character(),
    N_REQUESTED = integer(), N_VALID_P = integer(), N_WITHHELD = integer(), stringsAsFactors = FALSE)
  if (!nrow(table)) return(list(table = table, families = families))
  keys <- do.call(paste, c(table[groups], sep = " / "))
  for (key in unique(keys)) {
    index <- which(keys == key)
    valid <- index[table$INFERENCE_OK[index] & is.finite(table$P_VALUE[index])]
    table$BH_P_VALUE[valid] <- stats::p.adjust(table$P_VALUE[valid], method = "BH")
    families <- rbind(families, data.frame(COMPONENT = component, FAMILY = key,
      N_REQUESTED = length(index), N_VALID_P = length(valid), N_WITHHELD = length(index) - length(valid)))
  }
  list(table = table, families = families)
}

run_pwas_time <- function(pheno, omics, spec, preprocessing, n_cores = 1L, verbose = TRUE) {
  for (package in c("lme4", "lmerTest")) if (!requireNamespace(package, quietly = TRUE))
    stop("Install required R package: ", package, call. = FALSE)
  if (!is.numeric(n_cores) || length(n_cores) != 1L || is.na(n_cores) || !is.finite(n_cores) ||
      n_cores < 1 || n_cores != floor(n_cores) || n_cores > .Machine$integer.max)
    stop("n_cores must be an explicit positive integer.", call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("verbose must be TRUE or FALSE.", call. = FALSE)
  if (.Platform$OS.type != "unix" && n_cores > 1L)
    stop("Parallel execution requires Unix; use n_cores = 1 on Windows.", call. = FALSE)
  started <- Sys.time()
  code_metadata <- .pwas_code_metadata()
  spec <- .validate_pwas_spec(spec)
  validated <- validate_pwas_inputs(pheno, omics, spec, preprocessing)
  template <- validated$pheno[, c("TIME_YEARS", "AGE_C", "FEMALE", names(spec$covariates)), drop = FALSE]
  template$RESPONSE <- 0
  design_names <- colnames(stats::model.matrix(.pwas_formula(spec, random = FALSE), template))
  if (anyDuplicated(design_names))
    stop("Declared factor coding produces duplicate model coefficient names; rename covariates or levels.", call. = FALSE)
  rm(template)
  indices <- seq_len(nrow(validated$omics))
  n_workers <- min(as.integer(n_cores), length(indices))
  if (verbose) message("[PWAS_Time] Fitting ", length(indices), " analytes with ", n_workers, " worker(s).")
  worker <- function(i) .pwas_analyte(i, validated, spec)
  results <- if (n_workers == 1L) lapply(indices, worker) else
    parallel::mclapply(indices, worker, mc.cores = n_workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
  if (any(vapply(results, inherits, logical(1), what = "try-error")))
    stop("A parallel worker failed outside a model fit. No partial result was returned.", call. = FALSE)
  components <- c("coefficients", "model_qc", "visit_coverage", "exclusions")
  result <- stats::setNames(lapply(components, function(name) {
    table <- do.call(rbind, lapply(results, `[[`, name)); rownames(table) <- NULL; table
  }), components)
  result$covariance <- stats::setNames(lapply(results, `[[`, "covariance"), rownames(validated$omics))
  adjusted <- .pwas_adjust(result$coefficients, "TERM", "coefficients")
  result$coefficients <- adjusted$table
  result$multiplicity <- adjusted$families
  finished <- Sys.time()
  result$metadata <- list(pipeline_version = "0.2.1", specification = spec,
    specification_md5 = .pwas_object_md5(spec), preprocessing = preprocessing,
    input_qc = validated$qc, code = code_metadata,
    fixed_formula = paste(deparse(.pwas_formula(spec, random = FALSE)), collapse = " "),
    full_formula = paste(deparse(.pwas_formula(spec)), collapse = " "),
    confidence_level = spec$confidence_level,
    estimation = "ML; one model per eligible analyte",
    inference = list(coefficients = "Satterthwaite t"),
    multiplicity = "BH across valid proteins separately within each coefficient term; includes flagged singular fits; withheld p-values excluded",
    started_utc = format(started, tz = "UTC", usetz = TRUE),
    finished_utc = format(finished, tz = "UTC", usetz = TRUE),
    elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
    requested_cores = as.integer(n_cores), dispatched_workers = n_workers,
    observed_worker_pids = sort(unique(result$model_qc$WORKER_PID)),
    package_versions = vapply(c("lme4", "lmerTest", "Matrix"), function(p) as.character(utils::packageVersion(p)), character(1)),
    session_info = capture.output(utils::sessionInfo()))
  class(result) <- c("pwas_time_result", "list")
  if (verbose) message("[PWAS_Time] Finished in ", round(result$metadata$elapsed_seconds, 1), " seconds; ",
    sum(vapply(results, function(x) all(x$coefficients$INFERENCE_OK), logical(1L))),
    "/", nrow(result$model_qc), " analytes have complete inference.")
  result
}

summarize_pwas_time <- function(result) {
  if (!inherits(result, "pwas_time_result")) stop("Expected a pwas_time_result.", call. = FALSE)
  status <- as.data.frame(table(result$model_qc$STATUS), stringsAsFactors = FALSE)
  names(status) <- c("STATUS", "N_ANALYTES")
  list(input = result$metadata$input_qc, model_status = status,
       exclusions = stats::aggregate(N_SAMPLES ~ REASON, result$exclusions, sum),
       multiplicity = result$multiplicity,
       note = "Exclusions are summed across analytes; they are not unique participant/sample counts.")
}

.pwas_results_table <- function(result) {
  table <- result$model_qc[, c("ANALYTE_NAME", "STATUS", "N_OBS", "N_SUBJECTS",
                              "FULL_CONVERGED", "FULL_SINGULAR"), drop = FALSE]
  fields <- c("ESTIMATE", "SE", "DF", "CI_LOW", "CI_HIGH", "P_VALUE",
              "BH_P_VALUE", "INFERENCE_OK")
  for (term in unique(result$coefficients$TERM)) {
    rows <- result$coefficients[result$coefficients$TERM == term, , drop = FALSE]
    rows <- rows[match(table$ANALYTE_NAME, rows$ANALYTE_NAME), , drop = FALSE]
    for (field in fields) table[[paste0(term, "__", field)]] <- rows[[field]]
  }
  table
}

write_pwas_time <- function(result, output_dir) {
  if (!inherits(result, "pwas_time_result")) stop("Expected a pwas_time_result.", call. = FALSE)
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir))
    stop("output_dir must be one nonempty path.", call. = FALSE)
  if (file.exists(output_dir) && !dir.exists(output_dir)) stop("output_dir is a file.", call. = FALSE)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
    stop("output_dir is not empty; choose a new directory to preserve existing results.", call. = FALSE)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) stop("Cannot create output_dir.", call. = FALSE)
  path <- file.path(output_dir, "result.rds")
  temp <- tempfile(pattern = ".result-", tmpdir = output_dir)
  on.exit(unlink(temp), add = TRUE)
  saveRDS(result, temp)
  if (!file.rename(temp, path)) stop("Could not finalize result.rds.", call. = FALSE)
  paths <- c(result = path, results = file.path(output_dir, "results.csv"))
  utils::write.csv(.pwas_results_table(result), paths["results"], row.names = FALSE, na = "")
  invisible(paths)
}
