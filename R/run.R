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

run_pwas_time <- function(pheno, omics, spec, preprocessing, n_cores = 1L, verbose = TRUE,
                          output_dir = NULL, checkpoint_every = 50L) {
  for (package in c("lme4", "lmerTest")) if (!requireNamespace(package, quietly = TRUE))
    stop("Install required R package: ", package, call. = FALSE)
  if (!is.numeric(n_cores) || length(n_cores) != 1L || is.na(n_cores) || !is.finite(n_cores) ||
      n_cores < 1 || n_cores != floor(n_cores) || n_cores > .Machine$integer.max)
    stop("n_cores must be an explicit positive integer.", call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("verbose must be TRUE or FALSE.", call. = FALSE)
  if (.Platform$OS.type != "unix" && n_cores > 1L)
    stop("Parallel execution requires Unix; use n_cores = 1 on Windows.", call. = FALSE)
  if (!is.null(output_dir)) {
    if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir))
      stop("output_dir must be one nonempty path.", call. = FALSE)
    if (file.exists(output_dir) && !dir.exists(output_dir)) stop("output_dir is a file.", call. = FALSE)
    if (!is.numeric(checkpoint_every) || length(checkpoint_every) != 1L ||
        !is.finite(checkpoint_every) || checkpoint_every < 1 ||
        checkpoint_every != floor(checkpoint_every) || checkpoint_every > .Machine$integer.max)
      stop("checkpoint_every must be a positive integer.", call. = FALSE)
  }
  started <- Sys.time()
  code_metadata <- .pwas_code_metadata()
  package_versions <- vapply(c("lme4", "lmerTest", "Matrix"),
    function(p) as.character(utils::packageVersion(p)), character(1))
  spec <- .validate_pwas_spec(spec)
  validated <- validate_pwas_inputs(pheno, omics, spec, preprocessing)
  template <- validated$pheno[, c("TIME_YEARS", "AGE_C", "FEMALE", names(spec$covariates)), drop = FALSE]
  template$RESPONSE <- 0
  design_names <- colnames(stats::model.matrix(.pwas_formula(spec, random = FALSE), template))
  if (anyDuplicated(design_names))
    stop("Declared factor coding produces duplicate model coefficient names; rename covariates or levels.", call. = FALSE)
  rm(template)
  n_analytes <- nrow(validated$omics)
  components <- c("coefficients", "model_qc", "visit_coverage", "exclusions")
  result <- stats::setNames(vector("list", length(components)), components)
  result$covariance <- list()
  signature <- NULL
  n_completed <- 0L
  previous_elapsed <- 0
  started_utc <- format(started, tz = "UTC", usetz = TRUE)
  if (!is.null(output_dir)) {
    if (verbose) message("[PWAS_Time] Checking checkpoint inputs.")
    signature <- .pwas_object_md5(list(pheno = validated$pheno, omics = validated$omics,
      spec = spec, preprocessing = preprocessing, code = code_metadata$source_md5,
      package_versions = package_versions))
    checkpoint_path <- file.path(output_dir, "result.rds")
    if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
      if (!file.exists(checkpoint_path))
        stop("Output directory is not empty and has no result.rds checkpoint; choose a new directory.", call. = FALSE)
      saved <- readRDS(checkpoint_path)
      if (!inherits(saved, "pwas_time_result") || is.null(saved$metadata$checkpoint))
        stop("Existing result is not a resumable checkpoint; choose a new directory.", call. = FALSE)
      if (!identical(saved$metadata$checkpoint$signature, signature))
        stop("Checkpoint inputs, model, preprocessing, code, or package versions differ; choose a new directory.", call. = FALSE)
      n_completed <- nrow(saved$model_qc)
      if (n_completed < 1L || n_completed > n_analytes ||
          !identical(saved$model_qc$ANALYTE_NAME, rownames(validated$omics)[seq_len(n_completed)]) ||
          !identical(saved$metadata$checkpoint$complete, n_completed == n_analytes))
        stop("Checkpoint analyte inventory is inconsistent.", call. = FALSE)
      result <- saved
      previous_elapsed <- result$metadata$elapsed_seconds
      started_utc <- result$metadata$started_utc
      if (n_completed == n_analytes) {
        .pwas_write_files(result, output_dir)
        if (verbose) message("[PWAS_Time] All ", n_analytes, " analytes are already complete; no models refitted.")
        return(result)
      }
      if (verbose) message("[PWAS_Time] Resuming after ", n_completed, "/", n_analytes, " completed analytes.")
    }
  }
  indices <- seq.int(n_completed + 1L, n_analytes)
  batch_size <- if (is.null(output_dir)) length(indices) else as.integer(checkpoint_every)
  n_workers <- min(as.integer(n_cores), length(indices), batch_size)
  if (verbose) message("[PWAS_Time] Fitting ", length(indices), " analytes with ", n_workers, " worker(s).")
  worker <- function(i) .pwas_analyte(i, validated, spec)
  result$metadata <- list(pipeline_version = "0.3.0", specification = spec,
    specification_md5 = .pwas_object_md5(spec), preprocessing = preprocessing,
    input_qc = validated$qc, code = code_metadata,
    fixed_formula = paste(deparse(.pwas_formula(spec, random = FALSE)), collapse = " "),
    full_formula = paste(deparse(.pwas_formula(spec)), collapse = " "),
    confidence_level = spec$confidence_level,
    estimation = "ML; one model per eligible analyte",
    inference = list(coefficients = "Satterthwaite t"),
    multiplicity = "BH after all analytes finish, across valid proteins separately within each coefficient term; includes flagged singular fits; withheld p-values excluded",
    started_utc = started_utc,
    requested_cores = as.integer(n_cores), dispatched_workers = n_workers,
    package_versions = package_versions,
    session_info = capture.output(utils::sessionInfo()))
  class(result) <- c("pwas_time_result", "list")
  batches <- split(indices, ceiling(seq_along(indices) / batch_size))
  for (batch in batches) {
    results <- if (n_workers == 1L) lapply(batch, worker) else
      parallel::mclapply(batch, worker, mc.cores = min(n_workers, length(batch)),
                        mc.preschedule = TRUE, mc.set.seed = FALSE)
    if (any(vapply(results, inherits, logical(1), what = "try-error")))
      stop("A worker failed outside a model fit; any earlier checkpoint is preserved.", call. = FALSE)
    result$coefficients$BH_P_VALUE <- NULL
    for (name in components) {
      table <- rbind(result[[name]], do.call(rbind, lapply(results, `[[`, name)))
      rownames(table) <- NULL
      result[[name]] <- table
    }
    result$covariance <- c(result$covariance,
      stats::setNames(lapply(results, `[[`, "covariance"), rownames(validated$omics)[batch]))
    complete <- nrow(result$model_qc) == n_analytes
    adjusted <- .pwas_adjust(result$coefficients, "TERM", "coefficients")
    result$coefficients <- adjusted$table
    if (!complete) result$coefficients$BH_P_VALUE <- NA_real_
    result$multiplicity <- adjusted$families
    finished <- Sys.time()
    result$metadata$finished_utc <- format(finished, tz = "UTC", usetz = TRUE)
    result$metadata$elapsed_seconds <- previous_elapsed + as.numeric(difftime(finished, started, units = "secs"))
    result$metadata$observed_worker_pids <- sort(unique(result$model_qc$WORKER_PID))
    if (!is.null(output_dir)) {
      result$metadata$checkpoint <- list(signature = signature, n_completed = nrow(result$model_qc),
        n_total = n_analytes, complete = complete)
      .pwas_write_files(result, output_dir)
      if (verbose) message("[PWAS_Time] Saved ", nrow(result$model_qc), "/", n_analytes, " analytes.")
    }
  }
  if (verbose) message("[PWAS_Time] Finished in ", round(result$metadata$elapsed_seconds, 1), " seconds; ",
    sum(tapply(result$coefficients$INFERENCE_OK, result$coefficients$ANALYTE_NAME, all)),
    "/", n_analytes, " analytes have complete inference.")
  result
}

summarize_pwas_time <- function(result) {
  if (!inherits(result, "pwas_time_result")) stop("Expected a pwas_time_result.", call. = FALSE)
  status <- as.data.frame(table(result$model_qc$STATUS), stringsAsFactors = FALSE)
  names(status) <- c("STATUS", "N_ANALYTES")
  list(input = result$metadata$input_qc, progress = result$metadata$checkpoint[c("n_completed", "n_total", "complete")],
       model_status = status,
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
  .pwas_write_files(result, output_dir)
}

# Only the checkpoint runner may replace its existing, matching outputs.
.pwas_write_files <- function(result, output_dir) {
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) stop("Cannot create output_dir.", call. = FALSE)
  paths <- file.path(output_dir, c(result = "result.rds", results = "results.csv"))
  names(paths) <- c("result", "results")
  temp_rds <- tempfile(pattern = ".result-", tmpdir = output_dir)
  temp_csv <- tempfile(pattern = ".results-", tmpdir = output_dir)
  on.exit(unlink(c(temp_rds, temp_csv)), add = TRUE)
  saveRDS(result, temp_rds)
  utils::write.csv(.pwas_results_table(result), temp_csv, row.names = FALSE, na = "")
  if (!file.rename(temp_csv, paths["results"])) stop("Could not finalize results.csv.", call. = FALSE)
  if (!file.rename(temp_rds, paths["result"])) stop("Could not finalize result.rds.", call. = FALSE)
  invisible(paths)
}
