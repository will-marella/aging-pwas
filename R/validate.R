# Input validation is deliberately separate from analyte-specific exclusions.

.pwas_check_column_names <- function(data, label) {
  if (!.pwas_nonempty_strings(names(data)) || anyDuplicated(names(data))) {
    stop(paste0(label, " must have unique nonempty column names."), call. = FALSE)
  }
}

.pwas_constant_within_subject <- function(values, subject) {
  all(vapply(split(values, subject), function(x) {
    length(unique(x[!is.na(x)])) <= 1L
  }, logical(1L)))
}

.pwas_factor <- function(values, levels) {
  value <- factor(as.character(values), levels = levels)
  contrasts(value) <- stats::contr.treatment(levels, base = 1L)
  value
}

#' Validate and align analysis-ready longitudinal PWAS inputs.
#'
#' The original subject baseline is checked before excluding any observations.
#' No abundance transformation, imputation, batch correction, or row exclusion is
#' performed here. Validation errors never include sample or subject identifiers.
validate_pwas_inputs <- function(pheno, omics, spec, preprocessing) {
  spec <- .validate_pwas_spec(spec)
  preprocessing_fields <- c("abundance_scale", "normalization",
                            "batch_handling", "missing_values")
  if (missing(preprocessing) || !.pwas_named_list(preprocessing) ||
      !all(preprocessing_fields %in% names(preprocessing)) ||
      !all(vapply(preprocessing[preprocessing_fields], function(x) {
        .pwas_nonempty_strings(x) && length(x) == 1L
      }, logical(1L)))) {
    stop(paste0("preprocessing requires nonempty string declarations for ",
                "abundance_scale, normalization, batch_handling, and missing_values."),
         call. = FALSE)
  }
  if (!is.data.frame(pheno) || nrow(pheno) == 0L) {
    stop("pheno must be a nonempty data.frame.", call. = FALSE)
  }
  pheno <- as.data.frame(pheno, stringsAsFactors = FALSE)
  .pwas_check_column_names(pheno, "pheno")
  required <- c("SAMPLE_ID", "SUBJECT_ID", "TIME_YEARS", "BASELINE_AGE", "FEMALE")
  if (!all(c(required, names(spec$covariates)) %in% names(pheno))) {
    stop("pheno is missing required or declared covariate columns.",
         call. = FALSE)
  }
  for (field in c("SAMPLE_ID", "SUBJECT_ID")) {
    if (!.pwas_nonempty_strings(pheno[[field]])) {
      stop("Sample and subject identifiers must be nonmissing, nonempty character vectors.",
           call. = FALSE)
    }
  }
  if (anyDuplicated(pheno$SAMPLE_ID)) {
    stop("pheno contains duplicate samples.", call. = FALSE)
  }
  if (!.pwas_is_numeric(pheno$TIME_YEARS) ||
      any(!is.finite(pheno$TIME_YEARS)) || any(pheno$TIME_YEARS < 0)) {
    stop("TIME_YEARS must contain finite, nonnegative numeric elapsed years.",
         call. = FALSE)
  }
  if (!.pwas_is_numeric(pheno$BASELINE_AGE) ||
      any(!is.finite(pheno$BASELINE_AGE)) ||
      any(pheno$BASELINE_AGE < 0 | pheno$BASELINE_AGE > 120)) {
    stop("BASELINE_AGE must contain finite numeric ages between 0 and 120 years.",
         call. = FALSE)
  }
  if (anyDuplicated(pheno[c("SUBJECT_ID", "TIME_YEARS")])) {
    stop("pheno contains duplicate participant-time pairs.", call. = FALSE)
  }
  if (!.pwas_constant_within_subject(pheno$BASELINE_AGE, pheno$SUBJECT_ID)) {
    stop("BASELINE_AGE must remain constant within every participant.",
         call. = FALSE)
  }
  if (!all(vapply(split(pheno$TIME_YEARS, pheno$SUBJECT_ID), function(x) {
    any(x == 0)
  }, logical(1L)))) {
    stop("Every participant must have a defined TIME_YEARS = 0 sample in the input.",
         call. = FALSE)
  }
  if (.pwas_is_numeric(pheno$FEMALE)) {
    if (any(!is.na(pheno$FEMALE) & !pheno$FEMALE %in% c(0, 1))) {
      stop("FEMALE must be coded 0/1, with missing values allowed.", call. = FALSE)
    }
  } else if (is.factor(pheno$FEMALE)) {
    if (any(!levels(pheno$FEMALE) %in% c("0", "1")) ||
        any(!is.na(pheno$FEMALE) & !as.character(pheno$FEMALE) %in% c("0", "1"))) {
      stop("FEMALE factor labels must be 0/1, with missing values allowed.",
           call. = FALSE)
    }
  } else {
    stop("FEMALE must be numeric 0/1 or a factor with labels 0/1.", call. = FALSE)
  }
  if (!.pwas_constant_within_subject(pheno$FEMALE, pheno$SUBJECT_ID)) {
    stop("Observed FEMALE values must remain constant within every participant.",
         call. = FALSE)
  }
  pheno$FEMALE <- .pwas_factor(pheno$FEMALE, c("0", "1"))
  for (field in names(spec$covariates)) {
    declaration <- spec$covariates[[field]]
    values <- pheno[[field]]
    if (declaration$type == "numeric") {
      if (!.pwas_is_numeric(values) ||
          any(!is.na(values) & !is.finite(values))) {
        stop("A declared numeric covariate has a nonnumeric or infinite value.",
             call. = FALSE)
      }
    } else if (declaration$type == "logical") {
      if (!is.logical(values) || !is.null(dim(values))) {
        stop("A declared logical covariate is not a logical vector.", call. = FALSE)
      }
      # Explicit FALSE = 0 / TRUE = 1 coding avoids session-dependent factor
      # contrasts when model.matrix processes logical predictors.
      pheno[[field]] <- as.integer(values)
    } else {
      if ((!is.factor(values) && !is.character(values)) ||
          !is.null(dim(values)) ||
          (is.factor(values) && any(!levels(values) %in% declaration$levels)) ||
          any(!is.na(values) & !as.character(values) %in% declaration$levels)) {
        stop("A declared factor covariate has an invalid type or an undeclared label.",
             call. = FALSE)
      }
      pheno[[field]] <- .pwas_factor(values, declaration$levels)
    }
    if (declaration$timing == "baseline" &&
        !.pwas_constant_within_subject(values, pheno$SUBJECT_ID)) {
      stop("Observed baseline covariates must remain constant within every participant.",
           call. = FALSE)
    }
  }
  if (!is.data.frame(omics) || nrow(omics) == 0L || ncol(omics) < 2L) {
    stop("omics must be a nonempty data.frame with analyte and sample columns.",
         call. = FALSE)
  }
  omics <- as.data.frame(omics, stringsAsFactors = FALSE)
  .pwas_check_column_names(omics, "omics")
  if (names(omics)[1L] != "ANALYTE_NAME" ||
      !.pwas_nonempty_strings(omics[[1L]]) || anyDuplicated(omics[[1L]])) {
    stop("The first omics column must be ANALYTE_NAME with unique nonempty character values.",
         call. = FALSE)
  }
  if (!setequal(names(omics)[-1L], pheno$SAMPLE_ID)) {
    stop("pheno and omics sample sets must match exactly.", call. = FALSE)
  }
  if (!all(vapply(omics[-1L], .pwas_is_numeric, logical(1L)))) {
    stop("Every omics sample column must be numeric.", call. = FALSE)
  }
  if (any(vapply(omics[-1L], function(x) any(is.infinite(x)), logical(1L)))) {
    stop("omics contains infinite values; only finite values and missing values are allowed.",
         call. = FALSE)
  }
  pheno <- pheno[order(pheno$SUBJECT_ID, pheno$TIME_YEARS,
                       pheno$SAMPLE_ID, method = "radix"), , drop = FALSE]
  rownames(pheno) <- NULL
  pheno$AGE_C <- pheno$BASELINE_AGE - spec$age_center
  abundance <- as.matrix(omics[pheno$SAMPLE_ID])
  storage.mode(abundance) <- "double"
  rownames(abundance) <- omics$ANALYTE_NAME
  colnames(abundance) <- pheno$SAMPLE_ID
  covariate_fields <- c("FEMALE", names(spec$covariates))
  missing_covariate_rows <- !stats::complete.cases(pheno[covariate_fields])
  visit_counts <- as.integer(table(pheno$SUBJECT_ID))
  visit_distribution <- as.data.frame(table(visit_counts), stringsAsFactors = FALSE)
  names(visit_distribution) <- c("n_visits", "n_subjects")
  visit_distribution$n_visits <- as.integer(as.character(visit_distribution$n_visits))
  qc <- list(
    n_input_samples = nrow(pheno),
    n_input_subjects = length(unique(pheno$SUBJECT_ID)),
    n_analytes = nrow(abundance),
    n_missing_abundances = sum(is.na(abundance)),
    n_rows_missing_covariates = sum(missing_covariate_rows),
    n_subjects_with_missing_covariates = length(unique(pheno$SUBJECT_ID[missing_covariate_rows])),
    n_missing_by_covariate = vapply(pheno[covariate_fields], function(x) sum(is.na(x)), integer(1L)),
    n_subjects_one_visit = sum(visit_counts == 1L),
    visit_count_distribution = visit_distribution,
    time_range_years = range(pheno$TIME_YEARS),
    baseline_age_range_years = range(pheno$BASELINE_AGE)
  )
  list(pheno = pheno, omics = abundance, qc = qc)
}
