# Shared, versioned model specification. No cohort data preparation occurs here.

.pwas_is_numeric <- function(x) {
  is.numeric(x) && !is.complex(x) && is.null(dim(x))
}

.pwas_nonempty_strings <- function(x) {
  is.character(x) && is.null(dim(x)) && !anyNA(x) && all(nzchar(trimws(x)))
}

.pwas_named_list <- function(x) {
  is.list(x) && !is.data.frame(x) &&
    (length(x) == 0L ||
       (!is.null(names(x)) && .pwas_nonempty_strings(names(x)) &&
          !anyDuplicated(names(x))))
}

#' Construct the longitudinal PWAS model specification.
#'
#' Baseline age is centered at a required, shared constant. Additional covariates
#' have explicit numeric, logical, or factor types and baseline/time-varying
#' timing. Factor declarations require levels and a reference level.
pwas_time_spec <- function(age_center, covariates = list(),
                           random_effects = "intercept_slope",
                           min_subjects = 20L,
                           confidence_level = 0.95) {
  if (missing(age_center) || !.pwas_is_numeric(age_center) ||
      length(age_center) != 1L || !is.finite(age_center)) {
    stop("age_center must be one finite numeric value.", call. = FALSE)
  }
  if (!is.character(random_effects) || length(random_effects) != 1L ||
      is.na(random_effects) ||
      !random_effects %in% c("intercept_slope", "intercept")) {
    stop("random_effects must be 'intercept_slope' or 'intercept'.",
         call. = FALSE)
  }
  if (!.pwas_is_numeric(min_subjects) || length(min_subjects) != 1L ||
      !is.finite(min_subjects) || min_subjects < 2 ||
      min_subjects > .Machine$integer.max ||
      min_subjects != floor(min_subjects)) {
    stop("min_subjects must be an integer of at least two.", call. = FALSE)
  }
  if (!.pwas_is_numeric(confidence_level) || length(confidence_level) != 1L ||
      !is.finite(confidence_level) || confidence_level <= 0 ||
      confidence_level >= 1) {
    stop("confidence_level must lie strictly between zero and one.",
         call. = FALSE)
  }
  if (!.pwas_named_list(covariates)) {
    stop("covariates must be a uniquely named list of declarations.",
         call. = FALSE)
  }
  reserved <- c("SAMPLE_ID", "SUBJECT_ID", "TIME_YEARS", "BASELINE_AGE",
                "AGE_C", "FEMALE", "ANALYTE_NAME", "RESPONSE", "PROTEIN",
                "protein", ".protein", ".response", ".", "...")
  covariate_names <- names(covariates)
  if (length(covariates) &&
      (any(make.names(covariate_names) != covariate_names) ||
       any(covariate_names %in% reserved) ||
       any(grepl("^\\.\\.[0-9]+$", covariate_names)))) {
    stop("Covariate names must be syntactic R names and cannot use reserved fields.",
         call. = FALSE)
  }
  for (j in seq_along(covariates)) {
    declaration <- covariates[[j]]
    if (!.pwas_named_list(declaration) ||
        !all(c("type", "timing") %in% names(declaration))) {
      stop("Every covariate requires a named declaration with type and timing.",
           call. = FALSE)
    }
    if (!is.character(declaration$type) || length(declaration$type) != 1L ||
        is.na(declaration$type) ||
        !declaration$type %in% c("numeric", "factor", "logical")) {
      stop("Covariate type must be 'numeric', 'factor', or 'logical'.",
           call. = FALSE)
    }
    if (!is.character(declaration$timing) || length(declaration$timing) != 1L ||
        is.na(declaration$timing) ||
        !declaration$timing %in% c("baseline", "time_varying")) {
      stop("Covariate timing must be 'baseline' or 'time_varying'.",
           call. = FALSE)
    }
    allowed <- if (declaration$type == "factor") {
      c("type", "timing", "levels", "reference")
    } else {
      c("type", "timing")
    }
    if (!setequal(names(declaration), allowed)) {
      stop("Covariate declaration has missing or unknown fields.",
           call. = FALSE)
    }
    if (declaration$type == "factor") {
      if (!.pwas_nonempty_strings(declaration$levels) ||
          length(declaration$levels) < 2L || anyDuplicated(declaration$levels)) {
        stop("Factor levels must be at least two unique nonempty strings.",
             call. = FALSE)
      }
      if (!.pwas_nonempty_strings(declaration$reference) ||
          length(declaration$reference) != 1L ||
          !declaration$reference %in% declaration$levels) {
        stop("Each factor reference must be one declared level.",
             call. = FALSE)
      }
      declaration$levels <- c(declaration$reference,
                              setdiff(declaration$levels, declaration$reference))
    }
    covariates[[j]] <- declaration[allowed]
  }
  structure(list(
    version = "0.2.1",
    age_center = as.numeric(age_center),
    covariates = covariates,
    random_effects = random_effects,
    min_subjects = as.integer(min_subjects),
    confidence_level = as.numeric(confidence_level),
    time_unit = "years",
    min_visits_per_subject = 2L,
    fit_method = "ML",
    coefficient_inference = "Satterthwaite",
    singular_policy = "flag_only",
    convergence_policy = "withhold_inference",
    fallback = "none"
  ), class = "pwas_time_spec")
}

# Reconstruct caller-editable fields and reject mutations of fixed policies.
.validate_pwas_spec <- function(spec) {
  if (!inherits(spec, "pwas_time_spec") || !.pwas_named_list(spec)) {
    stop("spec must be constructed with pwas_time_spec().", call. = FALSE)
  }
  editable <- c("age_center", "covariates", "random_effects", "min_subjects",
                "confidence_level")
  fixed <- c("version", "time_unit", "min_visits_per_subject", "fit_method",
             "coefficient_inference", "singular_policy", "convergence_policy",
             "fallback")
  if (!setequal(names(spec), c(editable, fixed))) {
    stop("spec has missing or unknown fields; reconstruct it with pwas_time_spec().",
         call. = FALSE)
  }
  canonical <- do.call(pwas_time_spec, spec[editable])
  if (!all(vapply(fixed, function(field) {
    identical(spec[[field]], canonical[[field]])
  }, logical(1L)))) {
    stop("Fixed specification policies cannot be changed; use the versioned model specification.",
         call. = FALSE)
  }
  canonical
}
