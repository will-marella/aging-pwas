.pwas_formula <- function(spec, random = TRUE) {
  terms <- c("TIME_YEARS", "I(TIME_YEARS^2)", "AGE_C", "I(AGE_C^2)",
             "TIME_YEARS:AGE_C", "FEMALE", names(spec$covariates))
  if (random) terms <- c(terms, if (spec$random_effects == "intercept_slope")
    "(1 + TIME_YEARS | SUBJECT_ID)" else "(1 | SUBJECT_ID)")
  stats::as.formula(paste("RESPONSE ~", paste(terms, collapse = " + ")))
}

.pwas_capture <- function(expr) {
  warnings <- messages <- character()
  error <- NULL
  value <- tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }, message = function(m) {
      messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage")
    }), error = function(e) { error <<- conditionMessage(e); NULL })
  list(value = value, error = error, warnings = unique(warnings), messages = unique(messages))
}

.pwas_fit <- function(formula, data) {
  .pwas_capture(lmerTest::lmer(formula, data = data, REML = FALSE,
    na.action = stats::na.fail,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000L),
                              check.rankX = "stop.deficient")))
}

.pwas_diagnostics <- function(captured) {
  if (is.null(captured$value)) return(list(converged = FALSE, singular = NA,
    warnings = paste(captured$warnings, collapse = " | "), error = captured$error))
  fit <- captured$value
  opt_messages <- unlist(fit@optinfo$conv$lme4$messages, use.names = FALSE)
  # Boundary variance estimates are recorded separately from numerical convergence.
  numerical <- opt_messages[!grepl("singular", opt_messages, ignore.case = TRUE)]
  optimizer_code <- unlist(fit@optinfo$conv$opt, use.names = FALSE)
  warning_failure <- any(grepl("failed to converge|Hessian|negative eigenvalue|degenerate",
                              captured$warnings, ignore.case = TRUE))
  list(converged = all(optimizer_code == 0) && !length(numerical) && !warning_failure,
       singular = lme4::isSingular(fit, tol = 1e-4),
       warnings = paste(unique(c(captured$warnings, captured$messages, opt_messages)), collapse = " | "),
       error = if (is.null(captured$error)) "" else captured$error)
}

.pwas_coef_rows <- function(analyte, terms, status) {
  data.frame(ANALYTE_NAME = analyte, TERM = terms, ESTIMATE = NA_real_, SE = NA_real_,
    DF = NA_real_, CI_LOW = NA_real_, CI_HIGH = NA_real_, P_VALUE = NA_real_,
    INFERENCE_OK = FALSE, STATUS = status, stringsAsFactors = FALSE)
}

.pwas_analyte <- function(index, validated, spec) {
  pheno <- validated$pheno
  analyte <- rownames(validated$omics)[index]
  required <- c("TIME_YEARS", "BASELINE_AGE", "FEMALE", names(spec$covariates))
  pheno_ok <- stats::complete.cases(pheno[, required, drop = FALSE])
  abundance <- validated$omics[index, ]
  abundance_ok <- !is.na(abundance)
  d <- pheno[pheno_ok & abundance_ok, , drop = FALSE]
  d$RESPONSE <- abundance[pheno_ok & abundance_ok]
  counts <- table(d$SUBJECT_ID)
  repeated <- names(counts)[counts >= spec$min_visits_per_subject]
  singletons <- !(d$SUBJECT_ID %in% repeated)
  n_singletons <- sum(singletons)
  d <- d[!singletons, , drop = FALSE]
  d <- d[, c("SUBJECT_ID", "TIME_YEARS", "BASELINE_AGE", "AGE_C", "FEMALE",
              names(spec$covariates), "RESPONSE"), drop = FALSE]
  d$SUBJECT_ID <- factor(d$SUBJECT_ID)
  n_subjects <- nlevels(d$SUBJECT_ID)
  visits <- table(d$SUBJECT_ID)
  coverage <- data.frame(N_VISITS = integer(), N_SUBJECTS = integer())
  if (length(visits)) {
    coverage <- as.data.frame(table(as.integer(visits)), stringsAsFactors = FALSE)
    names(coverage) <- c("N_VISITS", "N_SUBJECTS")
    coverage$N_VISITS <- as.integer(as.character(coverage$N_VISITS))
  }
  coverage <- cbind(ANALYTE_NAME = rep(analyte, nrow(coverage)), coverage)
  exclusions <- data.frame(ANALYTE_NAME = analyte,
    REASON = c("missing_phenotype_covariate", "missing_abundance_after_phenotype_filter",
               "fewer_than_two_retained_visits"),
    N_SAMPLES = c(sum(!pheno_ok), sum(pheno_ok & !abundance_ok), n_singletons),
    stringsAsFactors = FALSE)
  fixed_formula <- .pwas_formula(spec, random = FALSE)
  # Retain declared factor levels so an absent level is a visible rank deficiency.
  design <- .pwas_capture(stats::model.matrix(fixed_formula, data = d))
  terms <- if (!is.null(design$value)) colnames(design$value) else character()
  if (!length(terms)) {
    template <- pheno[, c("TIME_YEARS", "AGE_C", "FEMALE", names(spec$covariates)), drop = FALSE]
    template$RESPONSE <- 0
    terms <- colnames(stats::model.matrix(fixed_formula, data = template))
  }
  status <- "ok"
  detail <- ""
  rank <- NA_integer_
  if (n_subjects < spec$min_subjects) {
    status <- "insufficient_subjects"
  } else if (is.null(design$value)) {
    status <- "design_error"; detail <- design$error
  } else {
    X <- design$value
    rank <- qr(X)$rank
    if (rank < ncol(X)) status <- "rank_deficient"
    else if (nrow(X) <= ncol(X)) status <- "insufficient_observations"
    else if (spec$random_effects == "intercept_slope" && nrow(d) <= 2L * n_subjects)
      status <- "insufficient_random_effect_observations"
    else if (length(unique(d$RESPONSE)) < 2L) status <- "constant_abundance"
  }
  qc <- data.frame(ANALYTE_NAME = analyte, STATUS = status, N_INPUT_SAMPLES = nrow(pheno),
    N_OBS = nrow(d), N_SUBJECTS = n_subjects,
    N_SUBJECTS_WITH_BASELINE = length(unique(d$SUBJECT_ID[d$TIME_YEARS == 0])),
    N_SUBJECTS_WITH_3PLUS_VISITS = sum(visits >= 3L),
    N_DISTINCT_TIMES = length(unique(d$TIME_YEARS)),
    TIME_MIN = if (nrow(d)) min(d$TIME_YEARS) else NA_real_,
    TIME_MAX = if (nrow(d)) max(d$TIME_YEARS) else NA_real_,
    BASELINE_AGE_MIN = if (nrow(d)) min(d$BASELINE_AGE) else NA_real_,
    BASELINE_AGE_MAX = if (nrow(d)) max(d$BASELINE_AGE) else NA_real_,
    FIXED_RANK = rank, N_FIXED = length(terms), FULL_CONVERGED = NA,
    FULL_SINGULAR = NA, RESIDUAL_SD = NA_real_,
    WARNINGS = "", ERROR = detail, WORKER_PID = Sys.getpid(), stringsAsFactors = FALSE)
  coefficients <- .pwas_coef_rows(analyte, terms, status)
  covariance <- NULL
  fit <- NULL
  if (status == "ok") {
    captured <- .pwas_fit(.pwas_formula(spec), d)
    fit_diag <- .pwas_diagnostics(captured)
    fit <- captured$value
    status <- if (is.null(fit)) "fit_error" else if (!fit_diag$converged) "nonconverged" else
      if (fit_diag$singular) "singular" else "ok"
    qc$FULL_CONVERGED <- fit_diag$converged; qc$FULL_SINGULAR <- fit_diag$singular
    qc$WARNINGS <- fit_diag$warnings; qc$ERROR <- if (is.null(fit_diag$error)) "" else fit_diag$error
    extracted_fit <- if (is.null(fit)) NULL else .pwas_capture(list(
      covariance = as.matrix(stats::vcov(fit)), beta = lme4::fixef(fit),
      residual_sd = stats::sigma(fit)))
    if (!is.null(extracted_fit) && is.null(extracted_fit$value)) {
      status <- "fit_summary_failed"
      qc$ERROR <- extracted_fit$error
      fit <- NULL
    }
    if (!is.null(fit)) {
      covariance <- extracted_fit$value$covariance
      beta <- extracted_fit$value$beta
      qc$RESIDUAL_SD <- extracted_fit$value$residual_sd
      if (any(!is.finite(covariance)) || any(diag(covariance) <= 0) || any(!is.finite(beta))) {
        status <- "invalid_covariance"
      }
      coefficients$ESTIMATE <- unname(beta[coefficients$TERM])
      coefficients$SE <- sqrt(diag(covariance))[coefficients$TERM]
      if (status == "ok") {
        extracted <- .pwas_capture(summary(fit, ddf = "Satterthwaite")$coefficients)
        if (is.null(extracted$value) || length(extracted$warnings)) {
          status <- "coefficient_inference_failed"
        } else {
          ct <- extracted$value[coefficients$TERM, , drop = FALSE]
          coefficients$DF <- ct[, "df"]
          coefficients$P_VALUE <- ct[, "Pr(>|t|)"]
          critical <- stats::qt((1 + spec$confidence_level) / 2, coefficients$DF)
          coefficients$CI_LOW <- coefficients$ESTIMATE - critical * coefficients$SE
          coefficients$CI_HIGH <- coefficients$ESTIMATE + critical * coefficients$SE
          coefficients$INFERENCE_OK <- is.finite(coefficients$P_VALUE) &
            coefficients$DF > 0 & is.finite(coefficients$CI_LOW)
          if (!all(coefficients$INFERENCE_OK)) status <- "coefficient_inference_failed"
        }
      }
    }
  }
  if (status != "ok") {
    coefficients$DF <- coefficients$CI_LOW <- coefficients$CI_HIGH <- coefficients$P_VALUE <- NA_real_
    coefficients$INFERENCE_OK <- FALSE
  }
  coefficients$STATUS <- status
  qc$STATUS <- status
  list(coefficients = coefficients, model_qc = qc,
       visit_coverage = coverage, exclusions = exclusions,
       covariance = covariance)
}
