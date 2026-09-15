fixture <- make_synthetic_pwas()
spec_fixture <- spec
fit_fixture <- run_pwas_time(fixture$pheno, fixture$omics, spec_fixture,
                             fixture$preprocessing, verbose = FALSE)

run_fixture <- function(pheno = fixture$pheno, omics = fixture$omics,
                        specification = spec_fixture, cores = 1L) {
  run_pwas_time(pheno, omics, specification, fixture$preprocessing,
                n_cores = cores, verbose = FALSE)
}

testthat::test_that("synthetic trajectories recover known temporal effects", {
  testthat::expect_true(all(fit_fixture$model_qc$STATUS == "ok"))
  estimate <- merge(fit_fixture$coefficients, fixture$truth,
    by = c("ANALYTE_NAME", "TERM"), suffixes = c("", "_TRUE"))
  for (term in c("TIME_YEARS", "I(TIME_YEARS^2)", "TIME_YEARS:AGE_C")) {
    x <- estimate[estimate$TERM == term, ]
    tolerance <- switch(term, TIME_YEARS = 0.03, "I(TIME_YEARS^2)" = 0.001,
                        "TIME_YEARS:AGE_C" = 0.0015)
    testthat::expect_lt(max(abs(x$ESTIMATE - x$ESTIMATE_TRUE)), tolerance)
  }
})

testthat::test_that("estimates, covariance, and coefficient inference agree with an independent direct fit", {
  d <- fixture$pheno
  d$AGE_C <- d$BASELINE_AGE - 50
  d$FEMALE <- factor(d$FEMALE, levels = c(0, 1))
  d$RESPONSE <- as.numeric(fixture$omics[6L, -1L])
  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000L))
  direct <- lmerTest::lmer(RESPONSE ~ TIME_YEARS + I(TIME_YEARS^2) + AGE_C +
    I(AGE_C^2) + FEMALE + SITE + TIME_YEARS:AGE_C + (1 + TIME_YEARS | SUBJECT_ID),
    data = d, REML = FALSE, control = ctrl)
  reported <- subset(fit_fixture$coefficients, ANALYTE_NAME == "synthetic_combined")
  testthat::expect_equal(reported$ESTIMATE, unname(lme4::fixef(direct)[reported$TERM]), tolerance = 1e-7)
  V <- fit_fixture$covariance$synthetic_combined
  testthat::expect_equal(V, as.matrix(stats::vcov(direct))[rownames(V), colnames(V)], tolerance = 1e-7)
  ct <- summary(direct)$coefficients[reported$TERM, ]
  testthat::expect_equal(reported$SE, unname(ct[, "Std. Error"]), tolerance = 1e-7)
  testthat::expect_equal(reported$DF, unname(ct[, "df"]), tolerance = 1e-7)
  testthat::expect_equal(reported$P_VALUE, unname(ct[, "Pr(>|t|)"]), tolerance = 1e-7)
  testthat::expect_equal(reported$CI_LOW,
    reported$ESTIMATE - qt(.975, reported$DF) * reported$SE)
  testthat::expect_equal(reported$CI_HIGH,
    reported$ESTIMATE + qt(.975, reported$DF) * reported$SE)
})

testthat::test_that("sample order and session contrast options cannot change the result", {
  set.seed(731)
  p <- fixture$pheno[sample(nrow(fixture$pheno)), ]
  o <- fixture$omics[, c(1L, sample(seq.int(2L, ncol(fixture$omics))))]
  old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old))
  reordered <- run_fixture(p, o)
  for (component in c("coefficients", "exclusions", "covariance"))
    testthat::expect_equal(reordered[[component]], fit_fixture[[component]], tolerance = 1e-10)
})

testthat::test_that("missingness is sequential, reconciled, and never resets baseline", {
  p <- fixture$pheno
  o <- fixture$omics[1L, ]
  o[1L, p$SAMPLE_ID[1L]] <- NA_real_ # First participant retains t=5,10,15.
  o[1L, p$SAMPLE_ID[5:7]] <- NA_real_ # Second participant's remaining singleton is excluded.
  p$FEMALE[9L] <- NA_integer_ # Third participant retains three visits.
  result <- run_fixture(p, o)
  testthat::expect_equal(result$model_qc$N_OBS, 634L)
  testthat::expect_equal(result$model_qc$N_SUBJECTS, 159L)
  testthat::expect_equal(result$model_qc$N_SUBJECTS_WITH_BASELINE, 157L)
  testthat::expect_equal(result$exclusions$N_SAMPLES, c(1L, 4L, 1L))
  testthat::expect_equal(sum(result$exclusions$N_SAMPLES) + result$model_qc$N_OBS, nrow(p))
  valid <- validate_pwas_inputs(p, o, spec_fixture, fixture$preprocessing)
  testthat::expect_equal(valid$pheno$AGE_C, valid$pheno$BASELINE_AGE - 50)
  # Direct fit on the explicitly retained rows verifies unchanged t=5/10/15.
  d <- valid$pheno[-c(1L, 5:9), ]
  d$RESPONSE <- as.numeric(o[1L, d$SAMPLE_ID])
  direct <- .pwas_fit(.pwas_formula(spec_fixture), d)$value
  testthat::expect_equal(result$coefficients$ESTIMATE,
    unname(lme4::fixef(direct)[result$coefficients$TERM]), tolerance = 1e-7)
})

testthat::test_that("all-missing and constant analytes remain in every fit ledger", {
  o <- fixture$omics[1:2, ]
  bad <- o
  bad$ANALYTE_NAME <- c("all_missing", "constant")
  bad[1L, -1L] <- NA_real_; bad[2L, -1L] <- 2
  result <- run_fixture(omics = rbind(o, bad))
  testthat::expect_equal(nrow(result$model_qc), 4L)
  testthat::expect_equal(tail(result$model_qc$STATUS, 2), c("insufficient_subjects", "constant_abundance"))
  testthat::expect_equal(length(result$covariance), 4L)
  testthat::expect_null(result$covariance$all_missing)
  testthat::expect_true(all(is.na(subset(result$coefficients, ANALYTE_NAME %in% c("all_missing", "constant"))$P_VALUE)))
  testthat::expect_equal(result$model_qc$N_OBS[3L], 0L)
  m <- result$multiplicity
  testthat::expect_equal(m$N_REQUESTED, rep(4L, nrow(m)))
  testthat::expect_equal(m$N_WITHHELD, rep(2L, nrow(m)))
  testthat::expect_equal(sum(subset(result$exclusions, ANALYTE_NAME == "all_missing")$N_SAMPLES), 640)
  out <- tempfile("pwas-failed-export-")
  on.exit(unlink(out, recursive = TRUE))
  write_pwas_time(result, out)
  exported <- read.csv(file.path(out, "results.csv"), check.names = FALSE)
  testthat::expect_equal(exported$ANALYTE_NAME, result$model_qc$ANALYTE_NAME)
  testthat::expect_equal(tail(exported$STATUS, 2), c("insufficient_subjects", "constant_abundance"))
  testthat::expect_true(all(is.na(tail(exported[["TIME_YEARS__P_VALUE"]], 2))))
})

testthat::test_that("two common visits cannot estimate the quadratic model", {
  p <- subset(fixture$pheno, TIME_YEARS %in% c(0, 10))
  o <- fixture$omics[1L, c("ANALYTE_NAME", p$SAMPLE_ID)]
  result <- run_fixture(p, o)
  testthat::expect_equal(result$model_qc$STATUS, "rank_deficient")
  testthat::expect_true(all(is.na(result$coefficients$P_VALUE)))
})

testthat::test_that("singular random slopes retain direct-fit inference and BH adjustment with flags", {
  p <- fixture$pheno
  t <- p$TIME_YEARS
  i <- rep(seq_len(160), each = 4)
  y <- sin(i) + .04 * t + rep(c(.2, -.2, -.2, .2), 160) * (-1)^i
  o <- fixture$omics[1:2, ]; o[1L, -1L] <- y
  result <- testthat::expect_message(
    run_pwas_time(p, o, spec_fixture, fixture$preprocessing),
    "2/2 analytes have complete inference")
  testthat::expect_equal(result$model_qc$FULL_SINGULAR, c(TRUE, FALSE))
  testthat::expect_equal(result$model_qc$STATUS, c("singular", "ok"))
  reported <- subset(result$coefficients, ANALYTE_NAME == o$ANALYTE_NAME[1L])
  testthat::expect_true(all(reported$STATUS == "singular"))
  testthat::expect_true(all(reported$INFERENCE_OK))
  testthat::expect_true(all(is.finite(as.matrix(reported[c("DF", "CI_LOW", "CI_HIGH", "P_VALUE")]))))

  d <- p
  d$AGE_C <- d$BASELINE_AGE - 50
  d$FEMALE <- factor(d$FEMALE, levels = c(0, 1))
  d$RESPONSE <- y
  direct <- suppressMessages(lmerTest::lmer(
    RESPONSE ~ TIME_YEARS + I(TIME_YEARS^2) + AGE_C + I(AGE_C^2) +
      TIME_YEARS:AGE_C + FEMALE + SITE + (1 + TIME_YEARS | SUBJECT_ID),
    data = d, REML = FALSE,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000L))))
  testthat::expect_true(lme4::isSingular(direct))
  ct <- summary(direct, ddf = "Satterthwaite")$coefficients[reported$TERM, ]
  testthat::expect_equal(reported$P_VALUE, unname(ct[, "Pr(>|t|)"]), tolerance = 1e-7)
  testthat::expect_equal(reported$DF, unname(ct[, "df"]), tolerance = 1e-7)
  testthat::expect_equal(reported$CI_LOW,
    unname(ct[, "Estimate"] - qt(.975, ct[, "df"]) * ct[, "Std. Error"]), tolerance = 1e-7)
  testthat::expect_equal(reported$CI_HIGH,
    unname(ct[, "Estimate"] + qt(.975, ct[, "df"]) * ct[, "Std. Error"]), tolerance = 1e-7)
  for (term in unique(result$coefficients$TERM)) {
    tab <- subset(result$coefficients, TERM == term)
    testthat::expect_equal(tab$BH_P_VALUE, p.adjust(tab$P_VALUE, "BH"))
  }
  testthat::expect_true(all(result$multiplicity$N_VALID_P == 2L))
  testthat::expect_true(all(result$multiplicity$N_WITHHELD == 0L))
  testthat::expect_equal(result$metadata$specification$random_effects, "intercept_slope")
  testthat::expect_equal(result$metadata$specification$singular_policy, "flag_only")

  out <- tempfile("pwas-singular-export-")
  on.exit(unlink(out, recursive = TRUE))
  write_pwas_time(result, out)
  exported <- read.csv(file.path(out, "results.csv"), check.names = FALSE)
  testthat::expect_equal(exported$FULL_SINGULAR, c(TRUE, FALSE))
  testthat::expect_equal(exported$STATUS, c("singular", "ok"))
  testthat::expect_true(all(is.finite(exported[["TIME_YEARS__P_VALUE"]])))
  testthat::expect_true(all(is.finite(exported[["TIME_YEARS__CI_LOW"]])))
  testthat::expect_true(all(is.finite(exported[["TIME_YEARS__CI_HIGH"]])))
})

testthat::test_that("each BH family uses exactly its valid p values", {
  for (term in unique(fit_fixture$coefficients$TERM)) {
    tab <- subset(fit_fixture$coefficients, TERM == term)
    testthat::expect_equal(tab$BH_P_VALUE, p.adjust(tab$P_VALUE, "BH"))
  }
})

with_fit_override <- function(replacement, expr) {
  env <- environment(run_pwas_time)
  original <- get(".pwas_fit", envir = env)
  assign(".pwas_fit", replacement, envir = env)
  on.exit(assign(".pwas_fit", original, envir = env))
  force(expr)
}

testthat::test_that("each eligible protein gets exactly one full-model fit", {
  original <- .pwas_fit
  calls <- 0L
  replacement <- function(formula, data) {
    calls <<- calls + 1L
    testthat::expect_equal(deparse(formula), deparse(.pwas_formula(spec_fixture)))
    original(formula, data)
  }
  result <- with_fit_override(replacement, run_fixture(omics = fixture$omics[1:2, ]))
  testthat::expect_equal(calls, 2L)
  testthat::expect_false(any(c("tests", "contrasts") %in% names(result)))
  testthat::expect_false(any(c("contrasts", "contrast_inference", "association_inference") %in%
                             names(result$metadata$specification)))
  testthat::expect_equal(result$metadata$inference, list(coefficients = "Satterthwaite t"))
})

testthat::test_that("nonconvergence remains distinct from singularity and withholds inference", {
  original <- .pwas_fit
  replacement <- function(formula, data) {
    captured <- original(formula, data)
    captured$value@optinfo$conv$lme4$messages <- "Model failed to converge (injected diagnostic)"
    captured
  }
  result <- with_fit_override(replacement, run_fixture(omics = fixture$omics[1L, ]))
  testthat::expect_equal(result$model_qc$STATUS, "nonconverged")
  testthat::expect_false(result$model_qc$FULL_CONVERGED)
  testthat::expect_false(result$model_qc$FULL_SINGULAR)
  testthat::expect_true(all(is.na(result$coefficients$P_VALUE)))
})

testthat::test_that("input errors are explicit and do not reveal identifiers", {
  check <- function(p = fixture$pheno, o = fixture$omics, pattern) {
    testthat::expect_error(run_fixture(p, o), pattern)
  }
  p <- fixture$pheno; p$SAMPLE_ID[2L] <- p$SAMPLE_ID[1L]
  check(p = p, pattern = "duplicate samples")
  check(o = fixture$omics[, -2L], pattern = "match exactly")
  p <- fixture$pheno; p$TIME_YEARS[1L] <- 1
  check(p = p, pattern = "TIME_YEARS = 0")
  p <- fixture$pheno; p$BASELINE_AGE[2L] <- p$BASELINE_AGE[2L] + 1
  check(p = p, pattern = "BASELINE_AGE must remain constant")
  p <- fixture$pheno; p$TIME_YEARS[2L] <- 0
  check(p = p, pattern = "duplicate participant-time")
  o <- fixture$omics; o[1L, 2L] <- Inf
  check(o = o, pattern = "infinite")
  p <- fixture$pheno; p$SITE <- as.character(p$SITE); p$SITE[1L] <- "UNKNOWN"
  check(p = p, pattern = "undeclared label")
  testthat::expect_error(pwas_time_spec(50, covariates = list(TIME_YEARS = list(type = "numeric", timing = "baseline"))), "reserved")
  for (name in c(".", "...", "..1")) {
    declaration <- stats::setNames(list(list(type = "numeric", timing = "time_varying")), name)
    testthat::expect_error(pwas_time_spec(50, covariates = declaration), "reserved")
  }
  s <- spec_fixture; s$fit_method <- "REML"
  testthat::expect_error(run_fixture(specification = s), "Fixed specification policies")
  testthat::expect_error(run_pwas_time(fixture$pheno, fixture$omics, spec_fixture, list()), "preprocessing requires")
  testthat::expect_error(run_fixture(cores = 0), "positive integer")
  p <- fixture$pheno; p$SITEB <- rep(seq_len(160), each = 4)
  s <- pwas_time_spec(50, covariates = c(spec_fixture$covariates,
    list(SITEB = list(type = "numeric", timing = "baseline"))))
  testthat::expect_error(run_fixture(p, specification = s), "duplicate model coefficient names")
  p <- fixture$pheno; p$SAMPLE_ID[1L] <- "secret_example_identifier"
  error <- tryCatch(run_fixture(p), error = conditionMessage)
  testthat::expect_false(grepl("secret_example_identifier", error, fixed = TRUE))
})

testthat::test_that("serialization is inspectable and cannot overwrite earlier results", {
  out <- tempfile("pwas-test-")
  on.exit(unlink(out, recursive = TRUE))
  # The CSV must align by analyte and term even if coefficient rows are reordered.
  shuffled <- fit_fixture
  shuffled$coefficients <- shuffled$coefficients[rev(seq_len(nrow(shuffled$coefficients))), ]
  paths <- write_pwas_time(shuffled, out)
  testthat::expect_true(all(file.exists(paths)))
  testthat::expect_setequal(list.files(out), c("result.rds", "results.csv"))
  restored <- readRDS(file.path(out, "result.rds"))
  testthat::expect_equal(restored, shuffled)
  exported <- read.csv(file.path(out, "results.csv"), check.names = FALSE)
  testthat::expect_equal(nrow(exported), nrow(fixture$omics))
  testthat::expect_equal(exported$ANALYTE_NAME, fit_fixture$model_qc$ANALYTE_NAME)
  testthat::expect_equal(exported$N_OBS, fit_fixture$model_qc$N_OBS)
  for (term in unique(fit_fixture$coefficients$TERM)) {
    rows <- subset(fit_fixture$coefficients, TERM == term)
    rows <- rows[match(exported$ANALYTE_NAME, rows$ANALYTE_NAME), ]
    for (field in c("ESTIMATE", "SE", "DF", "CI_LOW", "CI_HIGH", "P_VALUE", "BH_P_VALUE", "INFERENCE_OK"))
      testthat::expect_equal(exported[[paste0(term, "__", field)]], rows[[field]], tolerance = 1e-12)
  }
  testthat::expect_error(write_pwas_time(fit_fixture, out), "not empty")
  text <- paste(capture.output(dput(restored)), collapse = "\n")
  testthat::expect_false(grepl("synthetic_person_|synthetic_sample_", text))
  testthat::expect_true(all(c("lme4", "lmerTest", "Matrix") %in% names(restored$metadata$package_versions)))
  testthat::expect_equal(length(restored$metadata$code$source_md5), 5L)
})

testthat::test_that("parallel fitting uses distinct workers and reproduces serial output", {
  testthat::skip_if(.Platform$OS.type != "unix")
  result <- run_fixture(cores = 2L)
  testthat::expect_equal(length(result$metadata$observed_worker_pids), 2L)
  testthat::expect_false(Sys.getpid() %in% result$metadata$observed_worker_pids)
  for (component in c("coefficients", "covariance"))
    testthat::expect_equal(result[[component]], fit_fixture[[component]], tolerance = 1e-10)
})
