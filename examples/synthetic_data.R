# Synthetic inputs only; no cohort files are required.
make_synthetic_pwas <- function(n_subjects = 160L, seed = 20260905L) {
  if (length(n_subjects) != 1L || !is.numeric(n_subjects) ||
      !is.finite(n_subjects) || n_subjects < 4L ||
      n_subjects > .Machine$integer.max ||
      n_subjects != as.integer(n_subjects)) {
    stop("n_subjects must be an integer of at least four.", call. = FALSE)
  }
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("seed must be a nonnegative integer.", call. = FALSE)
  }
  set.seed(seed)
  n_subjects <- as.integer(n_subjects)
  times <- c(0, 5, 10, 15)
  subject_index <- rep(seq_len(n_subjects), each = length(times))
  baseline_age <- runif(n_subjects, 35, 65)
  female <- rep(c(0L, 1L), length.out = n_subjects)
  site <- rep(c("A", "B", "C"), length.out = n_subjects)

  pheno <- data.frame(
    SAMPLE_ID = sprintf("synthetic_sample_%04d", seq_along(subject_index)),
    SUBJECT_ID = sprintf("synthetic_person_%04d", subject_index),
    TIME_YEARS = rep(times, times = n_subjects),
    BASELINE_AGE = baseline_age[subject_index],
    FEMALE = female[subject_index],
    SITE = factor(site[subject_index], levels = c("A", "B", "C")),
    stringsAsFactors = FALSE
  )

  effects <- data.frame(
    ANALYTE_NAME = c("synthetic_null", "synthetic_linear_positive",
                     "synthetic_linear_negative", "synthetic_curvature",
                     "synthetic_age_interaction", "synthetic_combined"),
    intercept = c(2, 3, 4, 5, 6, 7),
    linear = c(0, 0.08, -0.08, -0.07, 0, 0.03),
    quadratic = c(0, 0, 0, 0.009, 0, 0.006),
    interaction = c(0, 0, 0, 0, 0.007, -0.005),
    stringsAsFactors = FALSE
  )
  age_centered <- pheno$BASELINE_AGE - 50
  baseline_terms <- 0.04 * age_centered + 0.001 * age_centered^2 +
    0.15 * pheno$FEMALE + 0.12 * (pheno$SITE == "B") -
    0.10 * (pheno$SITE == "C")
  abundance <- matrix(NA_real_, nrow = nrow(effects), ncol = nrow(pheno))
  truth_rows <- vector("list", nrow(effects))

  for (j in seq_len(nrow(effects))) {
    # Draw distinct participant effects for every protein.
    random_intercept <- rnorm(n_subjects, sd = 1)
    random_slope <- rnorm(n_subjects, sd = 0.05)
    abundance[j, ] <- effects$intercept[j] + baseline_terms +
      effects$linear[j] * pheno$TIME_YEARS +
      effects$quadratic[j] * pheno$TIME_YEARS^2 +
      effects$interaction[j] * pheno$TIME_YEARS * age_centered +
      random_intercept[subject_index] +
      random_slope[subject_index] * pheno$TIME_YEARS +
      rnorm(nrow(pheno), sd = 0.2)
    truth_rows[[j]] <- data.frame(
      ANALYTE_NAME = effects$ANALYTE_NAME[j],
      TERM = c("(Intercept)", "TIME_YEARS", "I(TIME_YEARS^2)",
               "AGE_C", "I(AGE_C^2)", "FEMALE1", "SITEB", "SITEC",
               "TIME_YEARS:AGE_C"),
      ESTIMATE = c(effects$intercept[j], effects$linear[j],
                   effects$quadratic[j], 0.04, 0.001, 0.15, 0.12, -0.10,
                   effects$interaction[j]),
      stringsAsFactors = FALSE
    )
  }
  colnames(abundance) <- pheno$SAMPLE_ID
  omics <- data.frame(ANALYTE_NAME = effects$ANALYTE_NAME,
                      abundance, check.names = FALSE, stringsAsFactors = FALSE)
  truth <- do.call(rbind, truth_rows)
  rownames(truth) <- NULL
  preprocessing <- list(
    abundance_scale = "Synthetic continuous abundance in arbitrary units",
    normalization = "None; simulated on a shared scale across visits",
    batch_handling = "No batch effects simulated; SITE is a baseline covariate",
    missing_values = "No missing values simulated"
  )
  list(pheno = pheno, omics = omics, truth = truth,
       preprocessing = preprocessing)
}
