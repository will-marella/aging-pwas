# Illustration only. Review the age center and adjustment set
# before defining a cohort configuration.
spec <- pwas_time_spec(
  age_center = 50,
  covariates = list(
    SITE = list(type = "factor", timing = "baseline",
                levels = c("A", "B", "C"), reference = "A")
  ),
  random_effects = "intercept_slope",
  min_subjects = 20L,
  confidence_level = 0.95
)

preprocessing <- list(
  abundance_scale = "Synthetic continuous abundance in arbitrary units",
  normalization = "None; simulated on a shared scale across visits",
  batch_handling = "No batch effects simulated; SITE is a baseline covariate",
  missing_values = "No missing values simulated"
)
