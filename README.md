# PWAS of the passage of time

A small R pipeline for testing longitudinal protein abundance trajectories. CARDIA and MESA use the same analysis engine and versioned model specification; cohort preparation runs separately in the secure computing environment.

The draft model fits one mixed model per protein:

```r
protein ~ TIME_YEARS + I(TIME_YEARS^2) +
  AGE_C + I(AGE_C^2) + TIME_YEARS:AGE_C +
  FEMALE + additional_covariates +
  (1 + TIME_YEARS | SUBJECT_ID)
```

`AGE_C` is baseline age minus a prespecified constant. Population curvature, differences in slopes by baseline age, and individual variation in slopes are distinct model components. This is a longitudinal association analysis; it does not construct a shared aging score across proteins.

## Run the synthetic example

Install the R packages `lme4` and `lmerTest` in your chosen environment, then run:

```sh
Rscript examples/run_synthetic.R
```

This generates six synthetic proteins measured at four visits, fits the models, and writes results and synthetic inputs under the ignored `outputs/synthetic/` directory. Existing output directories are refused. Supply a new directory and explicit core count if needed:

```sh
Rscript examples/run_synthetic.R outputs/synthetic-review-2 2
```

The example's age center of 50 and its covariate choices are illustrative. Scientific defaults remain subject to review before cohort analysis.

## Use the R interface

```r
source("R/pwas_time.R")

spec <- pwas_time_spec(
  age_center = 50,
  covariates = list(
    SITE = list(type = "factor", timing = "baseline",
                levels = c("A", "B", "C"), reference = "A")
  )
)

preprocessing <- list(
  abundance_scale = "Describe the supplied abundance scale",
  normalization = "Describe upstream normalization",
  batch_handling = "Describe upstream batch handling",
  missing_values = "Describe upstream missing-value handling"
)

result <- run_pwas_time(pheno, omics, spec, preprocessing, n_cores = 1L)
summarize_pwas_time(result)
write_pwas_time(result, "outputs/cohort-run")
```

## Command-line runner

Prepare `pheno` and `omics` as RDS data frames and a private R configuration defining `spec` and `preprocessing`. Run:

```sh
Rscript scripts/run_pwas.R PHENO.rds OMICS.rds CONFIG.R OUTPUT_DIR N_CORES
```

The command resolves repository code independently of the working directory. Input and output arguments are resolved from the working directory. See `INPUTS_OUTPUTS.md` for the input contract, coefficient inference, missingness rules, and output interpretation.

Only code, documentation, and synthetic examples belong in this repository. Cohort inputs and results remain on HPC. The ignore rules are a convenience; inspect changes before staging files.

## Outputs

Each run writes two result files:

- `result.rds`: coefficients, coefficient covariance matrices, QC, and run metadata.
- `results.csv`: one row per protein with all fixed-effect estimates, SEs, confidence intervals, p-values/BH corrections, counts, and fit status.

Inference is at the coefficient level. Reduced-model tests and generated change contrasts are not run. The retained covariance matrices support later combinations of coefficients and their standard errors. See `INPUTS_OUTPUTS.md` for the exact schema and limitations.

## Validate the draft

Install `testthat` in addition to the analysis dependencies, then run:

```sh
Rscript tests/run_tests.R
```

Tests compare estimates, covariance, and coefficient inference against direct R fits; verify one fit per eligible protein, known synthetic trajectories, sample ordering, missingness, failed fits, and serial/parallel agreement; and check the two-file export and output preservation. This is implementation validation, not a power study or proof that either cohort's preprocessing supports longitudinal inference.

Development checks used R 4.3.1, lme4 1.1.35.1, and lmerTest 3.1.3. Every run records its actual package versions and code hashes. No package lockfile is supplied in this draft.
