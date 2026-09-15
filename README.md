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

## CARDIA analysis-table runner

`scripts/run_cardia.R` reads the final CARDIA analysis CSV, prepares `pheno`/`omics`, and calls the shared engine. It also requires `dplyr`, `readr`, and `tidyr`. Run from the analysis working directory where `../aging-pwas` is the clone:

```sh
Rscript ../aging-pwas/scripts/run_cardia.R
```

It defaults to 20 proteins, two workers, and checkpointing every 50 proteins (also at the end of a shorter run). For the full run, set `n_proteins <- NULL`, choose a new output directory, and set the allocated core count. It uses each participant's first `VISIT_AGE_CALC` as baseline, `FEMALE = SEX - 1`, categorical race with reference code 5, and age center 50. The engine validates the prepared inputs. Only `result.rds` and `results.csv` are written. The script can also be sourced in R from the same working directory.

Repeat samples at the same participant-time are averaged per protein on the supplied NPX scale, retaining the first sample ID. Their baseline age, sex, and race must agree. Means use available values; proteins missing in all repeats remain `NA`.

### Checkpoints and tmux

The CARDIA runner updates `result.rds` and `results.csv` after each completed batch. Rerun the same script with the same output folder to resume: preparation runs again, but saved proteins are not refitted. The checkpoint must match the prepared inputs, model, preprocessing, engine code, and model-package versions. Core count and checkpoint interval may change. Completed fits include proteins flagged as singular or failed. Use one running process per output folder.

Each file is replaced through a temporary file; `result.rds` is the authoritative checkpoint. An interruption can lose the current unsaved batch, at most 50 proteins with the default setting. Partial CSVs contain completed proteins with raw p-values/CIs; BH values remain blank until all proteins finish. Console messages show `Saved X/Y analytes`. Earlier outputs without checkpoint metadata require a new folder.

From your analysis working directory, make a local run script:

```sh
cp ../aging-pwas/scripts/run_cardia.R run_cardia_full.R
```

Edit its settings to `n_proteins <- NULL`, your allocated `n_cores`, and `output_dir <- "../Results/PWAS_Time_full"`. Then start a persistent terminal and run the script inside it:

```sh
tmux new -s pwas
Rscript run_cardia_full.R
```

Detach with **Ctrl-b**, then **d**. Reconnect using `tmux attach -t pwas`. To stop the R process, attach and press **Ctrl-c**; restart with the same `Rscript` command to resume from the last checkpoint. tmux keeps the terminal session running through an SSH disconnect; checkpointing handles an interrupted R process or machine restart.

For a console log as well, replace the `Rscript` line with `Rscript run_cardia_full.R >> ../Results/PWAS_Time_full.log 2>&1` and monitor it from another terminal with `tail -f ../Results/PWAS_Time_full.log`.

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
