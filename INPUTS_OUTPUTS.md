# Input and output contract

This engine consumes analysis-ready abundance values. It does not read cohort-specific releases, transform abundance values, standardize within visits, impute missing observations, adjust batches, or infer covariates from column names. The same specification should run in discovery and replication.

## Input tables

`pheno` is a data frame with one row per sample:

| Column | Contract |
| --- | --- |
| `SAMPLE_ID` | Nonmissing, nonempty, unique character sample identifier. |
| `SUBJECT_ID` | Nonmissing, nonempty character participant identifier. |
| `TIME_YEARS` | Finite numeric elapsed years from the defined proteomic baseline; nonnegative. Each input participant has a baseline observation at zero. |
| `BASELINE_AGE` | Finite numeric age from 0 to 120 years at that baseline; constant within participant. |
| `FEMALE` | Numeric 0/1 or factor labels `"0"`/`"1"`; observed values constant within participant. Missing values follow covariate exclusions. Internally a factor with reference `"0"`. |
| Declared covariates | Named and typed in the specification. Missing covariate values are handled by the documented exclusion policy. |

Each participant–time pair must be unique. The time origin and baseline age are fixed before analyte-specific filtering: losing a baseline abundance must not reset time or age.

`omics` is a data frame with one row per analyte. The first column is `ANALYTE_NAME`, containing nonmissing, nonempty, unique character stable assay identifiers. Every remaining column is a numeric abundance vector named exactly as one `pheno$SAMPLE_ID`. Both tables must contain exactly the same sample set; input order may differ. Alignment uses names, never column position or a silent intersection. Preserve sample names when reading tables, for example with `check.names = FALSE`.

Keep assay identifiers and explicit protein/platform mappings upstream. Identical gene labels do not establish identical assay targets or comparable abundance scales across cohorts.

## Specification

Construct the versioned specification with:

```r
pwas_time_spec(
  age_center,
  covariates = list(),
  random_effects = "intercept_slope",
  min_subjects = 20L,
  confidence_level = 0.95
)
```

The fixed-effects model is:

```r
RESPONSE ~ TIME_YEARS + I(TIME_YEARS^2) +
  AGE_C + I(AGE_C^2) + TIME_YEARS:AGE_C +
  FEMALE + additional_covariates
```

`AGE_C = BASELINE_AGE - age_center`. Use a shared, prespecified age center rather than a cohort mean. Additional covariates enter as main effects; arbitrary caller-supplied formulas are not supported.

Each covariate is a named list with `type` (`"numeric"`, `"factor"`, or `"logical"`) and `timing` (`"baseline"` or `"time_varying"`). Factors also declare `levels` and `reference`. For example:

```r
covariates = list(
  SITE = list(type = "factor", timing = "baseline",
              levels = c("A", "B", "C"), reference = "A"),
  BASELINE_BMI = list(type = "numeric", timing = "baseline")
)
```

Baseline covariates must be consistent within participants. Logical covariates must be logical vectors and are modeled with explicit `FALSE = 0` and `TRUE = 1` coding. Factor covariates accept factors or character vectors containing the declared labels and use treatment coding with the declared reference. Time-varying adjustment changes the scientific interpretation: estimated time effects then condition on those variables. Covariate selection and timing require substantive justification.

The default random structure is `(1 + TIME_YEARS | SUBJECT_ID)`, with correlated participant intercepts and linear slopes. `random_effects = "intercept"` prespecifies `(1 | SUBJECT_ID)`. The engine does not switch structures after a failed or singular fit. Population curvature does not require participant-specific random quadratic terms.

The default minimum of 20 retained participants is a draft operational threshold, not a power calculation or a guarantee of reliable estimation.

Covariate names cannot be formula operators or reserved fields. Factor coding must produce unique coefficient names: for example, a factor `SITE` with level `B` cannot coexist with a numeric covariate named `SITEB`. Only declared modeling columns enter the fits.

### Preprocessing declaration

`preprocessing` is a named list of nonempty strings with these required fields:

```r
list(
  abundance_scale = "...",
  normalization = "...",
  batch_handling = "...",
  missing_values = "..."
)
```

Describe what was actually done before calling the engine. Processing must preserve the longitudinal changes being tested; a declaration does not establish measurement comparability. In particular, visit-specific normalization can remove population shifts, and batch perfectly confounded with visit cannot be separated using these measurements alone.

## Missingness and model eligibility

For each analyte, the engine excludes observations with missing abundance or required covariate values, then excludes participants with fewer than two retained observations. It reports aggregate exclusions and retained participant/sample counts. It does not impute observations or reset participant baselines after exclusions.

Models are screened for constant response, insufficient participants, and fixed-effect rank deficiency. Fits are assessed for convergence and singular random-effects structure. Failed models remain in the outputs with explicit status; they do not disappear from the analyte inventory. Numeric infinities are invalid inputs, not ordinary missing observations.

This available-observation policy does not resolve informative dropout, survivor selection, or missing-not-at-random abundances. Interpretation applies to participants contributing eligible repeated measurements under the modeled assumptions.

The random-intercept-and-slope specification also requires more retained observations than its two random effects per participant. If every retained participant has only two observations, that structure is flagged as unsupported; it is not replaced automatically.

At least three distinct times are needed to estimate separate linear and quadratic time effects under a common visit schedule. The fixed-effect rank check detects exact aliasing. Three or more times do not alone establish adequate precision or individual support for curvature.

## Inference

Each eligible protein gets one maximum-likelihood mixed-model fit. The engine reports fixed-effect estimates, standard errors, Satterthwaite degrees of freedom, confidence intervals, and p-values. BH correction runs across valid proteins separately for each coefficient term. There are no reduced-model comparisons or generated change contrasts.

`TIME_YEARS` describes the slope at time zero and the chosen age center; `I(TIME_YEARS^2)` describes curvature; `TIME_YEARS:AGE_C` describes how slopes differ by baseline age. Read estimates and uncertainty together. Coefficient magnitudes have different units, and p-values do not measure effect importance. Individual term tests do not provide an omnibus test of the full trajectory.

Converged singular fits retain computed confidence intervals and p-values, including BH adjustment, with `STATUS = "singular"` and `FULL_SINGULAR = TRUE`. `INFERENCE_OK` indicates that coefficient inference was computed successfully; it does not clear the singularity flag. The engine does not refit these proteins with a different random-effects structure.

Other failed fits remain in the results with explicit status. Available estimates, standard errors, and covariance may be retained while confidence intervals and p-values are withheld. Use `INFERENCE_OK`, `STATUS`, and `FULL_SINGULAR` when interpreting results; missing inference is not evidence of no association.

## Result object and files

`run_pwas_time(pheno, omics, spec, preprocessing, n_cores = 1L, verbose = TRUE, output_dir = NULL, checkpoint_every = 50L)` returns one result object. Supplying `output_dir` also enables incremental saving and automatic resume:

| Component | Contents |
| --- | --- |
| `coefficients` | One row per protein and fixed-effect term: `ANALYTE_NAME`, `TERM`, `ESTIMATE`, `SE`, `DF`, `CI_LOW`, `CI_HIGH`, `P_VALUE`, `BH_P_VALUE`, `INFERENCE_OK`, `STATUS`. |
| `covariance` | A named list of fixed-effect coefficient covariance matrices, one per protein when available. |
| `model_qc` | One row per protein: retained sample/participant counts, age/time ranges, rank, convergence, singularity, residual SD, warnings/errors, status, and worker PID. |
| `visit_coverage` | Retained visit-count distributions per protein. |
| `exclusions` | Sequential, nonoverlapping sample exclusion counts per protein and reason. |
| `multiplicity` | Requested, valid, and withheld coefficient p-value counts per term. |
| `metadata` | Model specification, preprocessing, input QC, execution settings, source/package versions, checksums, and session information. |

`TERM` holds the exact model coefficient name, such as `TIME_YEARS`, `I(TIME_YEARS^2)`, `TIME_YEARS:AGE_C`, or `FEMALE1`. `ANALYTE_NAME` identifies the protein wherever a protein key applies.

### Coefficient covariance

For one protein, the saved matrix is `Cov(beta_hat)`: its row and column names identify the fixed-effect coefficients. Diagonal entries are squared standard errors; off-diagonal entries describe covariance between coefficient estimates. This is neither covariance between proteins nor the random-intercept/slope variance matrix.

For a future combination with weights `L`, the coefficients give `sum(L * beta_hat)` and the covariance gives `SE = sqrt(t(L) %*% V %*% L)`. Exact Satterthwaite degrees of freedom for a new combination require additional fitted-model information, which is not saved. No change contrasts are computed or stored by this version.

### Saved files

`write_pwas_time(result, output_dir)` writes exactly two result files:

| File | Contents |
| --- | --- |
| `result.rds` | The complete result object described above. |
| `results.csv` | One row per protein, with fit status, counts, and all fixed-effect coefficient summaries. |

The CSV starts with `ANALYTE_NAME`, `STATUS`, `N_OBS`, `N_SUBJECTS`, `FULL_CONVERGED`, and `FULL_SINGULAR`. It then has `<TERM>__<FIELD>` columns for `ESTIMATE`, `SE`, `DF`, `CI_LOW`, `CI_HIGH`, `P_VALUE`, `BH_P_VALUE`, and `INFERENCE_OK`. For example, `TIME_YEARS__ESTIMATE` is the linear time estimate and `I(TIME_YEARS^2)__BH_P_VALUE` is the adjusted curvature coefficient p-value. Read with `check.names = FALSE` in R to preserve those names. Failed proteins remain in the CSV; missing values are empty cells.

The RDS retains detailed QC and provenance without exporting extra files or storing another copy of the wide CSV table. It contains no input tables, sample/participant identifiers, individual random effects, or fitted-model objects. Declared covariate labels and preprocessing descriptions remain in metadata and must not contain participant identifiers.

`summarize_pwas_time(result)` returns an aggregate console summary, including checkpoint progress when enabled. Its exclusion totals sum across proteins, not unique participants or samples. In a partial checkpoint, model status, exclusions, and multiplicity counts describe completed proteins; input QC describes the entire requested input. It does not write a separate summary file.

The standalone writer accepts a new or empty directory and refuses a nonempty directory. The generic `scripts/run_pwas.R` CLI and synthetic runner use this final-export mode. The synthetic runner additionally saves its generated `synthetic_pheno.rds`, `synthetic_omics.rds`, and `synthetic_truth.csv` for example reuse.

With `output_dir` supplied to `run_pwas_time`, results are saved every `checkpoint_every` proteins and after the final batch. The CARDIA runner enables this by default. Only the parent process writes, replacing each file through a temporary file in the output folder and committing `result.rds` last. This RDS is the authoritative checkpoint; rerunning the same analysis resumes after its completed proteins, including flagged/failed fits. An interruption may require refitting the current unsaved batch. Do not run two processes against the same output folder.

Resume requires an identical hash of the validated inputs, model specification, preprocessing declarations, engine source hashes, and model-package versions. Changing worker count or checkpoint interval is allowed. Old exports without checkpoint metadata and mismatched runs are refused. Data preparation and hashing run again on restart; no extra copy of the inputs is saved. `metadata$checkpoint` records `n_completed`, `n_total`, `complete`, and the signature. A completed checkpoint is returned without fitting again.

Partial files contain completed proteins and available raw p-values/CIs. `BH_P_VALUE` stays missing until all requested proteins have finished; final BH correction includes every valid p-value in the requested analysis. There are no extra checkpoint files beyond `result.rds` and `results.csv`. When checkpointing is enabled, the runner already saves the outputs, so do not call the standalone writer afterward.

`n_cores` is explicit and defaults to one. Parallel execution uses Unix fork workers. With `output_dir = NULL`, results are collected in memory and returned after fitting, as before.

Version 0.2.0 removes the earlier draft's likelihood-ratio tests, contrast configuration, and individual table exports. Version 0.2.1 retains inference for converged singular fits and flags them. Recreate configurations with the current `pwas_time_spec()` interface and rerun fits to obtain previously withheld inference; existing output directories are not modified.

Pipeline version 0.3.0 adds checkpoint/resume; the model specification remains version 0.2.1.

## Command-line runner

```sh
Rscript scripts/run_pwas.R PHENO.rds OMICS.rds CONFIG.R OUTPUT_DIR N_CORES
```

The configuration is executable R code loaded into a separate environment. It must define `spec` and `preprocessing`; `pwas_time_spec` is made available there. Use a trusted configuration. The CLI adds the two input-file MD5 checksums as `metadata$input_md5` and the configuration-file checksum as `metadata$config_md5`, without recording their filenames. The engine records the specification hash, source hashes, available git state, package versions, and session information. Retain this metadata with each output directory and retain the corresponding git commit in your execution records.

`specs/synthetic.R` demonstrates a configuration. Its values describe synthetic data and must be replaced with a reviewed cohort specification and accurate preprocessing declaration before a cohort run. Keep private configurations and cohort artifacts outside the repository, or in ignored private directories, and inspect the git diff before sharing code.
