library(dplyr)
library(readr)
library(tidyr)

source("../aging-pwas/R/pwas_time.R")

# Start small. For the full run, set n_proteins <- NULL and choose a new folder.
n_proteins <- 20L
n_cores <- 2L
output_dir <- "../Results/PWAS_Time_smoke"
stopifnot(!file.exists(output_dir))

proteomics <- read_csv(
  "../CARDIA_Proteomics/WTM_CARDIA_Proteomics_analysis_dataset.csv",
  show_col_types = FALSE
)

# One row per sample; baseline is the person's first available proteomic visit.
pheno <- proteomics %>%
  distinct(SAMPLE_ID, SUBJECT_ID, VISIT_AGE_CALC, SEX, RACE) %>%
  mutate(SUBJECT_ID = as.character(SUBJECT_ID),
         FEMALE = SEX - 1, RACE = as.character(RACE)) %>%
  group_by(SUBJECT_ID) %>%
  mutate(BASELINE_AGE = min(VISIT_AGE_CALC),
         TIME_YEARS = VISIT_AGE_CALC - BASELINE_AGE) %>%
  ungroup() %>%
  select(SAMPLE_ID, SUBJECT_ID, TIME_YEARS, BASELINE_AGE, FEMALE, RACE) %>%
  as.data.frame()

# One row per protein, one column per sample. Preserve missing measurements.
stopifnot(!anyDuplicated(proteomics[c("SAMPLE_ID", "OlinkID")]))
omics <- proteomics %>%
  select(ANALYTE_NAME = OlinkID, SAMPLE_ID, PCNormalizedNPX) %>%
  pivot_wider(names_from = SAMPLE_ID, values_from = PCNormalizedNPX,
              values_fill = NA_real_) %>%
  as.data.frame()

# Average repeat samples at the same participant-time, keeping the first ID.
stopifnot(!anyDuplicated(pheno$SAMPLE_ID))
repeat_visits <- pheno %>%
  group_by(SUBJECT_ID, TIME_YEARS) %>%
  filter(n() > 1) %>%
  group_split()

for (visit in repeat_visits) {
  stopifnot(nrow(unique(visit[c("BASELINE_AGE", "FEMALE", "RACE")])) == 1L)
  ids <- visit$SAMPLE_ID
  averaged <- rowMeans(omics[ids], na.rm = TRUE)
  averaged[is.nan(averaged)] <- NA_real_
  omics[[ids[1]]] <- averaged
  omics[ids[-1]] <- NULL
  pheno <- pheno[!pheno$SAMPLE_ID %in% ids[-1], ]
}

# Select proteins only after establishing every participant's baseline.
if (!is.null(n_proteins)) omics <- head(omics, n_proteins)
rm(proteomics)

spec <- pwas_time_spec(
  age_center = 50,
  covariates = list(
    RACE = list(type = "factor", timing = "baseline",
                levels = c("5", "4"), reference = "5")
  )
)

# Add the upstream normalization/QC details when available.
preprocessing <- list(
  abundance_scale = "Supplied PCNormalizedNPX",
  normalization = "Used as supplied; no additional normalization",
  batch_handling = "No additional batch correction or PlateID adjustment",
  replicate_handling = "Arithmetic mean of PCNormalizedNPX within participant-time and protein",
  missing_values = "Replicate means use available values; all-missing means remain NA; no imputation"
)

result <- run_pwas_time(pheno, omics, spec, preprocessing, n_cores = n_cores)
print(summarize_pwas_time(result))
write_pwas_time(result, output_dir)
