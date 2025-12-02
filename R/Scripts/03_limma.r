#*3: LIMMA
#+ 3.1: Prepare data for LIMMA
#- 3.1.0: Define synthetic rows to exclude from LIMMA
missing_samples <- tibble::tribble(
  ~Patient, ~Sample,
  "H5",     "S2",
  "H14",    "S2",
  "H50",    "S1"
)
#- 3.1.1: Clean up metadata into factors (Untargeted)
UFT_filtered_LIMMA <- UFT_filtered %>%
  anti_join(missing_samples, by = c("Patient", "Sample")) %>%
  select(Patient, Time, severe_PGD, everything(),-c(Sample_ID, Sample, PGD_grade_tier, any_PGD)) %>%
  mutate(
    Time = factor(Time, levels = c(12, 24), labels = c("12", "24")),
    Clinical_PGD = factor(severe_PGD, levels = c("N", "Y")),
    Patient = factor(Patient)
  )
#- 3.1.2: Clean up metadata into factors (Targeted)
TFT_combined_LIMMA <- TFT_combined %>%
  anti_join(missing_samples, by = c("Patient", "Sample")) %>%
  select(Patient, Time, severe_PGD, everything(), -c(Sample_ID, Sample, PGD_grade_tier, any_PGD)) %>%
  mutate(
    Time = factor(Time, levels = c(12, 24), labels = c("12", "24")),
    Clinical_PGD = factor(severe_PGD, levels = c("N", "Y")),
    Patient = factor(Patient)
  )
#+ 3.2: Run LIMMA
#- 3.2.1: Run on untargeted
limma_untarg <- run_limma_three_contrasts(
  data            = UFT_filtered_LIMMA,
  time_var        = "Time",
  group_var       = "severe_PGD",
  block_var       = "Patient",
  time_levels     = c("12", "24"),
  group_levels    = c("N", "Y"),
  meta_cols       = c("Patient", "Time", "severe_PGD"),
  metaboanalyst   = TRUE, # <- emit mummichog p-lists
  targeted_rename = FALSE # <- do NOT join/rename with any key
)
#- 3.2.2: Run on targeted
limma_targ <- run_limma_three_contrasts(
  data            = TFT_combined_LIMMA,
  time_var        = "Time",
  group_var       = "severe_PGD",
  block_var       = "Patient",
  time_levels     = c("12", "24"),
  group_levels    = c("N", "Y"),
  meta_cols       = c("Patient", "Time", "severe_PGD"), 
  metaboanalyst   = FALSE, # also emit mummichog p-lists
  targeted_rename = FALSE, # append Identified_Name (etc.)
  rename_key      = FALSE # must have a "Metabolite" column
)