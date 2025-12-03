#* 0d: Importing Clinical Metadata and Preprocess
#+ 0d.0: Set Vectors from Config
#- 0d.0.1: Generate patient list from config
patient_range <- config$analysis$cohort$patient_range
all_patients <- seq(patient_range$start, patient_range$end)
excluded_patients <- patient_range$excluded
#!!!!!!!!!!1
included_patients <- all_patients[!all_patients %in% excluded_patients]
analyzed_patients <- paste0("H", included_patients)
#- 0d.0.2: Extract table variables from config
table_vars <- config$analysis$cohort$table_variables
T1 <- table_vars$T1
T2 <- table_vars$T2
T3 <- table_vars$T3
#+ 0d.1: Import Clinical and Other metadata 
#- 0d.1.1: Clinical Data 
preop_i <- read_clinical_sheet("Preop", analyzed_patients, suppress_warnings = TRUE)
periop_i <- read_clinical_sheet("Periop", analyzed_patients)
outcomes_i <- read_clinical_sheet("Outcomes", analyzed_patients)
match_run_i <- read_clinical_sheet("Match Run", analyzed_patients)
#+ 0d.2: Combine clinical metadata into one tibble; format variables
clinical_metadata_i <- preop_i |>
  left_join(periop_i, by = "Patient") |>
  left_join(outcomes_i, by = "Patient") |>
  left_join(match_run_i, by = "Patient") |>
  select(-ends_with(".x"), -ends_with(".y")) |>
  mutate(across(where(is.character), as.factor))
#+ 0d.4: Compute RADIAL score, PHM, and ISHLT PGD Status
clinical_metadata <- clinical_metadata_i |>
  calc_radial(
    rap_col = "preop_RAP",
    age_col = "demographics_age_tpx",
    dm_col = "comorbidities_DM",
    inotrope_col = "rx_preop_inotrope",
    donor_age_col = "donor_age",
    ischemic_time_col = "operative_IT_minutes"
  ) |>
  calc_BiVAD() |>
  calc_ISHLT() |>
  calc_PHM() |>
  mutate(severe_PGD = if_else(
    postop_PGD_grade_ISHLT == "Severe", "Severe PGD", "No Severe PGD",
  )) |>
  mutate(PGD_grade_tier = case_when(
    postop_PGD_grade_ISHLT %in% c("Moderate", "Mild") ~ "Mild/Mod. PGD",
    postop_PGD_grade_ISHLT == "Severe" ~ "Severe PGD",
    postop_PGD_grade_ISHLT == "N" ~ "No PGD",
    TRUE ~ NA_character_
  ))
#+ 0d.4: Break into components for the final tables 
T1_data <- clinical_metadata |>
  select(Patient, severe_PGD, all_of(T1))
T2_data <- clinical_metadata |>
  select(Patient, severe_PGD, all_of(T2))
T3_data <- clinical_metadata |>
  select(Patient, severe_PGD, all_of(T3))
#+ 0d.5: Store PGD specifics for later joining 
PGD_specifics <- clinical_metadata |>
  select(
    Patient, postop_PGD_grade_ISHLT, severe_PGD, 
    PGD_grade_tier, postop_PGD_binary_ISHLT
  ) |>
  mutate(any_PGD = as.factor(if_else(
    postop_PGD_grade_ISHLT %in% c("Mild", "Moderate", "Severe"), "Y", "N", missing = "N"
  )))

