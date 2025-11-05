#*3: LIMMA
#+ 3.1: Prepare data for LIMMA
  #- 3.1.1: Drop fully missing rows... keep H19 which is half missing (no C18) (Untargeted)
    UFT_LIMMA <- UFT_metaboanalyst %>%
      filter(Sample_ID != "H2_12" & Sample_ID != "H5_12") %>%
      select(-c(Sex, Age, Sample_ID)) %>%
      mutate(
        Time = factor(Time, levels = c(12, 24), labels = c("12", "24")),
        Clinical_PGD = factor(Clinical_PGD, levels = c("N", "Y")),
        Patient = factor(Patient)
      )  
  #- 3.1.2: Drop fully missing rows... keep H19 which is half missing (no C18) (Targeted)
    TFT_LIMMA <- TFT %>%
      filter(Sample_ID != "H2_12" & Sample_ID != "H5_12") %>%
      select(-c(Sex, Age, Sample_ID)) %>%
      mutate(
        Time = factor(Time, levels = c(12, 24), labels = c("12", "24")),
        Clinical_PGD = factor(Clinical_PGD, levels = c("N", "Y")),
        Patient = factor(Patient)
      )
#+ 3.2: Run LIMMA
  #- 3.2.1: Run on untargeted
    limma_untarg <- run_limma_three_contrasts(
      data            = UFT_LIMMA,
      time_var        = "Time",
      group_var       = "Clinical_PGD",
      block_var       = "Patient",
      time_levels     = c("12", "24"),
      group_levels    = c("N", "Y"),
      meta_cols       = c("Patient", "Time", "Clinical_PGD"),
      export_prefix   = "UFT_limma",
      export_dir      = "Outputs/LIMMA",
      metaboanalyst   = TRUE, # <- emit mummichog p-lists
      targeted_rename = FALSE # <- do NOT join/rename with any key
    )
  #- 3.2.2: Run on targeted
    limma_targ <- run_limma_three_contrasts(
      data            = TFT_LIMMA,
      time_var        = "Time",
      group_var       = "Clinical_PGD",
      block_var       = "Patient",
      time_levels     = c("12", "24"),
      group_levels    = c("N", "Y"),
      meta_cols       = c("Patient", "Time", "Clinical_PGD"), 
      metaboanalyst   = FALSE, # also emit mummichog p-lists
      targeted_rename = TRUE, # append Identified_Name (etc.)
      rename_key      = combined_key_clean # must have a "Metabolite" column
    )


