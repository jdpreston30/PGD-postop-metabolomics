#* 1: Cleanup
#+ Import Data
  #- Path set Amshu
    raw_path <- "C://Users//amshu//OneDrive - Emory//Preston, Joshua's files - Amshu Josh PGD//Old Analysis"
  #- Path set Josh
    raw_path <- "/Users/jdp2019/Library/CloudStorage/OneDrive-Emory/Research/Manuscripts and Projects/Active Projects/Chan Lab/TPMO/Metabolomics - Postop PGD_Amshu/Old Analysis"
  #- Import Metadata
    metadata <- read_csv(paste0(raw_path, "/metadata.csv"))
  #- Import Untargeted
    Sys.setenv("VROOM_CONNECTION_SIZE" = 2^22)
    C18_unt <- read_csv(paste0(raw_path, "/C18_untargeted.csv"))
    HILIC_unt <- read_csv(paste0(raw_path, "/HILIC_untargeted.csv"))
  #- Import targeted_metabolite_comparison
    C18_targeted_raw <- read_csv(paste0(raw_path, "/C18_targeted.csv"))
    HILIC_targeted_raw <- read_csv(paste0(raw_path, "/HILIC_targeted.csv"))
#+ Skip and import objects
  #- Path set Amshu
    all_objects_path <- "C://Users//amshu//OneDrive - Emory//Preston, Joshua's files - Amshu Josh PGD//all_objects.rds"
  #- Path set Josh
    all_objects_path <- "/Users/jdp2019/Library/CloudStorage/OneDrive-Emory/Research/Manuscripts and Projects/Active Projects/Chan Lab/TPMO/Metabolomics - Postop PGD_Amshu/all_objects.rds"
  #- Load and list RDS
    all_objects <- readRDS(all_objects_path)
    list2env(all_objects, envir = .GlobalEnv)
#+ Set up untargeted data
  #- Preprocess and join
    C18_merged_data <- metadata %>%
      left_join(C18_unt, by = c("Sample_ID" = "Sample_Name"))
    HILIC_merged_data <- metadata %>%
      left_join(HILIC_unt, by = c("Sample_ID" = "Sample_Name"))
    common_cols <- setdiff(intersect(names(HILIC_merged_data), names(C18_merged_data)), "Sample_ID")
    combined_unt <- HILIC_merged_data %>%
      left_join(C18_merged_data %>% select(-all_of(common_cols)), by = "Sample_ID")
  #- Remove any columns with > 20% missing values
    feature_cols <- names(combined_unt)[!names(combined_unt) %in% c(common_cols, "Sample_ID")]
    zero_percentages <- combined_unt %>%
      select(all_of(feature_cols)) %>%
      summarise(across(everything(), ~ sum(.x == 0, na.rm = TRUE) / length(.x))) %>%
      pivot_longer(everything(), names_to = "feature", values_to = "zero_pct")
  #- Identify features with <= 20% zeros (keep these)
    features_to_keep <- zero_percentages %>%
      filter(zero_pct <= 0.20) %>%
      pull(feature)
  #- Filter UFT_metaboanalyst to keep only good features
    UFT_metaboanalyst_raw <- combined_unt %>%
      select(all_of(c(common_cols, "Sample_ID")), all_of(features_to_keep))
  #- Replace 0 values which remain with 1/2 the column minimum
    UFT_metaboanalyst_halfmin <- UFT_metaboanalyst_raw %>%
      mutate(across(all_of(features_to_keep), ~ ifelse(.x == 0, 0.5 * min(.x[.x > 0], na.rm = TRUE), .x)))
  #- Log2 transform dataset
    UFT_metaboanalyst_log2_i <- UFT_metaboanalyst_halfmin %>%
      mutate(across(all_of(features_to_keep), ~ log2(.x)))
#+ Set up targeted data
  #- HILIC
    C18_targeted <- C18_targeted_raw %>%
      filter(str_starts(Sample_ID, "H")) %>%
      mutate(across(where(is.numeric), ~.)) %>%
      select(Sample_ID, where(~ is.numeric(.) && mean(. == 0) <= 0.20)) %>%
      mutate(across(where(is.numeric), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .))) %>%
      mutate(across(where(is.numeric), log2)) %>%
      left_join(metadata %>% select(Sample_ID, Clinical_PGD, Time, Patient), by = "Sample_ID") %>%
      select(Patient, Sample_ID, Time, Clinical_PGD, everything()) %>%
      arrange(Clinical_PGD)
  #- C18
    HILIC_targeted <- HILIC_targeted_raw %>%
      filter(str_starts(Sample_ID, "H")) %>%
      mutate(across(where(is.numeric), ~.)) %>%
      select(Sample_ID, where(~ is.numeric(.) && mean(. == 0) <= 0.20)) %>%
      mutate(across(where(is.numeric), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .))) %>%
      # Log2 transform the data
      mutate(across(where(is.numeric), log2)) %>%
      left_join(metadata %>% select(Sample_ID, Clinical_PGD, Time, Patient), by = "Sample_ID") %>%
      select(Patient, Sample_ID, Time, Clinical_PGD, everything()) %>%
      arrange(Clinical_PGD)
  #- Read metabolite key
    feature_key <- read_csv(paste0(raw_path, "/targeted_feature_key.csv"))
    HILIC_key <- feature_key[grepl("HILIC", feature_key$Name), c(1, 3)]
    C18_key <- feature_key[grepl("C18", feature_key$Name), c(1, 3)]
  #- Pivot the data to long format for easier merging
    C18_targeted_proc <- C18_targeted %>%
      pivot_longer(cols = -c(Sample_ID, Patient, Time, Clinical_PGD), names_to = "Metabolite", values_to = "Value")
    HILIC_targeted_proc <- HILIC_targeted %>%
      pivot_longer(cols = -c(Sample_ID, Patient, Time, Clinical_PGD), names_to = "Metabolite", values_to = "Value")
  #- Merge the metabolite key with the C18 and HILIC data
    C18_targeted_named <- C18_targeted_proc %>%
      left_join(C18_key, by = c("Metabolite" = "Name")) %>%
      select(Patient, Sample_ID, Time, Clinical_PGD, everything())
    HILIC_targeted_named <- HILIC_targeted_proc %>%
      left_join(HILIC_key, by = c("Metabolite" = "Name")) %>%
      select(Patient, Sample_ID, Time, Clinical_PGD, everything())
  #- Combine C18 and HILIC targeted data
    combined_targeted <- rbind(C18_targeted_named, HILIC_targeted_named)
  #- Pull the combined name key
    combined_key <- combined_targeted %>%
      select(Metabolite, Identified_Name) %>%
      distinct()
  #- Create combined TFT, remove -12 and non-complete samples
    combined_TFT_i <- combined_targeted %>%
      select(Patient,Sample_ID,Time,Clinical_PGD,Metabolite,Value) %>%
      tidyr::pivot_wider(
        id_cols = c(Sample_ID),
        names_from = Metabolite,
        values_from = Value
      ) %>%
      left_join(metadata, by = "Sample_ID") %>%
      select(Patient, Time, Sex, Age, Clinical_PGD, Sample_ID, everything())
#+ Final Data Processing for Targeted
  #- Touch up TFT and remove NA sample (H19S2), H3, and all preop timepoints
    combined_TFT_ii <- combined_TFT_i %>%
      mutate(Clinical_PGD = as.factor(if_else(Clinical_PGD == "Y'", "Y", Clinical_PGD))) %>%
      filter(Patient != "H3") %>% #! Missing 12 and 24, removing
      filter(Time != -12) %>%
      filter(Sample_ID != "H19S2")
  #- Filter to Incomplete Cases
    incomplete_patients <- combined_TFT_ii %>%
      select(Patient, Time, Clinical_PGD) %>%
      group_by(Patient) %>%
      filter(!all(c(12, 24) %in% Time)) %>%
      ungroup()
    #! H2 is missing 12 hour (PGD Y)
    #! H5 is missing 12 hour (PGD Y)
    #! H19 is missing 24 hour (PGD Y)
  #- Define metabolite columns by exclusion (robust)
    meta_cols <- c("Patient", "Time", "Sex", "Age", "Clinical_PGD", "Sample_ID")
    met_cols <- setdiff(names(combined_TFT_ii), meta_cols)
  #- Build placeholder rows for the missing patient–time combos
    new_rows <- incomplete_patients %>%
      select(Patient, Time, Clinical_PGD) %>%
      mutate(Sex = NA_character_, Age = NA_real_) %>%
      tibble::add_column(!!!setNames(rep(list(NA_real_), length(met_cols)), met_cols)) %>%
      mutate(Time = if_else(Time == 12, 24, 12)) %>%
      mutate(
        Sample_ID = case_when(
          Patient == "H19" & Time == 24 ~ paste0(Patient, "S2"),
          Patient == "H2" & Time == 12 ~ paste0(Patient, "S1"),
          Patient == "H5" & Time == 12 ~ paste0(Patient, "S1"),
          TRUE ~ NA_character_
        )
      )
  #- Time-matched means for each metabolite (ignore PGD)
    time_means <- combined_TFT_ii %>%
      group_by(Time) %>%
      summarise(across(all_of(met_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  #- Fill metabolite NAs in new rows with the time means
    new_rows_filled <- new_rows %>%
      left_join(time_means, by = "Time", suffix = c("", "_mean")) %>%
      mutate(across(
        all_of(met_cols),
        ~ coalesce(.x, get(paste0(cur_column(), "_mean")))
      )) %>%
      select(-ends_with("_mean"))
  #- Append to the wide table
    combined_TFT <- bind_rows(combined_TFT_ii, new_rows_filled) %>%
      select(-c(Sex,Age)) %>%
      left_join(metadata %>% select(Patient, Sex, Age) %>% distinct(), by = "Patient") %>%
      arrange(Patient) %>%
      select(all_of(meta_cols), all_of(met_cols)) %>%
      mutate(
        Patient = factor(Patient),
        Time = factor(Time, levels = c(12, 24)),
        Clinical_PGD = factor(Clinical_PGD, levels = c("N", "Y")),
        Sex = factor(Sex, levels = c("M", "F"))
      )
#+ Final Data Processing for Untargeted
  #- Touch up TFT and remove NA sample (H19S2), H3, and all preop timepoints
    UFT_metaboanalyst_log2_ii <- UFT_metaboanalyst_log2_i %>%
      mutate(Clinical_PGD = as.factor(if_else(Clinical_PGD == "Y'", "Y", Clinical_PGD))) %>%
      filter(Patient != "H3") %>% # ! Missing 12 and 24, removing
      filter(Time != -12) %>%
      filter(Sample_ID != "H19S2")
  #- Filter to Incomplete Cases
    incomplete_patients_unt <- UFT_metaboanalyst_log2_ii %>%
      select(Patient, Time, Clinical_PGD) %>%
      group_by(Patient) %>%
      filter(!all(c(12, 24) %in% Time)) %>%
      ungroup()
  #- Define metabolite columns by exclusion (robust)
    met_cols_unt <- setdiff(names(UFT_metaboanalyst_log2_ii), meta_cols)
  #- Build placeholder rows for the missing patient–time combos
    new_rows_unt <- incomplete_patients_unt %>%
      select(Patient, Time, Clinical_PGD) %>%
      mutate(Sex = NA_character_, Age = NA_real_) %>%
      tibble::add_column(!!!setNames(rep(list(NA_real_), length(met_cols_unt)), met_cols_unt)) %>%
      mutate(Time = if_else(Time == 12, 24, 12)) %>%
      mutate(
        Sample_ID = case_when(
          Patient == "H19" & Time == 24 ~ paste0(Patient, "S2"),
          Patient == "H2" & Time == 12 ~ paste0(Patient, "S1"),
          Patient == "H5" & Time == 12 ~ paste0(Patient, "S1"),
          TRUE ~ NA_character_
        )
      )
  #- Time-matched means for each metabolite (ignore PGD)
  time_means_unt <- UFT_metaboanalyst_log2_ii %>%
    group_by(Time) %>%
    summarise(across(all_of(met_cols_unt), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  #- Fill metabolite NAs in new rows with the time means
    new_rows_filled_unt <- new_rows_unt %>%
      left_join(time_means_unt, by = "Time", suffix = c("", "_mean")) %>%
      mutate(across(
        all_of(met_cols_unt),
        ~ coalesce(.x, get(paste0(cur_column(), "_mean")))
      )) %>%
      select(-ends_with("_mean"))
  #- Append to the wide table
    combined_UFT <- bind_rows(UFT_metaboanalyst_log2_ii, new_rows_filled_unt) %>%
      select(-c(Sex, Age)) %>%
      left_join(metadata %>% select(Patient, Sex, Age) %>% distinct(), by = "Patient") %>%
      arrange(Patient) %>%
      select(all_of(meta_cols), all_of(met_cols_unt)) %>%
        mutate(
          Patient = factor(Patient),
          Time = factor(Time, levels = c(12, 24)),
          Clinical_PGD = factor(Clinical_PGD, levels = c("N", "Y")),
          Sex = factor(Sex, levels = c("M", "F"))
        )
