#* 1: Cleanup
#+ 1.1 Import Raw Data
  UFT_C18 <- read_csv("Raw_Data/C18neg_UFT_medsum.csv")
  UFT_HILIC <- read_csv("Raw_Data/HILICpos_UFT_medsum.csv")
  C18map <- read_csv("Raw_Data/Thyroid_mapping_c18neg.csv")
  HILICmap <- read_csv("Raw_Data/Thyroid_mapping_hilicpos.csv")
  tumor_IDs <- read_csv("Raw_Data/tumor_IDs.csv")
  TFT_C18 <- read_csv("Raw_Data/C18neg_TFT.csv")
  TFT_HILIC <- read_csv("Raw_Data/HILICpos_TFT.csv")  
#+ 1.2 Untargeted Processing
  #- 1.2.1: Name, make features, transpose for C18
    C18_U <- process_feature_table(UFT_C18, C18map,"C18") %>%
      left_join(tumor_IDs, by = "Sample_ID") %>%
      select(ID, everything(),-Sample_ID) %>%
      pivot_longer(-ID, names_to = "Feature", values_to = "Value") %>%
      pivot_wider(names_from = ID, values_from = Value) %>%
      separate(Feature, into = c("mode", "mz", "rt"), sep = "_", remove = FALSE) %>%
      relocate(mode, mz, rt, .after = Feature) %>%
      mutate(ESI = case_when(
        mode == "HILIC" ~ "pos",
        mode == "C18" ~ "neg",
        TRUE ~ NA_character_)) %>%
      relocate(ESI, .after = mode)
  #- 1.2.2: Name, make features, transpose for HILIC
    HILIC_U <- process_feature_table(UFT_HILIC, HILICmap, "HILIC") %>%
      left_join(tumor_IDs, by = "Sample_ID") %>%
      select(ID, everything(), -Sample_ID) %>%
      pivot_longer(-ID, names_to = "Feature", values_to = "Value") %>%
      pivot_wider(names_from = ID, values_from = Value) %>%
      separate(Feature, into = c("mode", "mz", "rt"), sep = "_", remove = FALSE) %>%
      relocate(mode, mz, rt, .after = Feature) %>%
      mutate(ESI = case_when(
        mode == "HILIC" ~ "pos",
        mode == "C18" ~ "neg",
        TRUE ~ NA_character_)) %>%
      relocate(ESI, .after = mode)
  #- 1.2.3: Combine HILIC and C18 UFTs
    UFT <- rbind(HILIC_U, C18_U)
#+ 1.3 Targeted Processing
  #- 1.3.1: Clean sample column names in C18 table
    sample_cols <- names(TFT_C18)[!(names(TFT_C18) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
    cleaned_names <- gsub("\\.mzXML$", "", sample_cols)
    rename_vector <- setNames(C18map$Sample_ID, C18map$File_Name)
    new_names <- rename_vector[cleaned_names]
    names(TFT_C18)[match(sample_cols, names(TFT_C18))] <- new_names
    id_map <- setNames(tumor_IDs$ID, tumor_IDs$Sample_ID)
  #- 1.3.2: Collapse technical replicates via median
    TFT_C18_median <- TFT_C18 %>%
      # Pivot to long format to work with sample columns
      pivot_longer(cols = -(1:5), names_to = "Sample_ID_full", values_to = "Value") %>%
      # Extract base sample ID (before the underscore)
      mutate(Sample_ID = str_extract(Sample_ID_full, "^[^_]+")) %>%
      # Group and summarize by median
      group_by(across(1:5), Sample_ID) %>%
      summarise(Value = median(Value, na.rm = TRUE), .groups = "drop") %>%
      # Pivot back to wide
      pivot_wider(names_from = Sample_ID, values_from = Value)
  #- 1.3.4: Strip off anything after first underscore for q4June2022 samples
    sample_cols <- names(TFT_C18_median)[!(names(TFT_C18_median) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
    new_sample_cols <- ifelse(
      str_starts(sample_cols, "q4June2022"),
      str_extract(sample_cols, "^[^_]+"),
      sample_cols
    )
    names(TFT_C18_median)[match(sample_cols, names(TFT_C18_median))] <- new_sample_cols
  #- 1.3.5: Final rename using tumor_IDs
    id_map <- setNames(tumor_IDs$ID, tumor_IDs$Sample_ID)
    sample_cols_median <- names(TFT_C18_median)[!(names(TFT_C18_median) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
    new_names_median <- id_map[sample_cols_median]
    names(TFT_C18_median)[match(sample_cols_median, names(TFT_C18_median))] <- new_names_median
  #- 1.3.6: Prep for joining
    TFT_C18_median <- TFT_C18_median %>%
      select(mz, time, mz.1, Name, delta_ppm_vec, everything()) %>%
      mutate(ESI = "neg") %>%
      relocate(ESI, .after = delta_ppm_vec)
  #- 1.3.7: Clean sample column names in HILIC table
    sample_cols_HILIC <- names(TFT_HILIC)[!(names(TFT_HILIC) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
    cleaned_names_HILIC <- gsub("\\.mzXML$", "", sample_cols_HILIC)
    rename_vector_HILIC <- setNames(HILICmap$Sample_ID, HILICmap$File_Name)
    new_names_HILIC <- rename_vector_HILIC[cleaned_names_HILIC]
    names(TFT_HILIC)[match(sample_cols_HILIC, names(TFT_HILIC))] <- new_names_HILIC
  #- 1.3.8: Collapse technical replicates via median
    TFT_HILIC_median <- TFT_HILIC %>%
      pivot_longer(cols = -(1:5), names_to = "Sample_ID_full", values_to = "Value") %>%
      mutate(Sample_ID = str_extract(Sample_ID_full, "^[^_]+")) %>%
      group_by(across(1:5), Sample_ID) %>%
      summarise(Value = median(Value, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Sample_ID, values_from = Value)
  #- 1.3.9: Strip off anything after first underscore for q4June2022 samples
    sample_cols_HILIC_median <- names(TFT_HILIC_median)[!(names(TFT_HILIC_median) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
    new_sample_cols_HILIC <- ifelse(
      str_starts(sample_cols_HILIC_median, "q4June2022"),
      str_extract(sample_cols_HILIC_median, "^[^_]+"),
      sample_cols_HILIC_median
    )
    names(TFT_HILIC_median)[match(sample_cols_HILIC_median, names(TFT_HILIC_median))] <- new_sample_cols_HILIC
  #- 1.3.10: Final rename using tumor_IDs
    id_map_HILIC <- setNames(tumor_IDs$ID, tumor_IDs$Sample_ID)
    sample_cols_HILIC_final <- names(TFT_HILIC_median)[!(names(TFT_HILIC_median) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
    new_names_HILIC_final <- id_map_HILIC[sample_cols_HILIC_final]
    names(TFT_HILIC_median)[match(sample_cols_HILIC_final, names(TFT_HILIC_median))] <- new_names_HILIC_final
  #- 1.3.11: Prep for joining
    TFT_HILIC_median <- TFT_HILIC_median %>%
      select(mz, time, mz.1, Name, delta_ppm_vec, everything()) %>%
      mutate(ESI = "pos") %>%
      relocate(ESI, .after = delta_ppm_vec)
#+ 1.4: Combine targeted datasets
  TFT <- rbind(TFT_C18_median, TFT_HILIC_median) %>%
    mutate(mode = case_when(
      ESI == "pos" ~ "HILIC",
      ESI == "neg" ~ "C18",
      TRUE ~ NA_character_
    )) %>%
    relocate(mode, .before = ESI)
#+ 1.5: Export processed datasets
  write.csv(UFT, "Raw_Data/Final/untargeted_FT.csv", row.names = FALSE)
  write.csv(TFT, "Raw_Data/Final/targeted_FT.csv", row.names = FALSE)
#+ 1.6: Prepare data for Metaboanalyst
  #- 1.6.1: Set up targeted data for metaboanalyst
    #_ 1.6.1.1: Transform targeted dataset
      TFT_metaboanalyst_i <- TFT %>%
        mutate(Feature = paste0(Name, "_", mode, "_", mz)) %>%
        select(Feature, everything(),-c(mode, ESI, mz, time, mz.1, Name, delta_ppm_vec, nist_1:q4_6)) %>%
        column_to_rownames("Feature") %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("Patient_ID") %>%
        mutate(Variant = case_when(
          str_detect(Patient_ID, "^FVPTC\\d+$") ~ "FV-PTC",
          str_detect(Patient_ID, "^F\\d+$") ~ "FTC",
          str_detect(Patient_ID, "^P\\d+$") ~ "PTC",
          TRUE ~ "Unknown"
        )) %>%
        relocate(Variant, .after = Patient_ID) %>%
        as_tibble()
    #_ 1.6.1.2: Remove any columns with > 20% missing values
      feature_cols_tft <- names(TFT_metaboanalyst_i)[!names(TFT_metaboanalyst_i) %in% c("Patient_ID", "Variant")]
      zero_percentages_tft <- TFT_metaboanalyst_i %>%
        select(all_of(feature_cols_tft)) %>%
        summarise(across(everything(), ~ sum(.x == 0, na.rm = TRUE) / length(.x))) %>%
        pivot_longer(everything(), names_to = "feature", values_to = "zero_pct")
      #_Identify features with <= 20% zeros (keep these)
      features_to_keep_tft <- zero_percentages_tft %>%
        filter(zero_pct <= 0.20) %>%
        pull(feature)
      #_Filter TFT_metaboanalyst to keep only good features
      TFT_metaboanalyst_raw <- TFT_metaboanalyst_i %>%
        select(Patient_ID, Variant, all_of(features_to_keep_tft))
    #_ 1.6.1.3: Replace 0 values which remain with 1/2 the column minimum
      TFT_metaboanalyst_halfmin <- TFT_metaboanalyst_raw %>%
        mutate(across(all_of(features_to_keep_tft), ~ ifelse(.x == 0, 0.5 * min(.x[.x > 0], na.rm = TRUE), .x)))
    #_ 1.6.1.4: Log2 transform dataset
      TFT_metaboanalyst_log2 <- TFT_metaboanalyst_halfmin %>%
        mutate(across(all_of(features_to_keep_tft), ~ log2(.x))) %>%
        mutate(Variant = as.factor(Variant))
  #- 1.6.2: Set up untargeted data for metaboanalyst
    #_ 1.6.2.1: Transform untargeted dataset
      UFT_metaboanalyst_i <- UFT %>%
        select(-c(mode, ESI, mz, rt, nist_1:q4_6)) %>%
        column_to_rownames("Feature") %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("Patient_ID") %>%
        mutate(Variant = case_when(
          str_detect(Patient_ID, "^FVPTC\\d+$") ~ "FV-PTC",
          str_detect(Patient_ID, "^F\\d+$") ~ "FTC",
          str_detect(Patient_ID, "^P\\d+$") ~ "PTC",
          TRUE ~ "Unknown"
        )) %>%
        relocate(Variant, .after = Patient_ID) %>%
        as_tibble()
    #_ 1.6.2.2: Remove any columns with > 20% missing values
      feature_cols <- names(UFT_metaboanalyst_i)[!names(UFT_metaboanalyst_i) %in% c("Patient_ID", "Variant")]
      zero_percentages <- UFT_metaboanalyst_i %>%
        select(all_of(feature_cols)) %>%
        summarise(across(everything(), ~ sum(.x == 0, na.rm = TRUE) / length(.x))) %>%
        pivot_longer(everything(), names_to = "feature", values_to = "zero_pct")
      #_Identify features with <= 20% zeros (keep these)
      features_to_keep <- zero_percentages %>%
        filter(zero_pct <= 0.20) %>%
        pull(feature)
      #_Filter UFT_metaboanalyst to keep only good features
      UFT_metaboanalyst_raw <- UFT_metaboanalyst_i %>%
        select(Patient_ID, Variant, all_of(features_to_keep))
    #_ 1.6.2.3: Replace 0 values which remain with 1/2 the column minimum
      UFT_metaboanalyst_halfmin <- UFT_metaboanalyst_raw %>%
        mutate(across(all_of(features_to_keep), ~ ifelse(.x == 0, 0.5 * min(.x[.x > 0], na.rm = TRUE), .x)))
    #_ 1.6.2.4: Log2 transform dataset
      UFT_metaboanalyst_log2 <- UFT_metaboanalyst_halfmin %>%
        mutate(across(all_of(features_to_keep), ~ log2(.x))) %>%
        mutate(Variant = as.factor(Variant))
    #_ 1.6.2.5: Split UFT into C18 and HILIC datasets
      hilic_cols <- names(UFT_metaboanalyst_log2)[str_starts(names(UFT_metaboanalyst_log2), "HILIC")]
      c18_cols <- names(UFT_metaboanalyst_log2)[str_starts(names(UFT_metaboanalyst_log2), "C18")]
      UFT_HILIC_metaboanalyst_log2 <- UFT_metaboanalyst_log2 %>%
        select(Patient_ID, Variant, all_of(hilic_cols))
      UFT_C18_metaboanalyst_log2 <- UFT_metaboanalyst_log2 %>%
        select(Patient_ID, Variant, all_of(c18_cols))
  #- 1.6.3: Export for metaboanalyst
    write.csv(TFT_metaboanalyst_log2, "Raw_Data/Metaboanalyst/TFT_metaboanalyst_log2.csv", row.names = FALSE)
    write.csv(UFT_metaboanalyst_raw, "Raw_Data/Metaboanalyst/UFT_metaboanalyst_raw.csv", row.names = FALSE)
    write.csv(UFT_C18_metaboanalyst_log2, "Raw_Data/Metaboanalyst/UFT_C18_metaboanalyst_log2.csv", row.names = FALSE)
    write.csv(UFT_HILIC_metaboanalyst_log2, "Raw_Data/Metaboanalyst/UFT_HILIC_metaboanalyst_log2.csv", row.names = FALSE)
    #! Performed HCA with heatmaps in metaboanalyst
#+ 1.7: Cleanup tumor pathology data
  #- 1.7.1: Read and preprocess tumor pathology data
    tumor_pathology <- read_excel("Raw_Data/Final/tumor_pathology.xlsx", sheet = "pathology") %>%
    mutate(LVI = ifelse(LVI == 2, NA, LVI)) %>%
    assign_T_stage(ld_col = "LD", ete_col = "ETE", units = "cm", out_col = "T_stage") %>%
    mutate(
      T_computed = factor(T_stage, levels = c("T1", "T2", "T3", "T4"), ordered = TRUE),
      T_computed_bin = factor(
        dplyr::case_when(
          T_computed %in% c("T1", "T2") ~ "T1-T2",
          T_computed %in% c("T3", "T4") ~ "T3-T4",
          TRUE ~ NA_character_
        ),
        levels = c("T1-T2", "T3-T4"), ordered = TRUE
      ),
      LVI = factor(
        dplyr::case_when(
          is.na(LVI) ~ NA_character_,
          LVI == 1 ~ "+LVI",
          LVI == 0 ~ "-LVI"
        ),
        levels = c("-LVI", "+LVI")
      ),
      Sex = factor(
        dplyr::case_when(
          is.na(Sex) ~ NA_character_,
          Sex == 1 ~ "Female",
          Sex == 0 ~ "Male"
        ),
        levels = c("Female", "Male")
      ),
      MFC = factor(
        dplyr::case_when(
          is.na(MFC) ~ NA_character_,
          MFC == 1 ~ "+MFC",
          MFC == 0 ~ "-MFC"
        ),
        levels = c("-MFC", "+MFC")
      ),
      Age = as.numeric(Age)
    ) %>%
    select(Patient_ID, Sex, Age, MFC, LD, LVI, ETE, T_computed, T_computed_bin, LD, LVI, T_computed)
  #- 1.7.2: Combine with UFT and TFT datasets
    UFT_metaboanalyst_log2_path <- UFT_metaboanalyst_log2 %>%
      left_join(tumor_pathology, by = "Patient_ID") %>%
      select(Patient_ID, Variant, colnames(tumor_pathology), everything()) %>%
      mutate(Variant = factor(Variant, levels = c("PTC", "FV-PTC", "FTC")))
    TFT_metaboanalyst_log2_path <- TFT_metaboanalyst_log2 %>%
      left_join(tumor_pathology, by = "Patient_ID") %>%
      select(Patient_ID, Variant, colnames(tumor_pathology), everything()) %>%
      mutate(Variant = factor(Variant, levels = c("PTC", "FV-PTC", "FTC")))
  #- 1.7.3: Define metadata and feature columns
    metadata_cols <- c("Patient_ID","Variant", "Sex","Age", "MFC", "LD", "LVI", "ETE", "T_computed", "T_computed_bin")
    feature_cols <- grep("^(C18|HILIC)", colnames(UFT_metaboanalyst_log2_path), value = TRUE)

