#* 1: Cleanup
#1.0 Set relative paths

  #- 1.6.2: Set up untargeted data for metaboanalyst

    #_ 1.6.2.2: Remove any columns with > 20% missing values
      feature_cols <- names(combined_unt)[!names(combined_unt) %in% c(common_cols, "Sample_ID")]
      zero_percentages <- combined_unt %>%
        select(all_of(feature_cols)) %>%
        summarise(across(everything(), ~ sum(.x == 0, na.rm = TRUE) / length(.x))) %>%
        pivot_longer(everything(), names_to = "feature", values_to = "zero_pct")
      #_Identify features with <= 20% zeros (keep these)
      features_to_keep <- zero_percentages %>%
        filter(zero_pct <= 0.20) %>%
        pull(feature)
      #_Filter UFT_metaboanalyst to keep only good features
      UFT_metaboanalyst_raw <- combined_unt %>%
        select(all_of(c(common_cols, "Sample_ID")), all_of(features_to_keep))
    #_ 1.6.2.3: Replace 0 values which remain with 1/2 the column minimum
      UFT_metaboanalyst_halfmin <- UFT_metaboanalyst_raw %>%
        mutate(across(all_of(features_to_keep), ~ ifelse(.x == 0, 0.5 * min(.x[.x > 0], na.rm = TRUE), .x)))
    #_ 1.6.2.4: Log2 transform dataset
      UFT_metaboanalyst_log2 <- UFT_metaboanalyst_halfmin %>%
        mutate(across(all_of(features_to_keep), ~ log2(.x)))
    #_ 1.6.2.5: Remove rows with false data or incomplete time points
      UFT_metaboanalyst_log2 <- UFT_metaboanalyst_log2[!((UFT_metaboanalyst_log2$Patient == 'H3')+(UFT_metaboanalyst_log2$Patient == 'H2')+(UFT_metaboanalyst_log2$Patient == 'H5')+(UFT_metaboanalyst_log2$Patient == 'H13'&UFT_metaboanalyst_log2$Time == -12)+(UFT_metaboanalyst_log2$Patient == 'H19'&UFT_metaboanalyst_log2$Time == 24)),]
      UFT_metaboanalyst_log2 <- UFT_metaboanalyst_log2[(UFT_metaboanalyst_log2$Time != -12),]
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

