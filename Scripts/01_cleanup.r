#* 1: Cleanup
#+ Import Data
  #- Path set
    raw_path <- "/Users/jdp2019/Library/CloudStorage/OneDrive-Emory/Research/Manuscripts and Projects/Active Projects/Chan Lab/TPMO/Metabolomics - Postop PGD_Amshu/raw_data/"
  #- Import Metadata and Sequence and Key
    metadata <- read_csv(paste0(raw_path, "metadata.csv"))
    sequence <- read_csv(paste0(raw_path, "sequence.csv")) %>%
      mutate(
        Sample_ID = Sample_ID %>%
          str_replace_all("NIST_1", "NIST1") %>%
          str_replace_all("NIST_2", "NIST2")
      )
    feature_key <- read_csv(paste0(raw_path, "targeted_feature_key.csv"))
  #- Import Untargeted
    Sys.setenv("VROOM_CONNECTION_SIZE" = 2^22)
    C18_unt_UP <- read_csv(paste0(raw_path, "/C18_neg_untargeted_FT.csv"))
    HILIC_unt_UP <- read_csv(paste0(raw_path, "/HILIC_pos_untargeted_FT.csv"))
  #- Import targeted
    C18_targeted_raw <- read_csv(paste0(raw_path, "Tidy/C18_targeted.csv"))
    HILIC_targeted_raw <- read_csv(paste0(raw_path, "Tidy/HILIC_targeted.csv"))
#+ Prep to append asterisks onto targeted with isomers for targeted
  #- Get list of targeted annotations with more than one identity
    collapsed_all_isomers <- feature_key %>%
      mutate(Mode = sub(".*_(HILIC|C18)$", "\\1", Name)) %>%
      select(Name, KEGGID, Mode, mz_time, Formula, Identified_Name) %>%
      mutate(.rowid = row_number()) %>% # preserve original order
      arrange(.rowid) %>%
      group_by(Mode, mz_time) %>%
      summarise(
        Metabolite = first(Name),
        primary_KEGGID = first(KEGGID),
        primary_Name = first(Identified_Name),
        other_KEGGIDs = if (n() > 1) paste(KEGGID[-1], collapse = ", ") else NA_character_,
        other_Names = if (n() > 1) paste(Identified_Name[-1], collapse = ", ") else NA_character_,
        n_candidates = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(n_candidates))
  #- Filter to duplicate metabolites
    dup_mets <- collapsed_all_isomers %>%
      filter(n_candidates > 1) %>%
      distinct(Metabolite) %>%
      pull()
  #- Compact lookup table of the duplicates
    dup_lookup <- collapsed_all_isomers %>%
      filter(n_candidates > 1) %>%
      select(Metabolite, primary_KEGGID, primary_Name, other_KEGGIDs, other_Names, n_candidates)
  #- Add asterisks to original feature key
    dup_mets <- collapsed_all_isomers %>%
      filter(n_candidates > 1) %>%
      pull(Metabolite)
    feature_key_isomark <- feature_key %>%
      mutate(Rename = if_else(Name %in% dup_mets, paste0(Name, "*"), Name))
  #- Create a clean key for rejoining
    combined_key <- feature_key_isomark %>%
      rename(Metabolite = Name) %>%
      select(Rename, Metabolite, Identified_Name)
#+ Preprocess Untargeted Feature Tables into Tidy format
  #! Targeted FTs were already processed into this format during annotation with MSMICA
  #- Set QC pattern
    qc_pattern <- "^(NIST1|NIST2|NIST_1|NIST_2|q[1-9]+)$" # adjust if needed
  #- Create synthetic rows for the two missing time points
    missing_pts <- metadata %>%
      filter(Patient %in% c("H2", "H5")) %>%
      distinct(Patient, Sex, Age, Clinical_PGD) %>%
      mutate(
        Time = 12,
        Sample_ID = paste0(Patient, "_", Time)
      ) %>%
      select(Sample_ID, Patient, Time, Sex, Age, Clinical_PGD)
  #- For C18
    C18_unt_proc <- process_feature_table(C18_unt_UP, sequence, "C18") %>%
      filter(Sample_ID != "Blank") %>%
      filter(!str_detect(Sample_ID, "(D0|S0)$")) %>%
      add_row(Sample_ID = "H19S2") %>% #! H19S2 missing for C18 only
      arrange(Sample_ID) 
  #- For HILIC
    HILIC_unt_proc <- process_feature_table(HILIC_unt_UP, sequence, "HILIC") %>%
      filter(Sample_ID != "Blank") %>%
      filter(!str_detect(Sample_ID, "(D0|S0)$")) %>%
      arrange(Sample_ID) 
  #- Join together, add metadata
    C18_HILIC_unt_proc <- full_join(C18_unt_proc, HILIC_unt_proc, by = "Sample_ID") %>%
      select(Sample_ID, starts_with("HILIC"), starts_with("C18")) %>%
      left_join(metadata, by = "Sample_ID") %>%
      select(Sample_ID, all_of(colnames(metadata)), everything()) %>%
        # rebuild Sample_ID only for non-QC rows; QC rows keep original
        mutate(
          .orig_id = Sample_ID,
          Sample_ID = if_else(
            str_detect(.orig_id, qc_pattern) | is.na(Patient) | is.na(Time),
            .orig_id,
            paste0(as.character(Patient), "_", as.character(Time))
          )
        ) %>%
      select(Sample_ID, everything(), -.orig_id) %>%
      bind_rows(missing_pts) %>%
      arrange(Patient, Time) %>%
      mutate(
        Sample_ID = factor(Sample_ID),
        Patient = factor(Patient),
        Time = factor(Time, levels = c(12, 24)),
        Clinical_PGD = factor(Clinical_PGD, levels = c("N", "Y")),
        Sex = factor(Sex, levels = c("M", "F"))
      )
#+ Clean up targeted
  #- For C18
    C18_targeted <- C18_targeted_raw %>%
      filter(Sample_ID != "Blank_1" & Sample_ID != "Blank_2") %>%
      mutate(
        Sample_ID = case_when(
          Sample_ID == "NIST_1" ~ "NIST1",
          Sample_ID == "NIST_2" ~ "NIST2",
          TRUE ~ Sample_ID
        )
      ) %>%
      filter(!str_detect(Sample_ID, "(D0|S0)$")) %>%
      add_row(Sample_ID = "H19S2") %>% #! H19S2 missing for C18 only
      arrange(Sample_ID)
  #- For HILIC
    HILIC_targeted <- HILIC_targeted_raw %>%
      filter(Sample_ID != "Blank_1" & Sample_ID != "Blank_2") %>%
      mutate(
        Sample_ID = case_when(
          Sample_ID == "NIST_1" ~ "NIST1",
          Sample_ID == "NIST_2" ~ "NIST2",
          TRUE ~ Sample_ID
        )
      ) %>%
      filter(!str_detect(Sample_ID, "(D0|S0)$")) %>%
      arrange(Sample_ID)
  #- Set metadata and feature columns
    meta_cols <- c("Patient", "Time", "Clinical_PGD", "Sex", "Age", "Sample_ID")
    feature_cols_targ <- setdiff(names(C18_HILIC_targeted), meta_cols)
  #- Join together, add metadata, rename columns with asterisks append
    C18_HILIC_targeted_i <- full_join(C18_targeted, HILIC_targeted, by = "Sample_ID") %>%
      left_join(metadata, by = "Sample_ID") %>%
      select(Sample_ID, all_of(colnames(metadata)), everything()) %>%
      mutate(
        .orig_id = Sample_ID,
        Sample_ID = if_else(
          str_detect(.orig_id, qc_pattern) | is.na(Patient) | is.na(Time),
          .orig_id,
          paste0(as.character(Patient), "_", as.character(Time))
        )
      ) %>%
      select(Sample_ID, everything(), -.orig_id) %>%
      bind_rows(missing_pts) %>%
      arrange(Patient, Time) %>%
      mutate(
        Sample_ID = factor(Sample_ID),
        Patient = factor(Patient),
        Time = factor(Time, levels = c(12, 24)),
        Clinical_PGD = factor(Clinical_PGD, levels = c("N", "Y")),
        Sex = factor(Sex, levels = c("M", "F"))
      )
  #- Prep name map
    name_map <- combined_key %>%
      distinct(Metabolite, Rename) %>%
      filter(Metabolite %in% names(C18_HILIC_targeted_i)) %>%
      tibble::deframe()
  #- Rename with asterisks
    C18_HILIC_targeted <- C18_HILIC_targeted_i %>%
      rename_with(~ unname(name_map[.x]), .cols = all_of(names(name_map)))
#+ Transformations and QC
  #- Run cleanup function on Targeted
    targeted_QCd <- FT_QC(
      df          = C18_HILIC_targeted,
      id_col      = "Sample_ID",
      common_cols = meta_cols,
      qc_remove   = TRUE # QC used for half-min & 20% check, then removed
    )
  #- Run cleanup function on Untargeted
    untargeted_QCd <- FT_QC(
      df          = C18_HILIC_unt_proc,
      id_col      = "Sample_ID",
      common_cols = meta_cols,
      qc_remove   = TRUE # QC used for half-min & 20% check, then removed
    )
#+ Final Assignments of feature tables and column names
  #- Untargeted
    UFT_metaboanalyst <- untargeted_QCd$full
    UFT <- untargeted_QCd$filtered
    UFT_nomiss <- UFT %>%
      filter(Sample_ID != "H2_12" & Sample_ID != "H5_12" & Sample_ID != "H19_24")
  #- Targeted
    TFT <- targeted_QCd$filtered
  #- Columns
    untargeted_features <- names(UFT) %>%
      purrr::keep(~ str_starts(.x, "HILIC") | str_starts(.x, "C18"))
