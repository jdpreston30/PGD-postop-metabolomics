#! ALL FROM OTHER REPO NEEDS TO BE UPDATED
#* 0d: Importing feature tables and preprocessing
#+ 0d.0: Pull PGD info
pgd_status <- clinical_metadata_i %>%
  select(Patient, PGD = postop_PGD_ISHLT)
#+ 0d.1:Import FTs, add Patient and Sample IDs, filter to S0 (preop); apply uniqueness filter
#- 0d.1.1: TFT annotated and preprocess; apply uniqueness filter
TFT_annot <- read_csv(config$paths$TFT_annot) %>%
  preprocess_FT(apply_unique_filter = TRUE, unique_threshold = 0.8)
#- 0d.1.1: UFT filtered and preprocess; NO uniqueness filter
UFT_filtered_i <- read_csv(config$paths$UFT_filtered) %>%
  preprocess_FT(apply_unique_filter = FALSE)
#- 0d.1.1: UFT filtered and preprocess; NO uniqueness filter
UFT_i <- read_csv(config$paths$UFT_full) %>%
  preprocess_FT(apply_unique_filter = FALSE)
# ! Filtering out H49 as this was a 'false start' where we collected early but then heart offer didn't go through
# ! Fixed typo in Sample_ID for H46SS0
#+ 0d.2: Import MSMICA feature key; join with relevant QC info
#- 0d.2.1: Read in QC Info
TFT_annot_QC <- read_xlsx(config$paths$TFT_annot_QC, sheet = "key") %>%
  select(short_name, include, human, microbial, plant, dietary, pharmacologic, environmental, chem_class, rationale, Name = name, isomer_names, abbreviated)
#- 0d.2.2: Read in pared/procured key
TFT_QC <- read_xlsx(config$paths$manual_QC, sheet = "pared") %>%
  arrange(p_value) %>%
  mutate(sig_ord = row_number()) %>%
  select(sig_ord, feature, display_name, long_name, note, adduct)
#- 0d.2.1: Read in key; join; clean column names
TFT_annot_key <- read_csv(config$paths$TFT_annot_key) %>%
  rename(
    feature = Feature,
    identified_name = `Identified Name`,
    rt = `Retention Time`,
    exact_mass = `Exact Mass`,
    ion_mode = `Ion Mode`,
    MMD = `Multi-Mode Detection`,
    mean_intensity = `Mean Intensity`,
    mean_intensity_w_zeros = `Mean Intensity w/ Zeros`,
    identification_method = `Identification Method`,
    annotation_probability = `Annotation Probability`,
    KEGG = `KEGG ID`,
    CID = `PubChem CID`,
    inchi_key = InChIKey,
    HMDB = `HMDB ID`,
    alternative_parent = `Alternative Parent`,
    isomer = Isomer
  ) %>%
  left_join(TFT_annot_QC, by = "Name")
#+ 0d.3: Bring in IROA IDX feature library
qstd_c18 <- read_tsv(config$paths$qstd_c18) %>%
  mutate(lib = "C18")
qstd_hilic <- read_tsv(config$paths$qstd_hilic) %>%
  mutate(lib = "HILIC")
qstd_quant <- bind_rows(qstd_c18, qstd_hilic)
idx_lib <- read_csv(config$paths$idx_library)
#+ 0d.4: Create Identified TFT_confirmed based on library
#- 0d.4.1: Run function to match features to library and create identified TFT_confirmed
identified <- create_identified_FT(
  feature_table = UFT_filtered_i,
  reference_library = idx_lib,
  mz_thresh_ppm = 5,
  time_thresh_sec = 30
)
#- 0d.4.2: Assign and apply unique filter fxn just in case
TFT_confirmed <- identified$TFT_confirmed
#- 0d.4.3: Build key and filter out entries without compound names
TFT_confirmed_key <- identified$matched_features %>%
  filter(!is.na(compound_name)) %>%  # Remove unmatched features
  select(identified_name = compound_name, isomer = library_isomer, everything()) %>%
  mutate(MMD = "")
#+ 0d.5: Create merged library/annotated TFT
#- 0d.5.1: Subset features from annotations with source tracking
TFT_annot_features <- TFT_annot_key %>%
  mutate(lib_conf = "N", source = "annotation") %>%
  select(feature, lib_conf, source) %>%
  unique()
#- 0d.5.2: Subset features from confirmed with source tracking
TFT_confirmed_features <- TFT_confirmed_key %>%
  mutate(lib_conf = "Y", source = "library") %>%
  select(feature, lib_conf, source) %>%
  unique()
#- 0d.5.3: Identify features present in both datasets
overlapping_features_lookup <- setdiff(intersect(TFT_annot_features$feature, TFT_confirmed_features$feature), "Sample")
#- 0d.5.4: Create source labels for overlapping and unique features
TFT_merged_features <- bind_rows(TFT_confirmed_features, TFT_annot_features) %>%
  arrange(feature, desc(lib_conf)) %>%
  group_by(feature) %>%
  summarise(
    lib_conf = first(lib_conf), # Prioritize "Y" over "N" due to arrange
    source = if_else(
      n() > 1,
      "both",
      dplyr::first(source)
    ),
    .groups = "drop"
  )
#- 0d.5.5: Identify overlapping feature columns to handle duplicates properly
annot_features <- names(TFT_annot)[!names(TFT_annot) %in% c("Patient", "Sample", "severe_PGD", "PGD_grade_tier", "any_PGD")]
confirmed_features <- names(TFT_confirmed)[!names(TFT_confirmed) %in% c("Patient", "Sample", "severe_PGD", "PGD_grade_tier", "any_PGD")]
overlapping_features <- intersect(annot_features, confirmed_features)
#- 0d.5.6: Create base from annotations, removing overlapping features
TFT_combined_base <- TFT_annot %>%
  select(-any_of(overlapping_features))
#- 0d.5.7: Prepare confirmed features only (Patient + feature columns)
TFT_confirmed_features_only <- TFT_confirmed %>%
  select(Patient, Sample, all_of(confirmed_features))
#- 0d.5.8: Merge datasets and apply uniqueness filter
TFT_combined_i <- TFT_combined_base %>%
  left_join(TFT_confirmed_features_only, by = c("Patient", "Sample")) %>%
  filter_unique_features(unique_threshold = 0.8)
#- 0d.5.9: Create combined key for TFT_combined
TFT_combined_key <- bind_rows(
  TFT_confirmed_key %>% 
    filter(!is.na(identified_name)) %>%
    mutate(source = "library", lib_conf = "Y"),
  TFT_annot_key %>% 
    filter(!is.na(identified_name)) %>%
    mutate(source = "annotation", lib_conf = "N")
) %>%
  left_join(TFT_merged_features %>% select(feature, source_final = source), by = "feature") %>%
  mutate(source = coalesce(source_final, source)) %>%
  select(-source_final) %>%
  arrange(feature, desc(lib_conf)) %>%
  group_by(feature) %>%
  dplyr::slice(1) %>%  # Keep library version when overlap
  ungroup()
#- 0d.5.10: Create clean key for downstream use (volcano/balloon plots)
combined_key_clean <- TFT_combined_key %>%
  group_by(identified_name) %>%
  mutate(
    is_isomer = n() > 1 | isomer %in% c("Y", "yes"),
    is_library = lib_conf == "Y",
    Identified_Name = identified_name,
    Identified_Name = if_else(is_isomer, paste0(Identified_Name, "*"), Identified_Name),
    Identified_Name = if_else(is_library, paste0(Identified_Name, "\u1D62"), Identified_Name)
  ) %>%
  ungroup() %>%
  select(Metabolite = feature, Identified_Name)
#+ 0d.6: Create synthetic rows for missing timepoints
#- 0d.6.1: Define missing samples metadata
missing_samples <- tibble::tribble(
  ~Patient, ~Sample, ~severe_PGD,       ~PGD_grade_tier,   ~any_PGD,
  "H5",     "S2",    "Severe PGD",      "Severe PGD",      "Y",
  "H14",    "S2",    "No Severe PGD",   "Mild/Mod. PGD",   "Y",
  "H50",    "S1",    "No Severe PGD",   "Mild/Mod. PGD",   "Y"
)
#- 0d.6.2: Apply synthetic row creation to all three feature tables; add Time and Sample_ID columns; recode severe_PGD
TFT_combined <- create_synthetic_rows(TFT_combined_i, missing_samples) %>%
  mutate(
    Sample_ID = paste0(Patient, Sample),
    Time = if_else(Sample == "S1", 12, 24),
    severe_PGD = if_else(severe_PGD == "Severe PGD", "Y", "N")
  ) %>%
  select(Patient, Sample, Sample_ID, Time, everything())
UFT_filtered <- create_synthetic_rows(UFT_filtered_i, missing_samples) %>%
  mutate(
    Sample_ID = paste0(Patient, Sample),
    Time = if_else(Sample == "S1", 12, 24),
    severe_PGD = if_else(severe_PGD == "Severe PGD", "Y", "N")
  ) %>%
  select(Patient, Sample, Sample_ID, Time, everything())
UFT <- create_synthetic_rows(UFT_i, missing_samples) %>%
  mutate(
    Sample_ID = paste0(Patient, Sample),
    Time = if_else(Sample == "S1", 12, 24),
    severe_PGD = if_else(severe_PGD == "Severe PGD", "Y", "N")
  ) %>%
  select(Patient, Sample, Sample_ID, Time, everything())
#- 0d.6.3: Take the feature names and make a vector
untargeted_features <- names(UFT) %>%
  purrr::keep(~ str_starts(.x, "HILIC") | str_starts(.x, "C18"))
