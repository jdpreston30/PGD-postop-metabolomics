#* 1: Importing feature tables and preprocessing
#+ 1.0: Pull PGD info
pgd_status <- clinical_metadata_i |>
  select(Patient, PGD = postop_PGD_ISHLT)
#+ 1.1:Import FTs, add Patient and Sample IDs, filter to S0 (preop); apply uniqueness filter
#- 1.1.1: TFT annotated and preprocess; apply uniqueness filter
TFT_annot <- read_csv(config$paths$TFT_annot) |>
  preprocess_FT(apply_unique_filter = TRUE, unique_threshold = 0.8)
#- 1.1.1: UFT filtered and preprocess; NO uniqueness filter
UFT_filtered_i <- read_csv(config$paths$UFT_filtered) |>
  preprocess_FT(apply_unique_filter = FALSE)
#- 1.1.1: UFT filtered and preprocess; NO uniqueness filter
UFT_i <- read_csv(config$paths$UFT_full) |>
  preprocess_FT(apply_unique_filter = FALSE)
# ! Filtering out H49 as this was a 'false start' where we collected early but then heart offer didn't go through
# ! Fixed typo in Sample_ID for H46SS0
#+ 1.2: Import MSMICA feature key; join with relevant QC info
#- 1.2.1: Read in pared/procured key
TFT_QC <- read_xlsx(config$paths$manual_QC, sheet = "pared") |>
  arrange(p_value) |>
  mutate(sig_ord = row_number()) |>
  select(sig_ord, feature, display_name, long_name, note, adduct)
#- 1.2.2: Read in key; join; clean column names
TFT_annot_key <- read_csv(config$paths$TFT_annot_key) |>
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
  )
#+ 1.3: Bring in IROA IDX feature library
idx_lib <- read_csv(config$paths$idx_library)
#+ 1.4: Create Identified TFT_confirmed based on library
#- 1.4.1: Run function to match features to library and create identified TFT_confirmed
identified <- create_identified_FT(
  feature_table = UFT_filtered_i,
  reference_library = idx_lib,
  mz_thresh_ppm = 5,
  time_thresh_sec = 30
)
#- 1.4.2: Assign and apply unique filter fxn just in case
TFT_confirmed <- identified$TFT_confirmed
#- 1.4.3: Build key and filter out entries without compound names
TFT_confirmed_key <- identified$matched_features |>
  filter(!is.na(compound_name)) |>  # Remove unmatched features
  select(identified_name = compound_name, isomer = library_isomer, everything()) |>
  mutate(MMD = "")
#+ 1.5: Create merged library/annotated TFT
#- 1.5.1: Subset features from annotations with source tracking
TFT_annot_features <- TFT_annot_key |>
  mutate(lib_conf = "N", source = "annotation") |>
  select(feature, lib_conf, source) |>
  unique()
#- 1.5.2: Subset features from confirmed with source tracking
TFT_confirmed_features <- TFT_confirmed_key |>
  mutate(lib_conf = "Y", source = "library") |>
  select(feature, lib_conf, source) |>
  unique()
#- 1.5.3: Identify features present in both datasets
overlapping_features_lookup <- setdiff(intersect(TFT_annot_features$feature, TFT_confirmed_features$feature), "Sample")
#- 1.5.4: Create source labels for overlapping and unique features
TFT_merged_features <- bind_rows(TFT_confirmed_features, TFT_annot_features) |>
  arrange(feature, desc(lib_conf)) |>
  group_by(feature) |>
  summarise(
    lib_conf = first(lib_conf), # Prioritize "Y" over "N" due to arrange
    source = if_else(
      n() > 1,
      "both",
      dplyr::first(source)
    ),
    .groups = "drop"
  )
#- 1.5.5: Identify overlapping feature columns to handle duplicates properly
annot_features <- names(TFT_annot)[!names(TFT_annot) %in% c("Patient", "Sample", "severe_PGD", "PGD_grade_tier", "any_PGD")]
confirmed_features <- names(TFT_confirmed)[!names(TFT_confirmed) %in% c("Patient", "Sample", "severe_PGD", "PGD_grade_tier", "any_PGD")]
overlapping_features <- intersect(annot_features, confirmed_features)
#- 1.5.6: Create base from annotations, removing overlapping features
TFT_combined_base <- TFT_annot |>
  select(-any_of(overlapping_features))
#- 1.5.7: Prepare confirmed features only (Patient + feature columns)
TFT_confirmed_features_only <- TFT_confirmed |>
  select(Patient, Sample, all_of(confirmed_features))
#- 1.5.8: Merge datasets and apply uniqueness filter
TFT_combined_i <- TFT_combined_base |>
  left_join(TFT_confirmed_features_only, by = c("Patient", "Sample")) |>
  filter_unique_features(unique_threshold = 0.8)
#+ 1.6: Create synthetic rows for missing timepoints
#- 1.6.1: Define missing samples metadata
missing_samples <- tibble::tribble(
  ~Patient, ~Sample, ~severe_PGD,       ~PGD_grade_tier,   ~any_PGD,
  "H5",     "S2",    "Severe PGD",      "Severe PGD",      "Y",
  "H14",    "S2",    "No Severe PGD",   "Mild/Mod. PGD",   "Y",
  "H50",    "S1",    "No Severe PGD",   "Mild/Mod. PGD",   "Y"
)
#- 1.6.2: Apply synthetic row creation to all three feature tables; add Time and Sample_ID columns; recode severe_PGD
TFT_combined <- create_synthetic_rows(TFT_combined_i, missing_samples) |>
  mutate(
    Sample_ID = paste0(Patient, Sample),
    Time = if_else(Sample == "S1", 12, 24),
    severe_PGD = if_else(severe_PGD == "Severe PGD", "Y", "N")
  ) |>
  select(Patient, Sample, Sample_ID, Time, everything())
UFT_filtered <- create_synthetic_rows(UFT_filtered_i, missing_samples) |>
  mutate(
    Sample_ID = paste0(Patient, Sample),
    Time = if_else(Sample == "S1", 12, 24),
    severe_PGD = if_else(severe_PGD == "Severe PGD", "Y", "N")
  ) |>
  select(Patient, Sample, Sample_ID, Time, everything())
UFT <- create_synthetic_rows(UFT_i, missing_samples) |>
  mutate(
    Sample_ID = paste0(Patient, Sample),
    Time = if_else(Sample == "S1", 12, 24),
    severe_PGD = if_else(severe_PGD == "Severe PGD", "Y", "N")
  ) |>
  select(Patient, Sample, Sample_ID, Time, everything())
#- 1.6.3: Take the feature names and make a vector
untargeted_features <- names(UFT) |>
  purrr::keep(~ str_starts(.x, "HILIC") | str_starts(.x, "C18"))
