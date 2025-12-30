#* 11 Supplementary Materials
#+ 11.1: Supplementary Material 1
#! This was a word document of the supplementary methods
#+ 11.2: Supplementary Material 2
#- 11.2.1: Structure confirmed key
SM2_L1_raw <- TFT_confirmed_key |>
  # Add the column type
  mutate(ion_mode = if_else(library_mode == "HILIC", "Positive", "Negative")) |>
  arrange(feature_mz, ion_mode) |>
  # Use function to clean up names
  mutate(identified_name = normalize_chem_names(identified_name)) |>
  # Mark if also in annotated key
  mutate(annotated = if_else(feature %in% TFT_annot_key$feature, "Y", "N")) |>
  # Create comma-delimited list of annotation names for this feature
  mutate(annotated_as = map_chr(feature, ~{
    matches <- TFT_annot_key |> filter(feature == .x) |> pull(identified_name)
    if(length(matches) > 0) paste(matches, collapse = ", ") else ""
  })) |>
  # Create comma-delimited list of OTHER identified names for this feature within Level 1
  mutate(other_L1_matches = map2_chr(feature, identified_name, ~{
    matches <- TFT_confirmed_key |> 
      filter(feature == .x, identified_name != .y) |> 
      pull(identified_name) |>
      normalize_chem_names()  # Apply normalization to the matches
    if(length(matches) > 0) paste(matches, collapse = ", ") else ""
  })) |>
  # Clean up isomer column
  mutate(isomer = case_when(
    is.na(isomer) ~ "N",
    isomer == "yes" ~ "Y", 
    str_detect(isomer, "yes -but different (time|RT)") ~ "Y, but different RT",
    TRUE ~ isomer
  ))
#- 11.2.2: Create clean version of SM2_L1 for export
SM2_L1 <- SM2_L1_raw |>
  select(
    "Feature Name" = feature,
    "Identified Name" = identified_name,
    "Isomer" = isomer,
    mz = feature_mz,
    "Library mz" = library_mz,
    "Mass Error (ppm)" = mz_error_ppm,
    "Retention Time" = feature_rt,
    "Library Retention Time" = library_rt,
    "Retention Time Error (sec)" = rt_error_sec,
    Column = library_mode,
    "Ion Mode" = ion_mode,
    Adduct = adduct,
    "Other Adducts Matched" = other_adducts_matched,
    "Annotated in Level 3" = annotated,
    "Annotated As" = annotated_as,
    Formula = formula,
    "Multi-Mode Detection" = MMD,
    "KEGG ID" = kegg_id,
    "HMDB ID" = hmdb_id,
    "Other Level 1 Matches" = other_L1_matches)
#- 11.2.3: Pull info from TFT_annot_key (level 3), rename cols
SM2_L3_raw <- TFT_annot_key |>
  arrange(mz, ion_mode) |>
  # Mark if also in confirmed key (Level 1)
  mutate(identified = if_else(feature %in% SM2_L1$`Feature Name`, "Y", "N")) |>
  # Create comma-delimited list of identified names for this feature from Level 1
  mutate(identified_as = map_chr(feature, ~{
    matches <- SM2_L1 |> filter(`Feature Name` == .x) |> pull(`Identified Name`)
    if(length(matches) > 0) paste(matches, collapse = ", ") else ""
  })) |>
  # Create comma-delimited list of OTHER identified names for this feature within Level 3
  mutate(other_L3_matches = map2_chr(feature, identified_name, ~{
    matches <- TFT_annot_key |> 
      filter(feature == .x, identified_name != .y) |> 
      pull(identified_name)
    if(length(matches) > 0) paste(matches, collapse = ", ") else ""
  }))
#- 11.2.4: Create clean version of SM2_L3 for export
SM2_L3 <- SM2_L3_raw |>
  select(
    "Feature Name" = feature,
    "Identified Name" = identified_name,
    "Isomer" = isomer,
    mz,
    "Retention Time" = rt,
    Column = Mode,
    "Ion Mode" = ion_mode,
    Adduct,
    Formula,
    "Multi-Mode Detection" = MMD,
    "Mean Intensity" = mean_intensity,
    "Identification Method" = identification_method,
    "Annotation Probability" = annotation_probability,
    "KEGG ID" = KEGG,
    "PubChem CID" = CID,
    SMILES,
    "InChIKey" = inchi_key,
    "HMDB ID" = HMDB,
    "Identified in Level 1" = identified,
    "Identified As" = identified_as,
    "Other Level 3 Matches" = other_L3_matches,
    Kingdom,
    Superclass,
    Class,
    Subclass,
    "Alternative Parent" = alternative_parent
  )
#- 11.2.5: Make small metadata joiner for SM2_L1
SM2_L1_meta <- SM2_L1_raw |>
  select(feature, identified_name, mz = feature_mz, rt = library_rt, column = library_mode, adduct, isomer, other_matches_in_level = other_L1_matches) |>
  group_by(feature) |>
  slice(1) |>
  ungroup() |>
  mutate("Identification Confidence" = "Level 1")
#- 11.2.6: Make small metadata joiner for SM2_L3
SM2_L3_meta <- SM2_L3_raw |>
  select(feature, identified_name, mz, rt, column = Mode, adduct = Adduct, isomer, other_matches_in_level = other_L3_matches) |>
  group_by(feature) |>
  slice(1) |>
  ungroup() |>
  mutate("Identification Confidence" = "Level 3")
#- 11.2.7: Join together; keep level 1 over level 3 for duplicates
SM2_meta <- bind_rows(SM2_L1_meta, SM2_L3_meta) |>
  distinct(feature, .keep_all = TRUE)
#- 11.2.8: Export as formatted excel sheet
create_SM2(SM2_L1, SM2_L3, "Outputs/Supplementary/Supplementary Material 2.xlsx")
#+ 11.3: Supplementary Material 3
#- 11.3.1: Combine 12h and 24h t-test results
SM3_ttest_combined <- targ_12PGD_ttest |>
  select(
    feature = Metabolite,
    p_value_12h = p_value,
    mean_difference_12h = mean_difference,
    mean_YPGD_12h = mean_YPGD,
    mean_NPGD_12h = mean_NPGD
  ) |>
  full_join(
    targ_24PGD_ttest |>
      select(
        feature = Metabolite,
        p_value_24h = p_value,
        mean_difference_24h = mean_difference,
        mean_YPGD_24h = mean_YPGD,
        mean_NPGD_24h = mean_NPGD
      ),
    by = "feature"
  ) |>
  # Calculate fold changes from mean differences (data already log2 transformed)
  # So fold change = 2^(mean_difference)
  mutate(
    fold_change_12h = 2^mean_difference_12h,
    fold_change_24h = 2^mean_difference_24h
  ) |>
  # Left join metadata
  left_join(SM2_meta, by = "feature") |>
  # Extract mz and rt directly from feature names (overwrite metadata values)
  mutate(
    # Parse feature names like "C18_85.0296_49.4" or "HILIC_86.06_22.8"
    feature_parts = str_split(feature, "_"),
    mz = as.numeric(map_chr(feature_parts, ~ .x[2])),
    rt = as.numeric(map_chr(feature_parts, ~ .x[3]))
  ) |>
  select(-feature_parts) |>  # Remove the temporary parsing column
  # Round numeric columns for display - use 7 decimal places to capture very small p-values
  mutate(
    p_value_12h = round(p_value_12h, 7),
    p_value_24h = round(p_value_24h, 7),
    fold_change_12h = round(fold_change_12h, 2),
    fold_change_24h = round(fold_change_24h, 2),
    mean_YPGD_12h = round(mean_YPGD_12h, 2),
    mean_NPGD_12h = round(mean_NPGD_12h, 2),
    mean_YPGD_24h = round(mean_YPGD_24h, 2),
    mean_NPGD_24h = round(mean_NPGD_24h, 2),
    mz = round(mz, 4)
    # rt - no rounding to preserve full precision
  ) |>
  # Sort by p-values
  arrange(p_value_12h, p_value_24h) |>
  # Select and rename columns for final output
  select(
    "P (12h)" = p_value_12h,
    "P (24h)" = p_value_24h,
    "Fold Change (12h)" = fold_change_12h,
    "Fold Change (24h)" = fold_change_24h,
    "Feature Name" = feature,
    "Identified Name" = identified_name,
    mz,
    "Retention Time" = rt,
    Column = column,
    Adduct = adduct,
    Isomer = isomer,
    "Identification Confidence",
    "Mean Yes PGD (12h)" = mean_YPGD_12h,
    "Mean No PGD (12h)" = mean_NPGD_12h,
    "Mean Yes PGD (24h)" = mean_YPGD_24h,
    "Mean No PGD (24h)" = mean_NPGD_24h,
    "Other Potential Matches" = other_matches_in_level
  )
#- 11.3.2: Pull LIMMA results to metabolite info
SM3_limma <- limma_targ$summary_final |>
  rename(feature = Metabolite) |>
  group_by(feature) |>
  fill(everything(), .direction = "downup") |> # fill chars + nums
  slice(1) |>
  ungroup() |>
  left_join(SM2_meta, by = "feature") |>
  # Extract mz and rt directly from feature names (overwrite metadata values)
  mutate(
    # Parse feature names like "C18_85.0296_49.4" or "HILIC_86.06_22.8"
    feature_parts = str_split(feature, "_"),
    mz = as.numeric(map_chr(feature_parts, ~ .x[2])),
    rt = as.numeric(map_chr(feature_parts, ~ .x[3]))
  ) |>
  select(-feature_parts) |>  # Remove the temporary parsing column
  arrange(p.value_interaction, p.value_group, p.value_time) |>
  # Round numeric columns for display
  mutate(
    # P-values: 3 decimals, but scientific if would round to 0.000
    p.value_interaction = round(p.value_interaction, 4),
    p.value_group = round(p.value_group, 4),
    p.value_time = round(p.value_time, 4),
    # Log2 fold changes: 2 decimals
    logFC_interaction = round(logFC_interaction, 2),
    logFC_group = round(logFC_group, 2), 
    logFC_time = round(logFC_time, 2),
    # Mean columns: 2 decimals
    Mean_N_12 = round(Mean_N_12, 2),
    Mean_N_24 = round(Mean_N_24, 2),
    Mean_Y_12 = round(Mean_Y_12, 2),
    Mean_Y_24 = round(Mean_Y_24, 2),
    # Other numeric columns: round as needed
    mz = round(mz, 4)
    # rt - removed rounding to preserve full precision 
  ) |>
  select(
      "P (PGD)" = p.value_group,
      "P (Time)" = p.value_time,
      "P (Interaction)" = p.value_interaction,
      "log2(FC PGD)" = logFC_group,
      "log2(FC Time)" = logFC_time,
      "log2(FC Interaction)" = logFC_interaction,
      "Feature Name" = feature,
      "Identified Name" = identified_name,
      mz,
      "Retention Time" = rt,
      Column = column,
      Adduct = adduct,
      Isomer = isomer,
      "Identification Confidence",
      "Mean No PGD (12h)" = Mean_N_12,
      "Mean No PGD (24h)" = Mean_N_24,
      "Mean Yes PGD (12h)" = Mean_Y_12,
      "Mean Yes PGD (24h)" = Mean_Y_24,
      "Other Potential Matches" = other_matches_in_level)
#- 11.3.3: Export as formatted excel sheet with conditional formatting
create_SM3(SM3_ttest_combined, SM3_limma, "Outputs/Supplementary/Supplementary Material 3.xlsx")
