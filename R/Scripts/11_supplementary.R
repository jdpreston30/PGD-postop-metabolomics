#* 9: Numbers for results section_columns
#+ 9.1: Supplementary Material 1
#! This was a word document of the supplementary methods
#+ 9.2: Supplementary Material 2
#- 9.2.2: Structure confirmed key
SM2_L1 <- TFT_confirmed_key |>
  # Add the column type
  mutate(ion_mode = if_else(library_mode == "HILIC", "Positive", "Negative")) |>
  arrange(feature_mz, ion_mode) |>
  # Mark if also in annotated key
  mutate(annotated = if_else(feature %in% TFT_annot_key$feature, "Y", "N")) |>
  # Create comma-delimited list of annotation names for this feature
  mutate(annotated_as = map_chr(feature, ~{
    matches <- TFT_annot_key |> filter(feature == .x) |> pull(identified_name)
    if(length(matches) > 0) paste(matches, collapse = ", ") else ""
  })) |>
  # Use function to clean up names
  mutate(identified_name = normalize_chem_names(identified_name)) |>
  # Clean up isomer column
  mutate(isomer = case_when(
    is.na(isomer) ~ "N",
    isomer == "yes" ~ "Y", 
    str_detect(isomer, "yes -but different (time|RT)") ~ "Y, but different RT",
    TRUE ~ isomer
  )) |>
  # Select final columns
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
    "HMDB ID" = hmdb_id)
#- 9.2.2: Pull info from TFT_annot_key (level 3), rename cols
SM2_L3 <- TFT_annot_key |>
  arrange(mz, ion_mode) |>
  # Mark if also in confirmed key (Level 1)
  mutate(identified = if_else(feature %in% SM2_L1$`Feature Name`, "Y", "N")) |>
  # Create comma-delimited list of identified names for this feature from Level 1
  mutate(identified_as = map_chr(feature, ~{
    matches <- SM2_L1 |> filter(`Feature Name` == .x) |> pull(`Identified Name`)
    if(length(matches) > 0) paste(matches, collapse = ", ") else ""
  })) |>
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
    Kingdom,
    Superclass,
    Class,
    Subclass,
    "Alternative Parent" = alternative_parent
  )
#- 9.2.2: Export as formatted excel sheet
create_SM2(SM2_L1, SM2_L3, "Outputs/Supplementary/Supplementary Material 2.xlsx")
#+ 9.3: Supplementary Material 3
#- 9.3.2: Join to LIMMA targeted summary_final, rename columns, organize
SM3 <- limma_targ$summary_final |>
  rename("Feature Name" = Metabolite) |>
  left_join(SM1_info, by = "Feature Name") |>
  select(
    "Feature Name (KEGGID_Mode)" = "Feature Name",
    "Identified Name" = Identified_Name,
    mz,
    "Retention Time",
    "Exact Mass",
    Isomer,
    "Multi-Mode Detection",
    "P (Interaction)" = p.value_interaction,
    "P (PGD)" = p.value_group,
    "P (Time)" = p.value_time,
    "log2(FC Interaction)" = logFC_interaction,
    "log2(FC PGD)" = logFC_group,
    "log2(FC Time)" = logFC_time,
    "Mean No PGD (12h)" = Mean_N_12,
    "Mean No PGD (24h)" = Mean_N_24,
    "Mean Yes PGD (12h)" = Mean_Y_12,
    "Mean Yes PGD (24h)" = Mean_Y_24) |>
  group_by(`Feature Name (KEGGID_Mode)`) |>
  fill(everything(), .direction = "downup") |> # fill chars + nums
  slice(1) |>
  ungroup() |>
  arrange(`P (Interaction)`, `P (PGD)`, `P (Time)`)
#- 9.3.3: Export as formatted excel sheet
create_SM3(SM3, "Outputs/Supplementary/Supplementary Material 23.xlsx")