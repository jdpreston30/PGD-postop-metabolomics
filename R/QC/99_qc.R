#* Cleaning up QC
#+ Bring in HMDB xenobiotics and KEGG drug metabolites
hmdb_xeno <- read_csv("Databases/QC/hmdb_metabolites_xenobiotics.csv") |>
  select(HMDB = HMDBID, KEGG = KEGG_ID)
kegg_pharm <- read_csv("Databases/QC/KEGG_drug_metabolite.csv") |>
  select(KEGG = metabolite_KEGG_ID)
#+ Clean up the combined key for relevant info
#- Clean up relevant confirmed key
conf_key_pared <- TFT_confirmed_key |>
  filter(!is.na(identified_name)) |>
  mutate(
    source = "library",
    lib_conf = "Y",
    # Fix capitalization if all caps
    identified_name = if_else(
      str_detect(identified_name, "^[A-Z0-9 ,'()\\-]+$") & str_detect(identified_name, "[A-Z]"),
      str_to_title(identified_name),
      identified_name
    ),
    isomer = case_when(
      is.na(isomer) ~ "N",
      tolower(isomer) == "yes" ~ "Y",
      TRUE ~ "Y-diff RT"
    )
  ) |>
  select(identified_name, feature, isomer, ref_mz = library_mz, ref_rt = library_rt, adduct, KEGG = kegg_id, source, lib_conf)
#- Clean up relevant annot key
annot_key_pared <- TFT_annot_key |>
  filter(!is.na(identified_name)) |>
  mutate(source = "annotation", lib_conf = "N") |>
  select(identified_name, feature, isomer, ref_mz = mz, ref_rt = rt, adduct = Adduct, KEGG, HMDB, source, lib_conf)
#- Combine the two keys; add symbols
combined_key_pared <- bind_rows(conf_key_pared, annot_key_pared) |>
  # Check if feature appears in both sources
  group_by(feature) |>
  mutate(
    lib_and_annot = if_else(n_distinct(source) > 1, "Y", "N"),
    # Check if feature maps to multiple unique KEGG IDs
    # BUT if lib_and_annot = Y with only 2 entries (library + annotation) and same KEGG, mark as N
    potential_isomer = case_when(
      # If lib_and_annot = Y, only 2 rows, and KEGG IDs match, then N
      lib_and_annot == "Y" & n() == 2 & n_distinct(KEGG[!is.na(KEGG)]) <= 1 ~ "N",
      # Otherwise apply original logic
      n_distinct(KEGG[!is.na(KEGG)]) > 1 | any(isomer %in% c("Y", "Y-diff RT")) ~ "Y",
      TRUE ~ "N"
    )
  ) |>
  ungroup() |>
  # Add xenobiotic and pharmaceutical flags
  mutate(
    hmdb_xeno_flag = if_else(
      (!is.na(KEGG) & KEGG %in% hmdb_xeno$KEGG) | (!is.na(HMDB) & HMDB %in% hmdb_xeno$HMDB),
      "Y", "N"
    ),
    kegg_pharm_flag = if_else(!is.na(KEGG) & KEGG %in% kegg_pharm$KEGG, "Y", "N")
  ) |>
  # Add symbols based on isomer and library confirmation
  mutate(
    is_isomer = isomer %in% c("Y", "Y-diff RT"),
    Identified_Name = case_when(
      KEGG == "C20913" ~ "Tabtoxin biosynthesis intermediate 3",
      KEGG == "C20915" ~ "Tabtoxin biosynthesis intermediate 5",
      TRUE ~ identified_name
    ),
    Identified_Name = if_else(is_isomer, paste0(Identified_Name, "*"), Identified_Name),
    Identified_Name = if_else(lib_conf == "Y", paste0(Identified_Name, "ⁱ"), Identified_Name)
  ) |>
  select(Identified_Name, KEGG, HMDB, everything(), -identified_name, -is_isomer, -lib_conf) |>
  arrange(feature) |>
  relocate(lib_and_annot, .after = source) |>
  relocate(potential_isomer, .after = isomer) |>
  relocate(hmdb_xeno_flag, kegg_pharm_flag, .after = potential_isomer)
#+ Joining with p-values
#- Pull 12 and 24h targeted t-test results for PGD from volcano analysis
ttest_12h_for_export <- targ_12PGD_ttest |>
  select(feature = Metabolite, p_value_12h = p_value, log2FC_YN_12h = mean_difference)
ttest_24h_for_export <- targ_24PGD_ttest |>
  select(feature = Metabolite, p_value_24h = p_value, log2FC_YN_24h = mean_difference)
#- Join 12h and 24h results by feature
ttest_combined <- ttest_12h_for_export |>
  full_join(ttest_24h_for_export, by = "feature")
#- Assign p values and FCs to all rows in combined_key_pared
combined_key_with_stats <- combined_key_pared |>
  left_join(ttest_combined, by = "feature")
#- Filter to significant features (p <= 0.05 at either timepoint)
combined_key_with_stats_sig <- combined_key_with_stats |>
  filter(p_value_12h <= 0.05 | p_value_24h <= 0.05)
#- Create coalesced version: prioritize library, collapse duplicates, capture alternate names
combined_key_with_stats_sig_coalesced <- combined_key_with_stats_sig |>
  group_by(feature) |>
  mutate(
    # Capture annotation names when both sources exist
    annot_names = if_else(
      lib_and_annot == "Y",
      paste(Identified_Name[source == "annotation"], collapse = "; "),
      NA_character_
    ),
    # Capture isomer names (all identified names except the first)
    isomer_names = if_else(
      any(isomer %in% c("Y", "Y-diff RT")) & n() > 1,
      paste(Identified_Name[-1], collapse = "; "),
      NA_character_
    )
  ) |>
  # Keep only: library version if lib_and_annot=Y, otherwise first row
  filter(
    (lib_and_annot == "Y" & source == "library") | 
    (lib_and_annot == "N" & row_number() == 1)
  ) |>
  ungroup()
#- Write to excel spreadsheets
write_xlsx(combined_key_with_stats, "combined_key_with_stats.xlsx")
write_xlsx(combined_key_with_stats_sig_coalesced, "combined_key_with_stats_sig_coalesced.xlsx")
