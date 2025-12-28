#* 9: Numbers for results section_columns
#+ 9.1: Supplementary Material 1
#- 9.1.1: Make a flag for multi-mode detection and isomers, rename cols
SM1 <- feature_key_isomark |>
  # flag if KEGGID appears in >1 ion mode
  group_by(KEGGID) |>
  mutate(multi_mode = if_else(n_distinct(ion_mode) > 1, "Y", "N")) |>
  ungroup() |>
  # flag isomers: exact same m/z AND RT within the same ion mode
  group_by(ion_mode, mz, time) |>
  mutate(Isomer = if_else(n() > 1, "Y", "N")) |>
  ungroup() |>
  arrange(mz, ion_mode) |>
  select(
    "Feature Name" = Rename,
    "Identified Name" = Identified_Name,
    Isomer,
    mz,
    "Retention Time" = time,
    "Ion Mode" = ion_mode,
    Adduct,
    "Exact Mass" = Exact_mass,
    Formula,
    "Multi-Mode Detection" = multi_mode,
    "Mean Intensity" = mean_intensity,
    "Identification Method" = identification_method,
    "Annotation Probability" = Probability,
    "KEGG ID" = KEGGID,
    "PubChem CID" = PUBCHEM_CID,
    SMILES,
    InChIKey,
    "HMDB ID" = HMDBID,
    Kingdom,
    Superclass,
    Class,
    Subclass,
    "Alternative Parent" = Alternative_parent)
#- 9.1.2: Export
write.xlsx(SM1, "Outputs/Supplemental/SM1.xlsx")
#+ 9.2: Supplementary Material 2
#- 9.2.1: Pull info from SM1 needed to join
SM1_info <- SM1 |>
  select("Feature Name", mz, "Retention Time", "Exact Mass", Isomer, "Multi-Mode Detection")
#- 9.2.2: Join to LIMMA targeted summary_final, rename columns, organize
SM2 <- limma_targ$summary_final |>
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
#- 9.2.3: Export
write.xlsx(SM2, "Outputs/Supplemental/SM2.xlsx")




#- 7.3.1: Organize data for balloon
balloon_data_export <- balloon_data_24 |>
  rbind(balloon_data_12) |>
  arrange(Time, Value) |>
  select(Time, Identified_Name, Metabolite, Value, p_value, mean_YPGD, mean_NPGD, mean_difference, Legend)
#- 7.3.2: Save data for balloon
write_csv(balloon_data_export, "Outputs/LIMMA/balloon_data_export.csv")
#- 7.3.3: Organize interaction data
interactions_data <- limma_targ$interaction |>
  arrange(logFC) |>
  filter(p.value < 0.05) |>
  select(Metabolite, Identified_Name, logFC, p.value) |>
  # drop explicit duplicates you don't want
  # C16309_C18 is detected in HILIC and stronger magnitude
  filter(!Metabolite %in% c("C11384_HILIC", "C11383_HILIC", "C20893_C18", "C16309_C18")) |>
  # start from the given name (fallback to ID if missing)
  mutate(label = coalesce(na_if(Identified_Name, ""), Metabolite)) |>
  mutate(
    label = dplyr::case_when(
      Metabolite == "C00079_HILIC" ~ "Phenylalanine",
      Metabolite == "C20892_C18*" ~ "3-(5-oxoisoxazolin-2-yl)-alanine",
      Metabolite == "C09848_HILIC*" ~ "Citronellal",
      Metabolite == "C01767_HILIC*" ~ "Carvone",
      Metabolite == "C22300_C18" ~ "3-Hydroxyisoleucine", # keep one version
      Metabolite == "C10462_HILIC" ~ "Gingerol",
      Metabolite == "C21087_HILIC" ~ "Epoxy-4S-H4HPP",
      Metabolite == "C20818_C18" ~ "Carbapenem biosynthesis intermediate 2",
      Metabolite == "C20324_HILIC" ~ "4-OH-TMCP acetate", # short name you chose
      TRUE ~ label
    )
  ) |>
  # append exactly one * if metabolite ID ends with * and label lacks one, then strip all * per your preference
  mutate(
    label = if_else(stringr::str_detect(Metabolite, "\\*$") & !stringr::str_detect(label, "\\*$"),
      paste0(label, "*"), label
    ),
    label = stringr::str_replace_all(label, "\\*", "")
  ) |>
  # FINAL ordering and factor levels (do this LAST and don’t touch `label` afterwards)
  arrange((logFC)) |> # top = most positive
  mutate(
    dir   = if_else(logFC >= 0, "Positive", "Negative"),
    label = factor(label, levels = label) # lock order
  )
#- 7.3.4: Subset feature table by these and compute means per group
TFT_interactions <- TFT |>
  select(Clinical_PGD, Time, all_of(intersect(interactions_data$Metabolite, names(TFT)))) |>
  pivot_longer(
    cols = -c(Clinical_PGD, Time),
    names_to = "Metabolite",
    values_to = "val"
  ) |>
    mutate(Group = paste0(Clinical_PGD, "_PGD_", Time)) |>
    group_by(Metabolite, Group) |>
    summarise(mean = mean(val, na.rm = TRUE), .groups = "drop") |>
    mutate(Group = factor(Group, levels = c("N_PGD_12", "Y_PGD_12", "N_PGD_24", "Y_PGD_24"))) |>
    arrange(Metabolite, Group) |>
    pivot_wider(names_from = Group, values_from = mean) |>
    select(Metabolite, N_PGD_12, Y_PGD_12, N_PGD_24, Y_PGD_24)
#- 7.3.5: Join means to interaction data
interactions_data_means <- interactions_data |>
  select(Metabolite, Identified_Name, logFC, p.value) |>
  left_join(TFT_interactions, by = "Metabolite")
#- 7.3.6: Save interaction data
write_csv(interactions_data_means, "Outputs/LIMMA/interactions_data_means.csv")
