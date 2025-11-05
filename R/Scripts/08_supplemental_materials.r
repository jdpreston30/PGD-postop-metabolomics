#* 8: Numbers for results section_columns
#+ 8.1: Supplementary Material 1
  #- 8.1.1: Make a flag for multi-mode detection and isomers, rename cols
    SM1 <- feature_key_isomark %>%
      # flag if KEGGID appears in >1 ion mode
      group_by(KEGGID) %>%
      mutate(multi_mode = if_else(n_distinct(ion_mode) > 1, "Y", "N")) %>%
      ungroup() %>%
      # flag isomers: exact same m/z AND RT within the same ion mode
      group_by(ion_mode, mz, time) %>%
      mutate(Isomer = if_else(n() > 1, "Y", "N")) %>%
      ungroup() %>%
      arrange(mz, ion_mode) %>%
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
  #- 8.1.2: Export
    write.xlsx(SM1, "Outputs/Supplemental/SM1.xlsx")
#+ 8.2: Supplementary Material 2
  #- 8.2.1: Pull info from SM1 needed to join
    SM1_info <- SM1 %>%
      select("Feature Name", mz, "Retention Time", "Exact Mass", Isomer, "Multi-Mode Detection")
  #- 8.2.2: Join to LIMMA targeted summary_final, rename columns, organize
    SM2 <- limma_targ$summary_final %>%
      rename("Feature Name" = Metabolite) %>%
      left_join(SM1_info, by = "Feature Name") %>%
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
        "Mean Yes PGD (24h)" = Mean_Y_24) %>%
      group_by(`Feature Name (KEGGID_Mode)`) %>%
      fill(everything(), .direction = "downup") %>% # fill chars + nums
      slice(1) %>%
      ungroup() %>%
      arrange(`P (Interaction)`, `P (PGD)`, `P (Time)`)
  #- 8.2.3: Export
    write.xlsx(SM2, "Outputs/Supplemental/SM2.xlsx")

