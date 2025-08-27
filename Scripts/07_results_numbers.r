#* 7: Numbers for results section_columns
#+ 7.1: Untargeted paragraph
  #- 7.1.1: Count features in complete UFT
    hilic_count <- C18_HILIC_unt_proc %>%
      select(starts_with("HILIC")) %>%
      ncol()
    c18_count <- C18_HILIC_unt_proc %>%
      select(starts_with("C18")) %>%
      ncol()
    combined_count <- hilic_count + c18_count
  #- 7.1.2: Remaining features following elimination of 20% detection
    UFT_filtered_count <- UFT %>%
      select(starts_with("C18"), starts_with("HILIC")) %>%
      ncol()
  #- 7.1.3: LIMMA significant
    sig_limma_untarg <- tibble(
      contrast = c("interaction", "time_main", "group_main"),
      n_sig = c(
        sum(limma_untarg$interaction$p.value < 0.05, na.rm = TRUE),
        sum(limma_untarg$time_main$p.value < 0.05, na.rm = TRUE),
        sum(limma_untarg$group_main$p.value < 0.05, na.rm = TRUE)
      )
    )
#+ 7.2: Untargeted paragraph
  #- 7.2.1: Count features in complete TFT
    HILIC_count_targ <- HILIC_targeted_raw %>%
      select(-Sample_ID) %>%
      ncol()
    C18_count_targ <- C18_targeted_raw %>%
      select(-Sample_ID) %>%
      ncol()
    combined_count_targ <- HILIC_count_targ + C18_count_targ
  #- 7.1.2: Remaining features following elimination of 20% detection
    targ_post_elim <- TFT %>%
      select(-meta_cols) %>%
      ncol()
  #- 7.1.3: Detection within both modes based on KEGG
    both_modes <- feature_key %>%
      group_by(KEGGID) %>%
      summarise(n_modes = n_distinct(ion_mode), .groups = "drop") %>%
      filter(n_modes == 2)
    nrow(both_modes)
  #- 7.1.4: Number of isomers
    collapsed_all_isomers_counts <- collapsed_all_isomers %>%
      filter(n_candidates > 1) %>%
      select(Mode, mz_time) %>%
      unique() %>%
      nrow()
  #- 7.1.5: Volcano counts (Up)
    balloon_data_12_up <- balloon_data_12 %>%
      filter(Legend == "Up in PGD") %>%
      nrow()
    balloon_data_24_up <- balloon_data_24 %>%
      filter(Legend == "Up in PGD") %>%
      nrow()
  #- 7.1.5: Volcano counts (Down)
    balloon_data_12_down <- balloon_data_12 %>%
      filter(Legend == "Down in PGD") %>%
      nrow()
    balloon_data_24_down <- balloon_data_24 %>%
      filter(Legend == "Down in PGD") %>%
      nrow()
  #- 7.1.5: LIMMA significant
    sig_limma_targ <- tibble(
      contrast = c("interaction", "time_main", "group_main"),
      n_sig = c(
        sum(limma_targ$interaction$p.value < 0.05, na.rm = TRUE),
        sum(limma_targ$time_main$p.value < 0.05, na.rm = TRUE),
        sum(limma_targ$group_main$p.value < 0.05, na.rm = TRUE)
      )
    )
#+ 7.3: Exports for Clay
  #- 7.3.1: Organize data for balloon
    balloon_data_export <- balloon_data_24 %>%
      rbind(balloon_data_12) %>%
      arrange(Time, Value) %>%
      select(Time, Identified_Name, Metabolite, Value, p_value, mean_YPGD, mean_NPGD, mean_difference, Legend)
  #- 7.3.2: Save data for balloon
    write_csv(balloon_data_export, "Outputs/LIMMA/balloon_data_export.csv")
  #- 7.3.3: Organize interaction data
    interactions_data <- limma_targ$interaction %>%
      arrange(logFC) %>%
      filter(p.value < 0.05) %>%
      select(Metabolite, Identified_Name, logFC, p.value) %>%
      # drop explicit duplicates you don't want
      # C16309_C18 is detected in HILIC and stronger magnitude
      filter(!Metabolite %in% c("C11384_HILIC", "C11383_HILIC", "C20893_C18", "C16309_C18")) %>%
      # start from the given name (fallback to ID if missing)
      mutate(label = coalesce(na_if(Identified_Name, ""), Metabolite)) %>%
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
      ) %>%
      # append exactly one * if metabolite ID ends with * and label lacks one, then strip all * per your preference
      mutate(
        label = if_else(stringr::str_detect(Metabolite, "\\*$") & !stringr::str_detect(label, "\\*$"),
          paste0(label, "*"), label
        ),
        label = stringr::str_replace_all(label, "\\*", "")
      ) %>%
      # FINAL ordering and factor levels (do this LAST and don’t touch `label` afterwards)
      arrange((logFC)) %>% # top = most positive
      mutate(
        dir   = if_else(logFC >= 0, "Positive", "Negative"),
        label = factor(label, levels = label) # lock order
      )
  #- 7.3.4: Subset feature table by these and compute means per group
    TFT_interactions <- TFT %>%
      select(Clinical_PGD, Time, all_of(intersect(interactions_data$Metabolite, names(TFT)))) %>%
      pivot_longer(
        cols = -c(Clinical_PGD, Time),
        names_to = "Metabolite",
        values_to = "val"
      ) %>%
        mutate(Group = paste0(Clinical_PGD, "_PGD_", Time)) %>%
        group_by(Metabolite, Group) %>%
        summarise(mean = mean(val, na.rm = TRUE), .groups = "drop") %>%
        mutate(Group = factor(Group, levels = c("N_PGD_12", "Y_PGD_12", "N_PGD_24", "Y_PGD_24"))) %>%
        arrange(Metabolite, Group) %>%
        pivot_wider(names_from = Group, values_from = mean) %>%
        select(Metabolite, N_PGD_12, Y_PGD_12, N_PGD_24, Y_PGD_24)
  #- 7.3.5: Join means to interaction data
    interactions_data_means <- interactions_data %>%
      select(Metabolite, Identified_Name, logFC, p.value) %>%
      left_join(TFT_interactions, by = "Metabolite")
  #- 7.3.6: Save interaction data
    write_csv(interactions_data_means, "Outputs/LIMMA/interactions_data_means.csv")
