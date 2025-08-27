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
  #- 7.3.1: Organize data
    balloon_data_export <- balloon_data_24 %>%
      rbind(balloon_data_12) %>%
      arrange(Time, Value) %>%
      select(Time, Identified_Name, Metabolite, Value, p_value, mean_YPGD, mean_NPGD, mean_difference, Legend)
  #- 7.3.2: Save data
    write_csv(balloon_data_export, "Outputs/LIMMA/balloon_data_export.csv")
