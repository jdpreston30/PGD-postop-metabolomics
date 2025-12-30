#+++++++++++++
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
