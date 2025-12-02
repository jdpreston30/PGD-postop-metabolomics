# Pull 12 and 24h targeted t-test results for PGD from volcano analysis
ttest_12h_for_export <- targ_12PGD_ttest %>%
  mutate(fold_change_12h = log2(mean_YPGD / mean_NPGD)) %>%
  select(Metabolite,
    p_value_12h = p_value,
    mean_YPGD_12h = mean_YPGD,
    mean_NPGD_12h = mean_NPGD,
    fold_change_12h
  )
ttest_24h_for_export <- targ_24PGD_ttest %>%
  mutate(fold_change_24h = log2(mean_YPGD / mean_NPGD)) %>%
  select(Metabolite,
    p_value_24h = p_value,
    mean_YPGD_24h = mean_YPGD,
    mean_NPGD_24h = mean_NPGD,
    fold_change_24h
  )
# Combine and join with name
ttest_combined_export <- ttest_12h_for_export %>%
  full_join(ttest_24h_for_export, by = "Metabolite") %>%
  left_join(combined_key_clean, by = "Metabolite") %>%
  select(
    Metabolite, Identified_Name,
    p_value_12h, mean_YPGD_12h, mean_NPGD_12h, fold_change_12h,
    p_value_24h, mean_YPGD_24h, mean_NPGD_24h, fold_change_24h
  ) %>%
  mutate(Identified_Name = str_replace_all(Identified_Name, "\u1D62", "ⁱ")) %>%
  arrange(Metabolite)

write_xlsx(ttest_combined_export, "Outputs/Balloon_and_Volcano/targeted_ttest_results.xlsx")