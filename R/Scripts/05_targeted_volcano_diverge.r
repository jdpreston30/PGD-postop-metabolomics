#* 5: Targeted Balloon and Volcano Plots
#+ 5.1: Set up data for calculating fold changes
targ12_PGD <- TFT_combined_LIMMA |>
  filter(Time == 12, severe_PGD == 'Y') |>
  select(-Patient, -Time, -severe_PGD, -Clinical_PGD)
targ12_NPGD <- TFT_combined_LIMMA |>
  filter(Time == 12, severe_PGD == 'N') |>
  select(-Patient, -Time, -severe_PGD, -Clinical_PGD)
targ24_PGD <- TFT_combined_LIMMA |>
  filter(Time == 24, severe_PGD == 'Y') |>
  select(-Patient, -Time, -severe_PGD, -Clinical_PGD)
targ24_NPGD <- TFT_combined_LIMMA |>
  filter(Time == 24, severe_PGD == 'N') |>
  select(-Patient, -Time, -severe_PGD, -Clinical_PGD)
#+ 5.2: Undo log2 transformation for the targeted data
targ_list <- list(targ12_PGD, targ12_NPGD, targ24_PGD, targ24_NPGD)
targ_list_raw <- lapply(targ_list, undo_log2_transform)
targ_list_final <- lapply(targ_list_raw, function(df) colMeans(df, na.rm = TRUE))
#+ 5.3: Calculate fold changes and conduct t tests for p-values
fc_targ12 <- data.frame(Value = log2(targ_list_final[[1]] / targ_list_final[[2]]))
fc_targ12$Time <- "12h"
fc_targ24 <- data.frame(Value = log2(targ_list_final[[3]] / targ_list_final[[4]]))
fc_targ24$Time <- "24h"
#+ 5.4: Conduct t-tests
ttest_result_12h <- lapply(seq_along(targ12_NPGD), function(i) {
    test <- t.test(targ12_NPGD[[i]], targ12_PGD[[i]])
    as_tibble(data.frame(
    Metabolite = colnames(targ12_PGD)[i],
    p_value = test$p.value,
    mean_difference = diff(test$estimate),
    mean_YPGD = mean(targ12_PGD[[i]]),
    mean_NPGD = mean(targ12_NPGD[[i]])
    ))
})
#+ 5.5: Combine results into a single table (12h)
targ_12PGD_ttest <- do.call(rbind, ttest_result_12h)
targ_12PGD_ttest$Time <- "12h"
#+ 5.6: 24h t test
ttest_result_24h <- lapply(seq_along(targ24_NPGD), function(i) {
    test <- t.test(targ24_NPGD[[i]], targ24_PGD[[i]])
    as_tibble(data.frame(
    Metabolite = colnames(targ24_PGD)[i],
    p_value = test$p.value,
    mean_difference = diff(test$estimate),
    mean_YPGD = mean(targ24_PGD[[i]]),
    mean_NPGD = mean(targ24_NPGD[[i]])
    ))
})
#+ 5.7: Combine results into a single table (24h)
targ_24PGD_ttest <- do.call(rbind, ttest_result_24h)
targ_24PGD_ttest$Time <- "24h"
#+ 5.8: Create volcano plots
#- 5.8.1: 12 Hours
volc_12 <- make_volcano(fc_targ12, targ_12PGD_ttest)
ggsave(
  "Outputs/Balloon_and_Volcano/volcano_plot_12h.png",
  volc_12$volcano_plot,
  width = 3.5,
  height = 3.5,
  units = "in",
  dpi = 600
)
#- 5.8.2: Create volcano plot based on volcano function
volc_24 <- make_volcano(fc_targ24, targ_24PGD_ttest)
ggsave(
  "Outputs/Balloon_and_Volcano/volcano_plot_24h.png",
  volc_24$volcano_plot,
  width = 3.5,
  height = 3.5,
  units = "in",
  dpi = 600
)
#+ 5.9: Prep data for diverging bars
#- 5.9.0: Subset to QC
TFT_QC_diverge <- TFT_QC |>
  filter(select_12h == "Y" | select_24h == "Y") |>
  select(display_name, feature, main_group, subgroup, select_12h, select_24h)
#- 5.9.1: Structure 12h data
targ_12PGD_diverging <- targ_12PGD_ttest |>
  select(feature = Metabolite, log2_fc = mean_difference, p_value)
#- 5.9.2: Subset the 12h
diverge_12h <- TFT_QC_diverge |>
  left_join(targ_12PGD_diverging, by = "feature") |>
  filter(select_12h == "Y") |>
  arrange(main_group, log2_fc)
#- 5.9.3: Structure 24h data
targ_24PGD_diverging <- targ_24PGD_ttest |>
  select(feature = Metabolite, log2_fc = mean_difference, p_value)
#- 5.9.4: Subset the 24h
diverge_24h <- TFT_QC_diverge |>
  left_join(targ_24PGD_diverging, by = "feature") |>
  filter(select_24h == "Y") |>
  arrange(main_group, log2_fc)
#- 5.9.1: 12 Hours
div_bars_12 <- plot_diverging_bars(diverge_12h, 
  group_ordering = TRUE, 
  add_group_labels = TRUE,
  max_features = 50, 
  fc_threshold = 0,
  x_max = 3.2
)
print_to_png(div_bars_12,
  "12h_test.png",
  width = 3.2,
  height = 5  
)
#- 5.9.1: 24 Hours
div_bars_24 <- plot_diverging_bars(diverge_24h, 
  group_ordering = TRUE, 
  add_group_labels = TRUE,
  max_features = 50,  
  fc_threshold = 0,
  x_max = 3.2
)
print_to_png(div_bars_24,
  "24h_test.png",
  width = 3.2,
  height = 5  
)
