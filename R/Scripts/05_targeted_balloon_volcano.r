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
#!!!!!!!!!!!!!!!!!!!!
# #+ 5.8: Create volcano and balloon plots for 12h
# #- 5.8.1: Create volcano plot based on volcano function
# volc_12 <- make_volcano(fc_targ12, targ_12PGD_ttest, combined_key_clean)
# ggsave(
#   "Outputs/Balloon_and_Volcano/volcano_plot_12h.png",
#   volc_12$volcano_plot,
#   width = 3.5,
#   height = 3.5,
#   units = "in",
#   dpi = 600
# )
# #- 5.8.2: Clean up data
# balloon_data_12 <- as_tibble(volc_12$volcano_data) |>
#   dplyr::mutate(
#     Identified_Name = ifelse(Identified_Name == "C20818",
#       "Carbapenem biosynthesis intermediate 2",
#       Identified_Name
#     ),
#     # Add asterisk to Identified_Name if Metabolite ends with "*"
#     Identified_Name = ifelse(stringr::str_detect(Metabolite, "\\*$"),
#       paste0(Identified_Name, "*"),
#       Identified_Name
#     )
#   ) |>
#   filter(Legend != "Not Significant")
# #- 5.8.3: Create a balloon plot based on volcano plot results
# balloon_12 <- balloon_plot(balloon_data_12)
# ggsave(
#   "Outputs/Balloon_and_Volcano/balloon_plot_12h.png",
#   balloon_12$ball_plot,
#   width = 5.5,
#   height = 8,
#   units = "in",
#   dpi = 600
# )
# #+ 5.9: Create volcano and balloon plots for 24h
# #- 5.9.1: Create volcano plot based on volcano function
# volc_24 <- make_volcano(fc_targ24, targ_24PGD_ttest, combined_key_clean)
# ggsave(
#   "Outputs/Balloon_and_Volcano/volcano_plot_24h.png",
#   volc_24$volcano_plot,
#   width = 3.5,
#   height = 3.5,
#   units = "in",
#   dpi = 600
# )
# #- 5.9.2: Touch up names
# balloon_data_24 <- as_tibble(volc_24$volcano_data) |>
#   dplyr::mutate(
#     Identified_Name = ifelse(Identified_Name == "C20818",
#       "Carbapenem biosynthesis intermediate 2",
#       Identified_Name
#     ),
#     # Add asterisk to Identified_Name if Metabolite ends with "*"
#     Identified_Name = ifelse(stringr::str_detect(Metabolite, "\\*$"),
#       paste0(Identified_Name, "*"),
#       Identified_Name
#     )
#   ) |>
#   filter(Legend != "Not Significant")
# #- 5.9.3: Create a balloon plot based on volcano plot results
# balloon_24 <- balloon_plot(balloon_data_24)
# ggsave(
#   "Outputs/Balloon_and_Volcano/balloon_plot_24h.png",
#   balloon_24$ball_plot,
#   width = 5.5,
#   height = 8,
#   units = "in",
#   dpi = 600
# )