#*2: PCA, PLSDA, and Heatmaps
#+ 2.1: Run Exploratory PLSDA
  #- 2.1.1: Define datasets
    UFT_12h <- UFT |> 
      filter(Time == 12) |>
      select(Sample_ID, Time, severe_PGD, all_of(untargeted_features))
    UFT_24h <- UFT |>
      filter(Time == 24) |>
      select(Sample_ID, Time, severe_PGD, all_of(untargeted_features))
    UFT_12and24 <- UFT |>
      select(Sample_ID, Time, severe_PGD, all_of(untargeted_features))
  #- 2.1.2: Define colors
    cluster_colors <- c(
      "Y" = "#94001E",
      "N" = "#03507D"
    )
  #- 2.1.3: Set Plot Specs
    plot_specs <- tribble(
      ~data, ~comp_x, ~comp_y, ~method, ~outpath,
      UFT_12h, 1, 2, "PLSDA", "Outputs/PLSDA/PLSDA_12h_ComparePGD.png",
      UFT_24h, 1, 2, "PLSDA", "Outputs/PLSDA/PLSDA_24h_ComparePGD.png",
      UFT_12and24, 1, 2, "PLSDA", "Outputs/PLSDA/PLSDA_combinedTime_ComparePGD.png",
      UFT_12h, 1, 2, "PCA", "Outputs/PCA/PCA_12h_ComparePGD.png",
      UFT_24h, 1, 2, "PCA", "Outputs/PCA/PCA_24h_ComparePGD.png",
      UFT_12and24, 1, 2, "PCA", "Outputs/PCA/PCA_combinedTime_ComparePGD.png"
    )
  #- 2.1.4: Run all the PLSDAs and PCAs in one sweep
    walk2(
      .x = pmap(plot_specs, ~ make_PCA(..1,
        comp_x = ..2, comp_y = ..3,
        group_var = "severe_PGD", method = ..4
      )),
      .y = plot_specs$outpath,
      ~ ggsave(.y, .x$plot, width = 5, height = 5, units = "in", dpi = 600, bg = "white")
    )
#+ 2.2: Create heatmaps
# #- 2.2.1: Create 12h heatmap (ttest 500)
# hmap_pgd_12 <- make_heatmap(
#   UFT_12h,
#   feature_selector = "variance",
#   top_features = 500,
#   file_path = "Outputs/Heatmaps/",
#   file_name = '12h_heatmap_500_ttest'
#   )
# #- 2.2.2: Create 24h heatmap (ttest 500)
# hmap_pgd_24 <- make_heatmap(
#   UFT_24h,
#   feature_selector = "variance",
#   top_features = 500,
#   file_path = "Outputs/Heatmaps/",
#   file_name = '24h_heatmap_500_ttest'
#   )
# #- 2.2.4: Create combined heatmap (ttest 500)
# hmap_combined <- make_heatmap(
#   UFT_12and24,
#   feature_selector = "variance",
#   top_features = 1000,
#   file_path = "Outputs/Heatmaps/",
#   file_name = 'combined_heatmap_500_ttest'
#   )
