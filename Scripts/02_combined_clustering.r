#* 2: Combined Clustering analysis
#+ 2.1: Run Exploratory PLSDA
  #- 2.1.1: Define datasets
    UFT_12h <- UFT_nomiss %>% 
      filter(Time == 12) %>%
      select(Sample_ID, Time, Clinical_PGD, all_of(untargeted_features))
    UFT_24h <- UFT_nomiss %>%
      filter(Time == 24) %>%
      select(Sample_ID, Time, Clinical_PGD, all_of(untargeted_features))
    UFT_12and24 <- UFT_nomiss %>%
      select(Sample_ID, Time, Clinical_PGD, all_of(untargeted_features))
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
        group_var = "Clinical_PGD", method = ..4
      )),
      .y = plot_specs$outpath,
      ~ ggsave(.y, .x$plot, width = 5, height = 5, units = "in", dpi = 600)
    )
#!!!!!!!!!!!!
#+ 2.2: Create heatmaps
  #- 2.2.1: Prepare data for heatmap
    heatmap_data_12h <- UFT_12h %>%
    select(Sample_ID, Clinical_PGD, Time, any_of(untargeted_features))
    heatmap_data_24h <- UFT_24h %>%
    select(Sample_ID, Clinical_PGD, Time, any_of(untargeted_features))
    heatmap_data_combTime <- UFT_12and24 %>%
    select(Sample_ID, Clinical_PGD, Time, any_of(untargeted_features))
  #- 2.2.2: Create heatmaps
    # Set file path for saving heatmaps
    save_path <- "C://Users//amshu//Desktop//PGD_time_course-1//Outputs//Heatmaps//"
    hmap_pgd_12 <- make_heatmap(
      heatmap_data_12h,
      feature_selector = "ttest",
      top_features = 250,
      file_path = save_path,
      file_name = '12h_heatmap'
      )

    hmap_pgd_24 <- make_heatmap(
      heatmap_data_24h,
      feature_selector = "ttest",
      top_features = 250,
      file_path = save_path,
      file_name = '24h_heatmap'
      )
  variance_combined_1000 <- make_heatmap(
    heatmap_data_combTime,
    feature_selector = "variance",
    top_features = FALSE,
    filename = "comb_variance"
  )

#+ 3.4: Run PERMANOVA analysis
  #- 3.4.1: Define feature columns
    permanova_features <- rownames(variance_combined_1000$M)
  #- 3.4.1: Extract features and prepare data
    variance_study <- combined_UFT %>%
      left_join(variance_combined_1000$cluster_df, by = "Sample_ID") %>%
      select(Sample_ID, Time, Age, Sex, Clinical_PGD, Patient, any_of(permanova_features)) %>%
      mutate(
        across(c(Sex, Time, Clinical_PGD, Age), as.factor),
        Age = as.numeric(Age),
        Clinical_PGD = factor(Clinical_PGD, levels = c("Y", "N")),
        Time = factor(Time, levels = c("12", "24")),
        Sex = factor(Sex, levels = c("M", "F"))
      )
  #- 3.4.2: Define PERMANOVA variables
    permanova_variables <-  c("Time", "Clinical_PGD", "Age", "Sex")
  #- 3.4.4: Extract feature data for analysis
    features_1000 <- variance_study %>%
      select(any_of(permanova_features))
  #- 3.4.5: Prepare metadata for individual tests
    meta_use <- variance_study %>%
      select(all_of(permanova_variables)) %>%
      mutate(
        across(c(Time, Clinical_PGD, Sex), as.factor),
        Age = as.numeric(Age)
      )
  #- 3.4.6: Run PERMANOVA for each variable
    permanova_results_1000 <- bind_rows(lapply(
      permanova_variables, 
      get_permanova,
      feature_data = features_1000,
      meta_data = meta_use,
      ctrl = permute::how(nperm = 9999),
      seed = 2025)) %>%
      arrange(p_value) %>%
      mutate(Variable = if_else(Variable == "T_stage", "T stage", Variable))
  #- 3.4.7: Create PERMANOVA visualiation data
    permanova_viz <- permanova_results_1000 %>%
      mutate(
        Significance = case_when(
          p_value < 0.001 ~ "p < 0.001",
          p_value < 0.01 ~ "p < 0.01", 
          p_value < 0.05 ~ "p < 0.05",
          p_value < 0.1 ~ "p < 0.1",
          TRUE ~ "n.s."
        ),
        Significance = factor(Significance, levels = c("p < 0.001", "p < 0.01", "p < 0.05", "p < 0.1", "n.s.")),
        p_label = case_when(
          p_value < 0.001 ~ "p < 0.001",
          p_value < 0.05 ~ paste0("p = ", round(p_value, 3)),
          TRUE ~ paste0("p = ", round(p_value, 2))
        )
      )
#- 3.4.8: Create PERMANOVA bar plot (vertical with horizontal bars)
  permanova_1B <- ggplot(permanova_viz, aes(y = reorder(Variable, -p_value), x = R2)) +
    geom_col(
      width = 0.72,
      fill = "black",
      color = "black",
      alpha = 1,
      linewidth = 0,
      na.rm = TRUE
    ) +
    geom_text(aes(label = p_label), hjust = -0.1, vjust = 0.5, size = 3.0, fontface = "bold") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
    coord_cartesian(xlim = c(0, 0.7), clip = "off") +
    theme_pub_simple() +
    labs(y = NULL, x = expression(bold("R"^2 * " (Explained Variance)")))
    ggsave(
      filename = "./Permanova.png",
      plot = last_plot(),
      device = 'png',
      width = 4,
      height = 4,
      unit = 'in',
      dpi = 600
    )
