#* 3: Combined Clustering analysis
#+ 3.1: Run Exploratory PCAs
  #- 3.1.1: Define datasets
   cluster_colors <- c(
      "Y" = "#94001E",
      "N" = "#03507D"
    )
    str(UFT_12h)
    combined_UFT$Clinical_PGD[combined_UFT$Clinical_PGD == "Y'"] <- 'Y'
    UFT_12h <- combined_UFT %>% 
      filter(Time == 12) %>%
      select(Sample_ID, Time, Clinical_PGD, all_of(features_to_keep))
    UFT_24h <- combined_UFT %>%
      filter(Time == 24) %>%
      select(Sample_ID, Time, Clinical_PGD, all_of(features_to_keep))
    UFT_12and24 <- combined_UFT %>%
      select(Sample_ID, Time, Clinical_PGD, all_of(features_to_keep))
    UFT_YPGD <- combined_UFT %>%
      filter(Clinical_PGD == 'Y') %>%
      select(Sample_ID, Time, all_of(features_to_keep))
    UFT_NPGD <- combined_UFT %>%
      filter(Clinical_PGD == 'N') %>%
      select(Sample_ID, Time, all_of(features_to_keep))
  #- 3.1.2: Create 12h only PCAs
    pca_12only <- make_PCA(
      data = UFT_12h,
      comp_x = , comp_y = 2, group_var = "Clinical_PGD", method = "PLSDA"
    )
    ggplot2::ggsave(paste0("C://Users//amshu//OneDrive - Emory//Preston, Joshua's files - Amshu Josh PGD//Final Analysis//PLSDA_12h_ComparePGD.png"), pca_12only$plot, width = 5, height = 5, units = "in", dpi = 600)
    
  #- 3.1.3: Create T-24h only PCAs
    pca_24only <- make_PCA(
      data = UFT_24h,
      comp_x = 1, comp_y = 2, group_var = "Clinical_PGD", method = "PLSDA"
    )
    ggplot2::ggsave(paste0("C://Users//amshu//OneDrive - Emory//Preston, Joshua's files - Amshu Josh PGD//Final Analysis//PLSDA_24h_ComparePGD.png"), pca_24only$plot, width = 5, height = 5, units = "in", dpi = 600)

  #- 3.1.4: Create combined 12h and 24h PCAs
    pca_time_combined <- make_PCA(
      data = UFT_12and24,
      comp_x = 1, comp_y = 2, group_var = "Clinical_PGD", method = "PLSDA"
    )    
    ggplot2::ggsave(paste0("C://Users//amshu//OneDrive - Emory//Preston, Joshua's files - Amshu Josh PGD//Final Analysis//PLSDA_combinedTime_ComparePGD.png"), pca_time_combined$plot, width = 5, height = 5, units = "in", dpi = 600)
    
#+ 3.2: Create heatmaps
  #- 3.2.1: Prepare data for heatmap
    heatmap_data_12h <- UFT_12h %>%
    select(Sample_ID, Clinical_PGD, Time, any_of(features_to_keep))
    heatmap_data_24h <- UFT_24h %>%
    select(Sample_ID, Clinical_PGD, Time, any_of(features_to_keep))
    heatmap_data_combTime <- UFT_12and24 %>%
    select(Sample_ID, Clinical_PGD, Time, any_of(features_to_keep))
  #- 3.2.2: Create heatmaps with different feature selections
    anova_1000 <- make_heatmap(
      heatmap_data_combTime,
      feature_selector = "ttest",
      top_features = 250,
      filename = "CombinedTime_TimeTtest"
      )
    anova_250 <- make_heatmap(
      heatmap_data_combTime,
      feature_selector = "anova",
      top_features = 250,
      filename = "CombinedTime_ANOVA"
    )
    hmap_pgd_12 <- make_heatmap(
      heatmap_data_12h,
      feature_selector = "ttest",
      top_features = 250,
      filename = "12h_ttest"
      )
    hmap_pgd_24 <- make_heatmap(
      heatmap_data_24h,
      feature_selector = "ttest",
      top_features = 250,
      filename = "24h_ttest"
      )
    variance_combined_1000 <- make_heatmap(
      heatmap_data_combTime,
      feature_selector = "variance",
      top_features = FALSE,
      filename = "comb_variance"
    )
  #- 3.2.3: Create cluster data for PCA analysis
    UFT_with_clusters <- combined_UFT %>%
      left_join(variance_combined_1000$cluster_df, by = "Sample_ID") %>%
      mutate(Cluster = factor(paste0("Cluster ", Cluster), levels = c("Cluster 1", "Cluster 2"))) %>%
      select(Patient_ID, Cluster, any_of(features_to_keep)) %>%
      arrange(Cluster)
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
#+ 3.5: Final Figure Creations
#- 3.5.1: Create heatmap for Figure 1
  heatmap_1A <- patchwork::wrap_elements(variance_1000$heatmap_plot$gtable)
#- 3.5.2 PCAs - Safe execution (uses unified theme)
  tryCatch(
    {
      pca_1D <- make_PCA(
        data = variant_data,
        ellipse_colors = variant_colors,
        point_size = 1
      )$plot + theme_pub_pca()

      pca_1F <- make_PCA(
        data = T_stage_data,
        ellipse_colors = T_stage_bin_colors,
        point_size = 1
      )$plot + theme_pub_pca()

      pca_1E <- make_PCA(
        data = Sex_data,
        ellipse_colors = sex_colors,
        point_size = 1
      )$plot + theme_pub_pca()

      pca_1C <- make_PCA(
        data = UFT_with_clusters,
        ellipse_colors = cluster_colors,
        point_size = 1
      )$plot + theme_pub_pca()

      cat("✅ All PCA plots created successfully!\n")
    },
    error = function(e) {
      cat("❌ Error creating PCA plots:", e$message, "\n")
      cat("Check your data and color definitions.\n")
    }
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
