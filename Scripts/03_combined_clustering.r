#* 3: Combined Clustering analysis
#+ 3.1: Run Exploratory PCAs
  #- 3.1.1: Define datasets
    variant_data <- UFT_metaboanalyst_log2_path %>% 
      select(Patient_ID, Variant, all_of(feature_cols))
    T_stage_data <- UFT_metaboanalyst_log2_path %>%
      select(Patient_ID, T_computed_bin, all_of(feature_cols))
    LVI_data <- UFT_metaboanalyst_log2_path %>%
      select(Patient_ID, LVI, all_of(feature_cols))
    Sex_data <- UFT_metaboanalyst_log2_path %>%
    select(Patient_ID, Sex, all_of(feature_cols))
  #- 3.1.2: Create variant PCAs
    variant_pca_12 <- make_PCA(
      data = variant_data,
      ellipse_colors = variant_colors,
      comp_x = 1, comp_y = 2
    )
    variant_pca_23 <- make_PCA(
      data = variant_data,
      ellipse_colors = variant_colors,
      comp_x = 2, comp_y = 3
    )
    variant_pca_34 <- make_PCA(
      data = variant_data,
      ellipse_colors = variant_colors,
      comp_x = 3, comp_y = 4
    )
  #- 3.1.3: Create T-stage PCAs
    T_stage_PCA_12 <- make_PCA(
      data = T_stage_data,
      ellipse_colors = T_stage_colors,
      comp_x = 1, comp_y = 2
    )
    T_stage_PCA_23 <- make_PCA(
      data = T_stage_data,
      ellipse_colors = T_stage_colors,
      comp_x = 2, comp_y = 3
    )
    T_stage_PCA_34 <- make_PCA(
      data = T_stage_data,
      ellipse_colors = T_stage_colors,
      comp_x = 3, comp_y = 4
    )
  #- 3.1.4: Create LVI PCAs
    LVI_PCA_12 <- make_PCA(
      data = LVI_data,
      ellipse_colors = LVI_colors,
      comp_x = 1, comp_y = 2
    )
    LVI_PCA_23 <- make_PCA(
      data = LVI_data,
      ellipse_colors = LVI_colors,
      comp_x = 2, comp_y = 3
    )
    LVI_PCA_34 <- make_PCA(
      data = LVI_data,
      ellipse_colors = LVI_colors,
      comp_x = 3, comp_y = 4
    )
  #- 3.1.5: Create Sex PCAs
    Sex_PCA_12 <- make_PCA(
      data = Sex_data,
      ellipse_colors = sex_colors,
      comp_x = 1, comp_y = 2
    )
    Sex_PCA_23 <- make_PCA(
      data = Sex_data,
      ellipse_colors = sex_colors,
      comp_x = 2, comp_y = 3
    )
    Sex_PCA_34 <- make_PCA(
      data = Sex_data,
      ellipse_colors = sex_colors,
      comp_x = 3, comp_y = 4
    )
#+ 3.2: Create heatmaps
  #- 3.2.1: Prepare data for heatmap
    heatmap_data_T <- UFT_metaboanalyst_log2_path %>%
    select(Patient_ID, Variant, T_computed_bin, any_of(feature_cols)) %>%
    rename("T_stage" = T_computed_bin)
  #- 3.2.2: Create heatmaps with different feature selections
    variance_1000 <- make_heatmap(
      heatmap_data_T,
      variant_colors = variant_colors,
      feature_selector = "variance",
      top_features = 1000,
      annotate_t_stage = TRUE,
      T_stage_colors = T_stage_bin_colors,
      cluster_colors = cluster_colors
    )
    anova_1000 <- make_heatmap(
      heatmap_data_T,
      variant_colors = variant_colors,
      feature_selector = "anova",
      top_features = 1000,
      annotate_t_stage = TRUE,
      T_stage_colors = T_stage_colors_heatmap,
      cluster_colors = cluster_colors
    )
    variance_full <- make_heatmap(
      heatmap_data_T,
      variant_colors = variant_colors,
      feature_selector = "variance",
      top_features = FALSE,
      annotate_t_stage = TRUE,
      T_stage_colors = T_stage_bin_colors,
      cluster_colors = cluster_colors
    )
  #- 3.2.3: Create cluster data for PCA analysis
    UFT_with_clusters <- UFT_metaboanalyst_log2_path %>%
      left_join(variance_1000$cluster_df, by = "Patient_ID") %>%
      mutate(Cluster = factor(paste0("Cluster ", Cluster), levels = c("Cluster 1", "Cluster 2"))) %>%
      select(Patient_ID, Cluster, any_of(feature_cols)) %>%
      arrange(Cluster)
#+ 3.4: Run PERMANOVA analysis
  #- 3.4.1: Define feature columns
    permanova_features <- rownames(variance_1000$M)
  #- 3.4.1: Extract features and prepare data
    variance_study <- UFT_metaboanalyst_log2_path %>%
      left_join(variance_1000$cluster_df, by = "Patient_ID") %>%
      select(Patient_ID, Cluster, Variant, T_computed_bin, Sex, Age, MFC, LVI, any_of(permanova_features)) %>%
      rename(T_stage = T_computed_bin) %>%
      mutate(
        across(c(Sex, LVI, Variant, MFC), as.factor),
        Age = as.numeric(Age),
        T_stage = factor(T_stage, levels = c("T1-T2", "T3-T4")),
        Cluster = factor(if_else(Cluster == 1, "Cluster 1", "Cluster 2"), 
                        levels = c("Cluster 1", "Cluster 2"))
      ) %>%
      arrange(Cluster)
  #- 3.4.2: Define PERMANOVA variables
    permanova_variables <-  c("T_stage", "Sex", "LVI", "Variant", "MFC", "Age", "Cluster")
  #- 3.4.3: Impute missing values
    #! IMPUTED for LVI variable, one sample missing data
    meta_i <- variance_study %>% 
      select(all_of(permanova_variables))
    meta <- complete(mice(
      meta_i,
      m = 1,
      method = replace(setNames(rep("", ncol(meta_i)), colnames(meta_i)), "LVI", "logreg"),
      predictorMatrix = {
        p <- matrix(0, ncol(meta_i), ncol(meta_i),
          dimnames = list(colnames(meta_i), colnames(meta_i))
        )
        p["LVI", c("T_stage", "Sex", "Variant", "MFC", "Age")] <- 1
        p
      },
      seed = 123,
      printFlag = FALSE
    ), 1)
  #- 3.4.4: Extract feature data for analysis
    features_1000 <- variance_study %>%
      select(any_of(permanova_features))
  #- 3.4.5: Prepare metadata for individual tests
    meta_use <- variance_study %>%
      select(all_of(permanova_variables)) %>%
      mutate(
        across(c(T_stage, Sex, LVI, Variant, MFC, Cluster), as.factor),
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
