#* 5: Targeted  
#+ 5.1: Merge TFT with path data
  #- 5.1.2: Join with TFT_metaboanalyst_log2
    TFT_metaboanalyst_log2_path_analysis <- TFT_metaboanalyst_log2_path %>%
      left_join(variance_1000$cluster_df, by = "Patient_ID") %>%
      select(Patient_ID, T_computed, Variant, Cluster, Variant, everything(), -c(Age, Sex, MFC, ETE, LVI, LD, T_computed_bin))
#+ 5.2: Statistical comparisons using targeted_metabolite_comparison function
  #- 5.2.1: ANOVA comparison between T_computed stages
    t_computed_results <- targeted_metabolite_comparison(
      data = TFT_metaboanalyst_log2_path_analysis,
      grouping_var = "T_computed",
      test_type = "anova",
      exclude_cols = c("Patient_ID", "Cluster", "Variant"),
      exclude_isotopes = TRUE
    )
  #- 5.2.3: ANOVA comparison between Variants
    variant_results <- targeted_metabolite_comparison(
      data = TFT_metaboanalyst_log2_path_analysis,
      grouping_var = "Variant",
      test_type = "anova",
      exclude_cols = c("Patient_ID", "T_computed", "Cluster"),
      exclude_isotopes = TRUE
    )
  #- 5.2.3: ANOVA comparison between Variants
    cluster_results <- targeted_metabolite_comparison(
      data = TFT_metaboanalyst_log2_path_analysis,
      grouping_var = "Cluster",
      test_type = "anova",
      exclude_cols = c("Patient_ID", "T_computed", "Variant"),
      exclude_isotopes = TRUE
    )
#+ 5.3: Filter metabolites and create top 10 dataset
  #- 5.3.1: Get top 10 metabolites from variant comparison (by FDR p-value)
    top_10_variant_metabolites <- variant_results %>%
      filter(significant_fdr == TRUE) %>%
      arrange(p_value_fdr) %>%
      slice_head(n = 10) %>%
      pull(Metabolite)
  #- 5.3.2: Create top 5 dataset with metadata for variants
    TFT_top5_dataset_variant <- TFT_metaboanalyst_log2_path_analysis %>%
      select(Variant, any_of(top_10_variant_metabolites)) %>%
      select(-c("Cinnamoylglycine_HILIC_206.0811588","L-Citrulline_HILIC_176.1024546", "Pyroglutamic acid_C18_128.0353782")) %>%
      rename_with(~ sub("_.*", "", .x)) %>%
      rename("Citrulline" = "L-Citrulline") %>%
      rename("Indolelactic Acid*" = "indolelactic acid") %>%
      rename("Oxoproline*" = "Oxoproline") %>%
      select(1:6)
  #- 5.3.3: Create top 10 metabolites from cluster comparison (by FDR p-value)
    top_10_cluster_metabolites <- cluster_results %>%
      filter(significant_fdr == TRUE) %>%
      arrange(p_value_fdr) %>%
      slice_head(n = 10) %>%
      pull(Metabolite)
  #- 5.3.4: Create top 5 dataset with metadata for clusters
    TFT_top5_dataset_cluster <- TFT_metaboanalyst_log2_path_analysis %>%
      select(Cluster, any_of(top_10_cluster_metabolites)) %>%
      select(-c("Niacinamide_HILIC_123.0552695")) %>%
      rename_with(~ sub("_.*", "", .x)) %>%
      rename("Nicotinamide*" = "Nicotinamide") %>%
      select(1:6)
#+ 5.4: Data prep for visualization
  #- 5.4.1: Prepare data for plotting
    #_Get metabolite columns maintaining current order
    metab_cols <- setdiff(names(TFT_top5_dataset_variant), "Variant")
    #_Convert to long format
    df_long <- TFT_top5_dataset_variant %>%
      mutate(Variant = factor(Variant, levels = names(variant_colors))) %>%
      pivot_longer(
        cols = all_of(metab_cols),
        names_to = "Metabolite",
        values_to = "Value"
      ) %>%
      mutate(Metabolite = factor(Metabolite, levels = metab_cols)) %>%
      mutate(Metabolite = forcats::fct_recode(
        Metabolite,
        "Indolelactic\nAcid*" = "Indolelactic Acid*"
      ))
    #_Compute means per metabolite × variant for bar heights
    sum_df <- df_long %>%
      group_by(Metabolite, Variant) %>%
      summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")
  #- 5.4.4: Cluster distribution bar plot
    #_Prepare cluster dataset for plotting (like 5.4.1)
    metab_cols_cluster <- setdiff(names(TFT_top5_dataset_cluster), "Cluster")
    df_long_cluster <- TFT_top5_dataset_cluster %>%
      mutate(Cluster = factor(paste0("Cluster ", Cluster), levels = names(cluster_colors))) %>%
      pivot_longer(
        cols = all_of(metab_cols_cluster),
        names_to = "Metabolite",
        values_to = "Value"
      ) %>%
      mutate(Metabolite = factor(Metabolite, levels = metab_cols_cluster)) %>%
      mutate(Metabolite = forcats::fct_recode(
        Metabolite,
        "Pantothenic\nAcid" = "Pantothenic acid"
      ))
    #_Compute means per metabolite × cluster for bar heights
    sum_df_cluster <- df_long_cluster %>%
      group_by(Metabolite, Cluster) %>%
      summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")      
#+ 5.5: Create dot plot with floating mean bars
#- 5.5.1: Clusters based targeted top 5
  clusters_targ_3B <- ggplot() +
    # mean bars (no SD)
    geom_col(
      data = sum_df_cluster,
      aes(
        x = Metabolite, y = mean_value,
        fill = Cluster, color = Cluster, group = Cluster
      ),
      position = position_dodge2(width = 0.8, preserve = "single"),
      width = 0.72,
      linewidth = 0.5,
      na.rm = TRUE
    ) +
    # mean dots on bar centers
    geom_point(
      data = sum_df_cluster,
      aes(x = Metabolite, y = mean_value, color = Cluster, group = Cluster),
      position = position_dodge(width = 0.8),
      size = 1.2
    ) +
    # jittered raw points
    geom_jitter(
      data = df_long_cluster,
      aes(x = Metabolite, y = Value, color = Cluster),
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      size = 0.5, alpha = 0.8, show.legend = FALSE
    ) +
    # palettes (use lighter fills, solid outlines)
    scale_fill_manual(values = cluster_light, drop = FALSE) +
    scale_color_manual(values = cluster_colors, drop = FALSE) +
    # keep bars alive with 0 baseline; zoom view for readability
    scale_y_pub(lims = c(0, 35), expand = c(0, 0)) +
    coord_cartesian(ylim = c(8, 35), clip = "on") +
    theme_pub_dotbar() +
    theme_pub_barplot() +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      fill  = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    labs(
      x = NULL,
      y = expression(bold(log[2] ~ "(Spectral Intensity)")),
      fill = "Cluster", color = "Cluster"
    )
#- 5.5.2: Variant based targeted top 5
  variant_targ_3C <- ggplot() +
    # === MEAN BARS (no SD) ===
    geom_col(
      data = sum_df,
      aes(
        x = Metabolite, y = mean_value,
        fill = Variant, color = Variant, group = Variant
      ),
      position = position_dodge2(width = 0.8, preserve = "single"),
      width = 0.72,
      linewidth = 0.5, # solid outline
      na.rm = TRUE
    ) +
    # === MEAN DOTS (on bar centers) ===
    geom_point(
      data = sum_df,
      aes(x = Metabolite, y = mean_value, color = Variant, group = Variant),
      position = position_dodge(width = 0.8),
      size = 1.2
    ) +
    # === JITTERED RAW POINTS ===
    geom_jitter(
      data = df_long,
      aes(x = Metabolite, y = Value, color = Variant),
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      size = 0.5, alpha = 0.8, show.legend = FALSE
    ) +
    # === COLOR SCALES (lighter fills, dark outlines) ===
    scale_fill_manual(values = variant_light, drop = FALSE) +
    scale_color_manual(values = variant_colors, drop = FALSE) +
    # === AXES & VIEW ===
    scale_y_pub(lims = c(0, 35), expand = c(0, 0)) +
    coord_cartesian(ylim = c(8, 35), clip = "on") +
    # === THEMES ===
    theme_pub_dotbar() +
    theme_pub_barplot() +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      fill  = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    labs(
      x = NULL,
      y = expression(bold(log[2] ~ "(Spectral Intensity)")),
      fill = "Variant", color = "Variant"
    )
