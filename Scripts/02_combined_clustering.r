#* 3: Combined Clustering analysis
#+ 3.1: Run Exploratory PLSDA
  #- 3.1.1: Define datasets
    UFT_12h <- UFT_nomiss %>% 
      filter(Time == 12) %>%
      select(Sample_ID, Time, Clinical_PGD, all_of(untargeted_features))
    UFT_24h <- UFT_nomiss %>%
      filter(Time == 24) %>%
      select(Sample_ID, Time, Clinical_PGD, all_of(untargeted_features))
    UFT_12and24 <- UFT_nomiss %>%
      select(Sample_ID, Time, Clinical_PGD, all_of(untargeted_features))
  #- 3.1.2: Define colors
    cluster_colors <- c(
      "Y" = "#94001E",
      "N" = "#03507D"
    )
  #- 3.1.3: Set Plot Specs
    plot_specs <- tribble(
      ~data, ~comp_x, ~comp_y, ~method, ~outpath,
      UFT_12h, 1, 2, "PLSDA", "Outputs/PLSDA/PLSDA_12h_ComparePGD.png",
      UFT_24h, 1, 2, "PLSDA", "Outputs/PLSDA/PLSDA_24h_ComparePGD.png",
      UFT_12and24, 1, 2, "PLSDA", "Outputs/PLSDA/PLSDA_combinedTime_ComparePGD.png",
      UFT_12h, 1, 2, "PCA", "Outputs/PCA/PCA_12h_ComparePGD.png",
      UFT_24h, 1, 2, "PCA", "Outputs/PCA/PCA_24h_ComparePGD.png",
      UFT_12and24, 1, 2, "PCA", "Outputs/PCA/PCA_combinedTime_ComparePGD.png"
    )
  #- Run all the PLSDAs and PCAs in one sweep
    walk2(
      .x = pmap(plot_specs, ~ make_PCA(..1,
        comp_x = ..2, comp_y = ..3,
        group_var = "Clinical_PGD", method = ..4
      )),
      .y = plot_specs$outpath,
      ~ ggsave(.y, .x$plot, width = 5, height = 5, units = "in", dpi = 600)
    )
#+ 3.2: Create heatmaps
  #- 3.2.1: Prepare data for heatmap
    heatmap_data_12h <- UFT_12h %>%
    select(Sample_ID, Clinical_PGD, Time, any_of(untargeted_features))
    heatmap_data_24h <- UFT_24h %>%
    select(Sample_ID, Clinical_PGD, Time, any_of(untargeted_features))
    heatmap_data_combTime <- UFT_12and24 %>%
    select(Sample_ID, Clinical_PGD, Time, any_of(untargeted_features))
  #- 3.2.2: Create heatmaps with different feature selections
#* HCA with heatmaps, PLS-DA, and volcano plots
# Divide data into groups based on clinical PGD status and time
unt_24 <- UFT_nomiss %>%
  filter(Time == 24) %>%
  select(-Patient, -Time, -Sex, -Age)
unt_12 <- UFT_nomiss %>%
  filter(Time == 12) %>%
  select(-Patient, -Time, -Sex, -Age)
# Metabolite values averaged between 12h and 24h
metcols <- names(unt_24)[1:2]
geomean <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  exp(mean(log(x + 1))) - 1
}
unt_avg <- bind_rows(unt_12, unt_24) %>%
  group_by(Clinical_PGD) %>%
  summarise(across(-Sample_ID, geomean, na.rm = TRUE), .groups = "drop")
#+ Heatmap of metabolite clusters
# Function to create heatmap with hierarchical clustering
create_heatmap <- function(data, title) {
  # Extract metabolite data and set row names
  df_sub <- data

  # Extract metabolite matrix only
  met_cols <- c("Sample_ID", "Clinical_PGD")
  data_cols <- setdiff(colnames(df_sub), met_cols)
  met_data <- as.matrix(df_sub[, data_cols])
  is.numeric(met_data)
  met_scaled <- as.matrix(scale(met_data))
  row_order <- order(df_sub$Clinical_PGD)
  rownames(met_scaled) <- df_sub$Sample_ID
  met_scaled <- met_scaled[, colSums(is.na(met_scaled)) == 0]
  met_scaled <- met_scaled[row_order, ]
  annotation <- data.frame(Clinical_PGD = df_sub$Clinical_PGD)
  rownames(annotation) <- df_sub$Sample_ID
  annotation <- annotation[row_order, , drop = FALSE]

  # Define colors
  ann_colors <- list(
    Clinical_PGD = c("Y" = "#800017", "N" = "#113d6a")
  )
  pheatmap(
    met_scaled,
    cluster_rows = FALSE, # hierarchical clustering of samples
    cluster_cols = TRUE, # hierarchical clustering of metabolites
    annotation_row = annotation,
    annotation_colors = ann_colors,
    scale = "row", # optional: scale metabolites to mean=0, sd=1
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 8,
    main = title,
    filename = paste0("./Outputs/Heatmaps/Heatmap_", title, ".png"),
  )
}
# Create heatmaps for each group
heatmap_24h <- create_heatmap(unt_24, "PGD vs Non-PGD at 24h")








    anova_1000 <- make_heatmap(
      UFT_nomiss,
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
      UFT_12h,
      feature_selector = "ttest",
      top_features = 250
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
