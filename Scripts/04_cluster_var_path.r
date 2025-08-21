#* 4: Pathology analysis by cluster
#+ 4.1 Prepare Data
  cluster_path_analysis <- variance_study %>%
    select(-any_of(feature_cols)) %>%
    left_join(tumor_pathology %>% select(Patient_ID, T_computed), by = "Patient_ID") %>%
    select(T_computed, Cluster, LVI, MFC, T_stage, Variant)
#+ 4.2: Faceted Stacked Bar Chart
  #- 4.2.1 Prepare Data for Faceted Stacked Bar Plot
    stacked_bar_data_facet <- cluster_path_analysis %>%
      count(Cluster, Variant, T_computed, name = "count") %>%
      mutate(
        Cluster = factor(Cluster, levels = c("Cluster 1", "Cluster 2")),
        Variant = factor(Variant, levels = c("PTC", "FV-PTC", "FTC")),
        T_computed = factor(T_computed, levels = c("T1", "T2", "T3", "T4"))
      )
  #- 4.2.2 Create Faceted Stacked Bar Chart
    facet_3A <- ggplot(stacked_bar_data_facet, aes(x = Variant, y = count, fill = T_computed)) +
      geom_col(width = 0.7, color = "black", linewidth = 0.4, position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = T_stage_cluster_colors, name = "T Stage") +
      guides(fill = guide_legend(reverse = TRUE)) +
      scale_y_continuous(limits = c(0, 21), breaks = seq(0, 20, by = 5), expand = c(0, 0)) +
      facet_grid(~Cluster, space = "free_x", scales = "free_x") +
      labs(
        x = NULL,
        y = "Patients (n)"
      ) +
      theme_publication_bars() +
      theme(
        plot.margin = grid::unit(c(8, 5, 5, 5), "pt"),
        axis.line = ggplot2::element_line(linewidth = 0.5, lineend = "square"),
        axis.ticks = ggplot2::element_line(linewidth = 0.3),
        axis.ticks.length = grid::unit(0.2, "cm"),
        axis.title.y = ggplot2::element_text(size = 9, face = "bold", family = "Arial"),
        axis.title.x = ggplot2::element_text(size = 9, face = "bold", family = "Arial"),
        axis.text.y  = ggplot2::element_text(size = 8, face = "bold", family = "Arial"),
        axis.text.x  = ggplot2::element_text(size = 8, face = "bold", family = "Arial", angle = 45, hjust = 1, vjust = 0.97),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.ontop = FALSE,
        strip.text = ggplot2::element_text(size = 13, face = "bold", family = "Arial"),
        legend.title = ggplot2::element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.key.size = grid::unit(0.4, "cm"),
        legend.key.width = grid::unit(0.4, "cm"),
        legend.key.height = grid::unit(0.4, "cm"),
        legend.text = ggplot2::element_text(size = 8, face = "bold", family = "Arial"),
        plot.title = ggplot2::element_text(size = 12, face = "bold", family = "Arial", hjust = 0.5)
      )
