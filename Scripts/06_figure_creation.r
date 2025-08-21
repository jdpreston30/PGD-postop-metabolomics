#* 6: Figure Creation
  #+ 6.1: Figure 1 - Clustering Analysis Layout
    #- 6.1.1: Fix the border widths
      permanova_1B <- permanova_1B + theme_pub_simple(border_linewidth = 0.5) # Half the size of PCAs
      pca_1C <- pca_1C + theme_pub_pca(border_linewidth = 0.5) 
      pca_1D <- pca_1D + theme_pub_pca(border_linewidth = 0.5)
      pca_1E <- pca_1E + theme_pub_pca(border_linewidth = 0.5)
      pca_1F <- pca_1F + theme_pub_pca(border_linewidth = 0.5)
    #- 6.1.4: Assemble Figure 1 (Heatmap top 50%, vertical PERMANOVA + 4 PCAs bottom 50%)
      Figure_1 <- patchwork::wrap_plots(
        # Top 50%: Heatmap spanning full width
        heatmap_1A,
        # Bottom 50%: Vertical PERMANOVA (left) + 4 PCAs (2x2 grid, right)
        permanova_1B, pca_1C, pca_1D, pca_1E, pca_1F,
        design = "AAAA\nAAAA\nBBCD\nBBEF",  # Heatmap top 50%, vertical bar + 2x2 PCAs bottom 50%
        heights = c(1, 1, 0.5, 0.5)  # 50% heatmap, 50% bottom split
      ) +
        patchwork::plot_annotation(
          title = "Figure 1",
          tag_levels = "A",
          theme = ggplot2::theme(
            plot.title.position = "plot",
            plot.title = ggplot2::element_text(hjust = 0, face = "bold", family = "Arial", size = 16),
            plot.margin = grid::unit(c(0.3, 0.5, 0.3, 0.5), "in")
          )
        ) &
        ggplot2::theme(
          plot.tag.position = c(0, 0.98),
          plot.tag = ggplot2::element_text(size = 14, face = "bold", vjust = 0, hjust = 0, family = "Arial", color = "black")
        )
    #- 6.1.6: Save Figure 1 to PDF and PNG
      print_to_png(Figure_1, "Figure_1_preview")
  #+ 6.3: Figure 3 - Targeted Metabolomics
  Figure_3 <- patchwork::wrap_plots(
    facet_3A, clusters_targ_3B,                # Top 25%: A, B
    variant_targ_3C,
    design = "A\nB\nC",
    heights = c(0.25, 0.25, 0.25, 0.25)
  ) +
    patchwork::plot_annotation(
      title = "Figure 3",
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(hjust = 0, face = "bold", family = "Arial", size = 16),
        plot.margin = grid::unit(c(0.1, 1.9, 0.1, 1.9), "in")
      )
    ) &
    ggplot2::theme(
      plot.tag.position = c(0, 0.98),
      plot.tag = ggplot2::element_text(size = 14, face = "bold", vjust = 0, hjust = 0, family = "Arial", color = "black")
    )
  # Save Figure 3 to PNG (add PDF if needed)
  print_to_png(Figure_3, "Figure_3_preview")
   