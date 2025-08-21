#* Base theme (never hides legend)
  theme_pub_base <- function(base_size = 12, base_family = "Arial",
                             border_linewidth = 0.8,
                             major_frac = 0.5, # major grid = 50% of border
                             minor_frac = 0.5, # minor grid = 50% of major (i.e., 25% of border)
                             show_major_x = TRUE, show_major_y = FALSE,
                             show_minor_x = FALSE, show_minor_y = FALSE) {
    major_lwd <- border_linewidth * major_frac
    minor_lwd <- major_lwd * minor_frac

    ggprism::theme_prism(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        plot.margin = grid::unit(c(8, 8, 6, 6), "pt"),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_line(linewidth = 0.3),
        axis.ticks.length = grid::unit(0.2, "cm"),
        axis.title.x = ggplot2::element_text(size = 9, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 9, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 8, face = "bold"),

        # Border
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = border_linewidth),
        panel.background = ggplot2::element_blank(),

        # Major grids
        panel.grid.major.x = if (isTRUE(show_major_x)) {
          ggplot2::element_line(color = "gray80", linewidth = major_lwd)
        } else {
          ggplot2::element_blank()
        },
        panel.grid.major.y = if (isTRUE(show_major_y)) {
          ggplot2::element_line(color = "gray80", linewidth = major_lwd)
        } else {
          ggplot2::element_blank()
        },

        # Minor grids
        panel.grid.minor.x = if (isTRUE(show_minor_x)) {
          ggplot2::element_line(color = "gray90", linewidth = minor_lwd)
        } else {
          ggplot2::element_blank()
        },
        panel.grid.minor.y = if (isTRUE(show_minor_y)) {
          ggplot2::element_line(color = "gray90", linewidth = minor_lwd)
        } else {
          ggplot2::element_blank()
        },
        legend.title = ggplot2::element_text(size = 9, face = "bold"),
        legend.text = ggplot2::element_text(size = 8, face = "bold")
      )
  }
#* Simple (no legend)
  theme_pub_simple <- function(base_size = 12, base_family = "Arial",
                               border_linewidth = 0.8,
                               major_frac = 0.5, minor_frac = 0.5) {
    theme_pub_base(base_size, base_family,
      border_linewidth = border_linewidth,
      major_frac = major_frac, minor_frac = minor_frac,
      show_major_x = TRUE, show_major_y = FALSE,
      show_minor_x = FALSE, show_minor_y = FALSE
    )
  }
#* Dot Bar
  theme_pub_dotbar <- function(base_size = 12, base_family = "Arial",
                              border_linewidth = 0.8,
                              bar_border_scale = 1.5) {
    theme_pub_base(
      base_size, base_family,
      border_linewidth = border_linewidth * bar_border_scale,
      major_frac = 0.25, minor_frac = 0.25,
      show_major_x = FALSE, show_major_y = FALSE,
      show_minor_x = FALSE, show_minor_y = FALSE # ⟵ only horizontal minor grid
    ) +
      ggplot2::theme(
        panel.ontop       = FALSE, # gridlines behind bars/dots
        axis.text.x       = ggplot2::element_text(size = base_size * 0.75, angle = 0, vjust = 1),
        legend.position   = "top",
        legend.direction  = "horizontal",
        legend.box        = "horizontal",
        legend.margin     = ggplot2::margin(t = 0, r = 0, b = -4, l = 0),
        legend.key.width  = grid::unit(0.25, "cm"),
        legend.key.height = grid::unit(0.25, "cm"),
        legend.key.size   = grid::unit(0.25, "cm"),
        legend.spacing.x  = grid::unit(0.05, "cm"),
        legend.title      = ggplot2::element_blank()
      )
  }
#* Grouped (keeps legend; choose position) ----
  theme_pub_grouped <- function(base_size = 12, base_family = "Arial",
                                legend_position = "right",
                                legend_direction = "vertical") {
    theme_pub_base(base_size, base_family) +
      ggplot2::theme(
        legend.position = legend_position,
        legend.direction = legend_direction,
        # match bar plot frame explicitly
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.major.x = ggplot2::element_line(color = "gray80", linewidth = 0.3),
        panel.grid.major.y = ggplot2::element_blank()
      )
  }
#* Bar plot
  theme_pub_barplot <- function() {
    theme_pub_dotbar() +
      theme(
        panel.ontop        = FALSE,
        # kill vertical gridlines:
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        # keep only horizontal major grid:
        panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        axis.text.x        = element_text(size = 7, face = "bold"),
        legend.key.height  = grid::unit(0.28, "cm"),
        legend.key.width   = grid::unit(0.28, "cm"),
        legend.key.size    = grid::unit(0.28, "cm"),
        legend.spacing.x   = grid::unit(0.05, "cm"),
        legend.text        = element_text(size = 7, face = "bold")
      )
  }
#* PCA (square aspect + top horizontal legend)
  theme_pub_pca <- function(base_size = 12, base_family = "Arial",
                            border_linewidth = 0.8,
                            major_frac = 0.25, minor_frac = 0.25,
                            pca_border_scale = 2.0) {
    theme_pub_base(base_size, base_family,
      border_linewidth = border_linewidth * pca_border_scale,
      major_frac = major_frac, minor_frac = minor_frac,
      show_major_x = TRUE, show_major_y = TRUE,
      show_minor_x = FALSE, show_minor_y = FALSE
    ) +
      ggplot2::theme(
        aspect.ratio = 1,
        plot.margin = grid::unit(c(2, 8, 8, 8), "pt"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        legend.margin = ggplot2::margin(t = 0, r = 0, b = -3, l = 0),

        # 🔽 controls symbol size
        legend.key.width = grid::unit(0.35, "cm"),
        legend.key.height = grid::unit(0.35, "cm"),
        legend.key.size = grid::unit(0.35, "cm"),

        # 🔽 distance between symbol and text
        legend.spacing.x = grid::unit(0, "cm"),

        # 🔽 fine-tunes text alignment (0 = left, 1 = right; <0 moves text closer to key)
        legend.text.align = -10,
        legend.title = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(
          margin = ggplot2::margin(r = 0), hjust = 0.5
        )
      )
  }
#* Helpers
#+ Scale helpers (offset axes, tight expand)
  scale_x_pub <- function(lims = NULL, breaks = ggplot2::waiver(), minor_breaks = NULL,
                          expand = c(0, 0), offset_guides = TRUE) {
    ggplot2::scale_x_continuous(
      limits = lims, breaks = breaks, minor_breaks = minor_breaks, expand = expand,
      guide = if (isTRUE(offset_guides)) ggprism::guide_prism_offset() else "none"
    )
  }
  scale_y_pub <- function(lims = NULL, breaks = ggplot2::waiver(), minor_breaks = NULL,
                          expand = c(0, 0), offset_guides = TRUE) {
    ggplot2::scale_y_continuous(
      limits = lims, breaks = breaks, minor_breaks = minor_breaks, expand = expand,
      guide = if (isTRUE(offset_guides)) ggprism::guide_prism_offset() else "none"
    )
  }
#* Palettes and Colors
  #+ Color/Fill manual scales
    scale_color_variant <- function(...) ggplot2::scale_color_manual(values = variant_colors, ...)
    scale_fill_variant <- function(...) ggplot2::scale_fill_manual(values = variant_colors, ...)
  #+ Definitions
    variant_colors <- c(
      "PTC"    = "#DF8D0A",
      "FV-PTC" = "#23744E",
      "FTC"    = "#194992"
    )
    variant_light <- c(
      "PTC"    = "#e9bb71",
      "FV-PTC" = "#80ac95",
      "FTC"    = "#7a92bc"
    )
    LVI_colors <- c(
      "-LVI" = "#9C27B0", # purple
      "+LVI" = "#FFC107" # amber
    )
    cluster_colors <- c(
      "Cluster 1" = "#94001E",
      "Cluster 2" = "#03507D"
    )
    cluster_light <- c(
      "Cluster 1" = "#bb697b",
      "Cluster 2" = "#6f96b1"
    )
    sex_colors <- c(
      "Male"   = "#0a9af3",
      "Female" = "#d55e70"
    )
    T_stage_bin_colors <- c(
      "T1-T2" = "#dfba37",
      "T3-T4" = "#72061c"
    )
    T_stage_colors_heatmap <- c(
      "T1" = "#f2a3b3",
      "T2" = "#f2a3b3",
      "T3" = "#72061c",
      "T4" = "#72061c"
    )
    T_stage_cluster_colors <- c(
      "T1" = "#FFF096",
      "T2" = "#FDC586",
      "T3" = "#F45158",
      "T4" = "#72061c"
    )
  #+ List of palettes
    .palettes <- list(
      variant   = variant_colors,
      variant_light = variant_light,
      LVI       = LVI_colors,
      cluster   = cluster_colors,
      cluster_light = cluster_light,
      sex       = sex_colors,
      T_bin     = T_stage_bin_colors,
      T_heatmap = T_stage_colors_heatmap,
      T_cluster = T_stage_cluster_colors
    )
  #+ Get palette helper function
    get_palette <- function(name) {
      if (!name %in% names(.palettes)) stop(sprintf("Palette '%s' not found.", name))
      .palettes[[name]]
    }
  #+ Scale helpers (one-liners you can use in plots)
    scale_color_variant <- function(...) ggplot2::scale_color_manual(values = variant_colors, ...)
    scale_fill_variant <- function(...) ggplot2::scale_fill_manual(values = variant_colors, ...)
    scale_color_LVI <- function(...) ggplot2::scale_color_manual(values = LVI_colors, ...)
    scale_fill_LVI <- function(...) ggplot2::scale_fill_manual(values = LVI_colors, ...)
    scale_color_cluster <- function(...) ggplot2::scale_color_manual(values = cluster_colors, ...)
    scale_fill_cluster <- function(...) ggplot2::scale_fill_manual(values = cluster_colors, ...)
    scale_color_sex <- function(...) ggplot2::scale_color_manual(values = sex_colors, ...)
    scale_fill_sex <- function(...) ggplot2::scale_fill_manual(values = sex_colors, ...)
    scale_color_T_bin <- function(...) ggplot2::scale_color_manual(values = T_stage_bin_colors, ...)
    scale_fill_T_bin <- function(...) ggplot2::scale_fill_manual(values = T_stage_bin_colors, ...)
    scale_color_T_cluster <- function(...) ggplot2::scale_color_manual(values = T_stage_cluster_colors, ...)
    scale_fill_T_cluster <- function(...) ggplot2::scale_fill_manual(values = T_stage_cluster_colors, ...)
    scale_fill_T_heatmap <- function(...) ggplot2::scale_fill_manual(values = T_stage_colors_heatmap, ...)
  #+ Safety check: warn if data has levels with no color =====
    check_palette_levels <- function(x, palette, palette_name = deparse(substitute(palette))) {
      lv <- if (is.factor(x)) levels(x) else sort(unique(as.character(x)))
      missing <- setdiff(lv, names(palette))
      if (length(missing)) {
        warning(sprintf(
          "Palette '%s' missing %d level(s): %s",
          palette_name, length(missing), paste(missing, collapse = ", ")
        ), call. = FALSE)
      }
      invisible(TRUE)
    }