balloon_plot <- function(balloon_dat,
                         name_col = "Identified_Name",
                         value_col = "Value", # assumed to be log2FC
                         label_col = "Label") {
  stopifnot(all(c(name_col, value_col, label_col) %in% names(balloon_dat)))

  bd <- balloon_dat %>%
    dplyr::filter(.data[[label_col]] == TRUE) %>%
    dplyr::select(dplyr::all_of(c(name_col, value_col))) %>%
    dplyr::mutate(
      log2FC_clamp = pmax(pmin(.data[[value_col]], 2.5), -2.5),
      magn         = abs(log2FC_clamp),
      direction    = ifelse(log2FC_clamp > 0, "Up in PGD", "Down in PGD")
    ) %>%
    dplyr::arrange(log2FC_clamp) %>%
    dplyr::mutate(
      !!name_col := factor(.data[[name_col]], levels = .data[[name_col]])
    )

  p <- ggplot2::ggplot(bd, ggplot2::aes(x = 1, y = .data[[name_col]])) +
    ggplot2::geom_point(
      ggplot2::aes(size = magn, color = direction),
      alpha = 0.95,
      show.legend = TRUE
    ) +
    ggplot2::scale_size_continuous(
      name   = expression("|log"[2] * "FC|"),
      range  = c(3, 10),
      limits = c(0, 2.5),
      breaks = c(0.5, 1, 1.5, 2, 2.5),
      labels = c("0.5", "1.0", "1.5", "2.0", "2.5")
    ) +
    ggplot2::scale_color_manual(
      values = c("Up in PGD" = "#800017", "Down in PGD" = "#113d6a"),
      guide  = "none" # remove color legend
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        override.aes = list(color = "black", shape = 16) # black dots in legend
      )
    ) +
    ggplot2::theme_minimal(base_family = "Arial") +
    ggplot2::theme(
      # remove all gridlines
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # axis labels off
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),

      # keep Y axis text for metabolite names
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 12, face = "bold", color = "black"),

      # legend formatting
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 11, color = "black", family = "Arial"),
      legend.title = ggplot2::element_text(size = 12, face = "bold", color = "black", family = "Arial")
    )

  list(
    balloon_dat = bd,
    ball_plot   = p
  )
}
