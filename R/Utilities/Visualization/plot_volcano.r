make_volcano <- function(FC_list, ttest_res) {
  # ---- 1) Build FC table robustly ----
  # Accept either a data.frame/tibble with Metabolite/Time/Value
  # or a list/data.frame with rownames as Metabolite and cols "Value","Time".
  if (all(c("Metabolite", "Time", "Value") %in% names(FC_list))) {
    FC_table <- FC_list |>
      dplyr::select(Metabolite, Time, Value)
  } else {
    # Fallback to old structure
    FC_table <- data.frame(
      Metabolite = rownames(FC_list),
      Value = FC_list$Value,
      Time = FC_list$Time,
      row.names = NULL,
      check.names = FALSE
    )
  }

  # Ensure character Time for a clean merge
  FC_table$Time <- as.character(FC_table$Time)
  ttest_res$Time <- as.character(ttest_res$Time)

  # ---- 2) Merge FC with test results ----
  volcano_data <- dplyr::full_join(
    FC_table,
    ttest_res,
    by = c("Metabolite", "Time")
  )

  # Expect a column named p_value in ttest_res
  if (!"p_value" %in% names(volcano_data)) {
    stop("`ttest_res` must contain a column named `p_value`.")
  }

  # ---- 3) Color classification ----
  thr <- log2(1.5)

  volcano_data <- volcano_data |>
    dplyr::mutate(
      Legend = dplyr::case_when(
        p_value < 0.05 & Value >= thr ~ "Up in sPGD",
        p_value < 0.05 & Value <= -thr ~ "Down in sPGD",
        TRUE ~ "Not Significant"
      )
    )

  # ---- 4) Pick top/bottom 10 among the significant for labels ----
  top10sig <- volcano_data |>
    dplyr::filter(p_value < 0.05) |>
    dplyr::arrange(Value) |>
    dplyr::slice_head(n = 10)

  bottom10sig <- volcano_data |>
    dplyr::filter(p_value < 0.05) |>
    dplyr::arrange(dplyr::desc(Value)) |>
    dplyr::slice_head(n = 10)

  label_set <- dplyr::bind_rows(top10sig, bottom10sig) |>
    dplyr::distinct(Metabolite) |>
    dplyr::pull(Metabolite)

  volcano_data <- volcano_data |>
    dplyr::mutate(Label = Metabolite %in% label_set)

  # ---- 5) Plot ----
  volcano_plot <- ggplot2::ggplot(
    volcano_data,
    ggplot2::aes(x = Value, y = -log10(p_value), color = Legend)
  ) +
    ggplot2::geom_point(size = 1.5, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c(
        "Not Significant" = "gray70",
        "Up in sPGD"       = "#800017",
        "Down in sPGD"     = "#113d6a"
      ),
      breaks = c("Down in sPGD", "Up in sPGD"),
      name = NULL
    ) +
    ggplot2::theme_light(base_family = "Arial") +
    ggplot2::labs(
      x = expression(bold(log[2]("Fold Change"))),
      y = expression(bold(-log[10](p)))
    ) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-thr, thr), linetype = "dashed", color = "black") +
    ggplot2::scale_x_continuous(
      limits = c(-3, 3),
      breaks = seq(-3, 3, 1),
      minor_breaks = seq(-3, 3, 1)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-0.25, 4),
      breaks = seq(0, 4, 1),
      minor_breaks = seq(-0.5, 4, 0.5)
    ) +
    ggplot2::theme(
      # Axis titles and labels
      axis.title.x = ggplot2::element_text(size = 15, face = "bold", color = "black"),
      axis.title.y = ggplot2::element_text(size = 15, face = "bold", color = "black"),
      axis.text.x = ggplot2::element_text(size = 12, face = "bold", color = "black"),
      axis.text.y = ggplot2::element_text(size = 12, face = "bold", color = "black"),

      # Legend styling - centered at top, universal spacing
      legend.position = "top",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.margin = ggplot2::margin(t = 0, r = 0, b = -3, l = 0),
      legend.key.width = grid::unit(0.35, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.key.size = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(0.1, "cm"),
      legend.text = ggplot2::element_text(size = 8, face = "bold", color = "black", margin = ggplot2::margin(l = 4, r = 4)),
      legend.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),

      # General
      strip.text = ggplot2::element_text(size = 12, face = "bold", color = "black"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.2),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = grid::unit(0.15, "cm")
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(shape = 16, size = 3) # legend dots smaller
    ))

  list(
    volcano_data = volcano_data,
    volcano_plot = volcano_plot
  )
}
