make_volcano <- function(FC_list, ttest_res, key) {
  # ---- 1) Build FC table robustly ----
  # Accept either a data.frame/tibble with Metabolite/Time/Value
  # or a list/data.frame with rownames as Metabolite and cols "Value","Time".
  if (all(c("Metabolite", "Time", "Value") %in% names(FC_list))) {
    FC_table <- FC_list %>%
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
  ) %>%
    # Attach human-readable names (compatible with combined_key_clean)
    dplyr::left_join(key, by = "Metabolite")

  # Expect a column named p_value in ttest_res
  if (!"p_value" %in% names(volcano_data)) {
    stop("`ttest_res` must contain a column named `p_value`.")
  }

  # ---- 3) Color classification ----
  thr <- log2(1.5)

  volcano_data <- volcano_data %>%
    dplyr::mutate(
      Legend = dplyr::case_when(
        p_value < 0.05 & Value >= thr ~ "Up in PGD",
        p_value < 0.05 & Value <= -thr ~ "Down in PGD",
        TRUE ~ "Not Significant"
      )
    )

  # ---- 4) Pick top/bottom 10 among the significant for labels ----
  top10sig <- volcano_data %>%
    dplyr::filter(p_value < 0.05) %>%
    dplyr::arrange(Value) %>%
    dplyr::slice_head(n = 10)

  bottom10sig <- volcano_data %>%
    dplyr::filter(p_value < 0.05) %>%
    dplyr::arrange(dplyr::desc(Value)) %>%
    dplyr::slice_head(n = 10)

  label_set <- dplyr::bind_rows(top10sig, bottom10sig) %>%
    dplyr::distinct(Metabolite) %>%
    dplyr::pull(Metabolite)

  volcano_data <- volcano_data %>%
    dplyr::mutate(Label = Metabolite %in% label_set)

  # ---- 5) Plot ----
  volcano_plot <- ggplot2::ggplot(
    volcano_data,
    ggplot2::aes(x = Value, y = -log10(p_value), color = Legend)
  ) +
    ggplot2::geom_point(size = 2, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c(
        "Not Significant" = "gray70",
        "Up in PGD"       = "#800017",
        "Down in PGD"     = "#113d6a"
      ),
      name = NULL
    ) +
    ggplot2::theme_light(base_family = "Arial") +
    ggplot2::labs(
      x = expression(log[2]("Fold Change")),
      y = expression(-log[10](p))
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

      legend.position = "none", # ✅ hides all legends

      # # Legend
      # legend.position = c(0.05, 0.95), # top-left inside plot
      # legend.justification = c("left", "top"),
      # legend.background = ggplot2::element_rect(fill = alpha("white", 0.7), color = NA),
      # legend.key = ggplot2::element_blank(),
      # legend.title = ggplot2::element_blank(),
      # legend.text = ggplot2::element_text(size = 10, face = "bold", color = "black"),

      # General
      strip.text = ggplot2::element_text(size = 12, face = "bold", color = "black"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.2),
      axis.ticks = ggplot2::element_blank()
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(shape = 16, size = 3) # legend dots smaller
    ))

  list(
    volcano_data = volcano_data,
    volcano_plot = volcano_plot
  )
}
