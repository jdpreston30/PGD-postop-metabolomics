#' Create Line Plot with Error Bars for Summary Data
#'
#' This function creates a line plot showing mean values over time with error bars (SEM or SD)
#' for two groups, with styling consistent with other plots.
#'
#' @param data Tibble with columns: time_var (x-axis), mean_var (y-axis), error_var (error bars), group_var (grouping)
#' @param time_var Character; name of the time variable column (default: "Time")
#' @param mean_var Character; name of the mean value column (default: "mean_succinate")
#' @param error_var Character; name of the error (SE or SD) column (default: "se_succinate")
#' @param group_var Character; name of the grouping variable column (default: "severe_PGD")
#' @param x_label Character; x-axis label (default: "Time (hours)")
#' @param y_label Character; y-axis label (default: "Succinate (mean ± SE)")
#' @param base_family Character; font family (default: "Arial")
#' @param text_scale Numeric; scaling factor for all text elements (default: 1.0)
#' @param line_size Numeric; line width (default: 1.2)
#' @param point_size Numeric; point size (default: 3)
#' @param error_width Numeric; error bar cap width (default: 0.1)
#' @param error_linewidth Numeric; error bar line width (default: 0.8)
#'
#' @return ggplot object
#'
#' @examples
#' \dontrun{
#'   p <- plot_line_summary(
#'     data = summary_data,
#'     time_var = "Time",
#'     mean_var = "mean_value",
#'     error_var = "se_value"
#'   )
#'   print(p)
#' }
#'
#' @export
plot_line_summary <- function(data,
                              time_var = "Time",
                              mean_var = "mean_succinate",
                              error_var = "se_succinate",
                              group_var = "severe_PGD",
                              x_label = "Time (hours)",
                              y_label = "Succinate (mean ± SE)",
                              base_family = "Arial",
                              text_scale = 1.0,
                              line_size = 1.2,
                              point_size = 3,
                              error_width = 0.1,
                              error_linewidth = 0.8) {
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  
  # Rename columns to standard names for plotting
  plot_data <- data %>%
    rename(
      time_plot = !!sym(time_var),
      mean_plot = !!sym(mean_var),
      error_plot = !!sym(error_var),
      group_plot = !!sym(group_var)
    )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = time_plot, y = mean_plot, group = group_plot, color = group_plot)) +
    geom_line(linewidth = line_size) +
    geom_point(size = point_size) +
    geom_errorbar(
      aes(ymin = mean_plot - error_plot, ymax = mean_plot + error_plot), 
      width = error_width, 
      linewidth = error_linewidth
    ) +
    scale_color_manual(
      values = c("N" = "#03507D", "Y" = "#94001E"),
      labels = c("N" = "No sPGD", "Y" = "sPGD"),
      name = NULL,
      guide = guide_legend(
        override.aes = list(
          linewidth = 1.2,
          linetype = 1,
          shape = 19,
          size = 2.5
        )
      )
    ) +
    scale_x_discrete(
      expand = expansion(add = c(0.3, 0.3))
    ) +
    labs(
      x = x_label,
      y = expression(bold("log")[2]*bold("(Peak Area)")),
      title = "Succinate"
    ) +
    theme_minimal(base_family = base_family) +
    theme(
      # Panel styling - matching diverging bars
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.545*2),
      
      # Axis styling - matching diverging bars
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text.x = element_text(size = 12 * text_scale, face = "bold", color = "black"),
      axis.text.y = element_text(size = 12 * text_scale, face = "bold", color = "black"),
      axis.title = element_text(size = 15 * text_scale, face = "bold", color = "black"),
      axis.title.x = element_text(size = 15 * text_scale, face = "bold", color = "black"),
      
      # Plot title - centered on top
      plot.title = element_text(size = 15 * text_scale, face = "bold", color = "black", hjust = 0.5),
      
      # Legend styling - stacked in bottom right corner inside plot
      legend.position = c(0.94, 0.02),
      legend.justification = c("right", "bottom"),
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 8 * text_scale, face = "bold", color = "black"),
      legend.key.size = unit(0.8, "cm"),
      legend.key.width = unit(0.8, "cm"),
      legend.key.height = unit(0.35, "cm"),
      legend.spacing.y = unit(0.3, "cm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      
      # Panel margins - matching diverging bars
      plot.margin = margin(t = 20, r = 10, b = 20, l = 10)
    )
  
  return(p)
}
