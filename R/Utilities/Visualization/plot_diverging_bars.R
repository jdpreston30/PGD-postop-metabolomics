#' Create Diverging Bar Plot for Fold Changes
#'
#' This function creates a horizontal diverging bar plot showing log2 fold changes
#' for significant metabolomic features, with styling consistent with volcano plots.
#'
#' @param results_tibble Tibble with columns: display_name, fold_change
#' @param base_family Font family for plots (default: "Arial")
#' @param max_features Maximum number of features to display (default: 20)
#' @param order_by Column to order features by for selection (default: "p_value")
#' @param title Plot title (default: "Log2 Fold Change (Severe vs No Severe PGD)")
#' @param text_scale Scaling factor for all text elements (default: 1.0, use 0.8 for 80% size, etc.)
#' @param fc_threshold Minimum absolute fold change threshold to include features (default: 1.5)
#' @param group_ordering Logical; if TRUE and main_group column exists, orders by group then by FC within each group (positive FC descending, then negative FC ascending) (default: FALSE)
#' @param add_group_labels Logical; if TRUE and main_group column exists, adds group brackets and labels on y-axis (default: FALSE)
#' @param x_max Maximum value for x-axis (default: 3.5)
#'
#' @return ggplot object
#'
#' @examples
#' \dontrun{
#'   # Create diverging bar plot
#'   p <- plot_diverging_bars(
#'     results_tibble = inspect,
#'     max_features = 15,
#'     fc_threshold = 2.0
#'   )
#'   print(p)
#' }
#'
#' @export
plot_diverging_bars <- function(results_tibble,
                                base_family = "Arial",
                                max_features = 20,
                                order_by = "p_value",
                                text_scale = 1.0,
                                fc_threshold = 1.5,
                                group_ordering = FALSE,
                                add_group_labels = FALSE,
                                x_max = 3.5,
                                title = NULL) {
  
  # Load required libraries
  library(dplyr)
  library(ggplot2)
  library(stringr)
  
  # Prepare data
  plot_data <- results_tibble %>%
    # Use existing log2_FC/log2FC column if available, otherwise calculate from fold_change
    {
      if ("log2_fc" %in% names(.)) {
        # log2_fc already exists, use as-is
        .
      } else if ("log2FC" %in% names(.)) {
        # log2FC exists, rename to log2_fc
        rename(., log2_fc = log2FC)
      } else {
        # Neither exists, calculate from fold_change
        mutate(., log2_fc = log2(fold_change))
      }
    } %>%
    # Apply fold change threshold filter
    filter(abs(log2_fc) >= log2(fc_threshold)) %>%
    # Create color category based on fold change direction
    mutate(
      fc_direction = ifelse(log2_fc >= 0, "positive", "negative"),
      # Clean feature names for display and add +/- prefix
      display_name = str_replace_all(display_name, "_", " "),
      display_name = str_wrap(display_name, width = 40),  # Increased width to reduce line breaks
      # Convert to absolute values for plotting (positive x-axis only)
      log2_fc_abs = abs(log2_fc)
    ) %>%
    # Take top features first (by order_by if specified, otherwise all)
    {if(order_by != "log2_fc") arrange(., !!sym(order_by)) else .} %>%
    slice_head(n = max_features) %>%
    # Ordering logic: by group or by log2FC
    {
      if (group_ordering && "main_group" %in% names(.)) {
        # Group ordering: within each main_group, positive FCs descending (largest to smallest),
        # then negative FCs ascending (most negative to least negative)
        # Reverse at the end so reds appear on top
        arrange(., main_group, desc(fc_direction == "positive"), desc(log2_fc)) %>%
          mutate(display_name = factor(display_name, levels = rev(display_name)))
      } else {
        # Sort by log2FC: lowest to highest (original values, not absolute)
        arrange(., log2_fc) %>%
          mutate(display_name = factor(display_name, levels = display_name))
      }
    }
  
  # Color scheme: dark red for positive FC, dark blue for negative FC
  colors <- c(
    "positive" = "#800017",  # Dark red for positive fold changes
    "negative" = "#113d6a"   # Dark blue for negative fold changes
  )
  
  # Publication-style theme matching your existing plots
  theme_pub_diverging <- function() {
    theme_minimal(base_family = "Arial") +
      theme(
        # Panel styling - clean background with boxed border
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.545),
        
        # Axis styling - matching volcano plot exactly
        axis.ticks = element_line(color = "black", linewidth = 0.6),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 12 * text_scale, face = "bold", color = "black"),
        axis.text.y = element_text(size = 6 * text_scale, face = "bold", color = "black", margin = margin(r = 2)),
        axis.title = element_text(size = 15 * text_scale, face = "bold", color = "black"),
        axis.title.x = element_text(size = 15 * text_scale, face = "bold", color = "black"),
        axis.title.y = element_blank(),
        
        # Plot title
        plot.title = element_blank(),
        
        # Legend styling - centered at top, universal spacing
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.title = element_blank(),
        legend.text = element_text(size = 8 * text_scale, face = "bold", color = "black", margin = margin(l = 4, r = 4)),
        legend.key.size = unit(0.35, "cm"),
        legend.key.width = unit(0.35, "cm"),
        legend.key.height = unit(0.15, "cm"),
        legend.margin = margin(t = 0, r = 0, b = -3, l = 0),
        legend.box.margin = margin(0, 0, 0, 0),
        
        # Panel margins - give more space on the left for y-axis labels and right for group labels
        plot.margin = margin(t = 20, r = 10, b = 20, l = 10)
      )
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = log2_fc_abs, y = display_name, fill = fc_direction)) +
    # Add bars
    geom_col(width = 0.7, color = "black", linewidth = 0.3) +
    # Color and fill scales with descriptive labels
    scale_fill_manual(
      values = colors,
      labels = c("negative" = "Lower in sPGD", "positive" = "Higher in sPGD")
    ) +
    # Y-axis spacing control
    scale_y_discrete(expand = expansion(add = c(0.75, 0.75))) +
    # X-axis scale
    scale_x_continuous(
      limits = c(0, x_max),
      breaks = seq(0, floor(x_max), by = 1),
      expand = expansion(mult = c(0.0037, 0.0037))
    ) +
    # Axis labels
    labs(
      title = title,
      x = expression(bold("|log")[2]*bold("(Fold Change)|")),
      # y = "Metabolite"
    ) +
    # Apply theme
    theme_pub_diverging()
  
  # Add group dividing lines and labels if requested and main_group column exists
  if (add_group_labels && "main_group" %in% names(plot_data)) {
    # Calculate group positions
    group_info <- plot_data %>%
      mutate(y_pos = as.numeric(display_name)) %>%
      group_by(main_group) %>%
      summarise(
        y_min = min(y_pos),
        y_max = max(y_pos),
        y_mid = mean(c(y_min, y_max)),
        .groups = "drop"
      ) %>%
      arrange(y_min)
    
    # Add horizontal dividing lines between groups (but not before first or after last)
    if (nrow(group_info) > 1) {
      for (i in 1:(nrow(group_info) - 1)) {
        divider_y <- group_info$y_max[i] + 0.5
        p <- p +
          annotate("segment",
            x = 0, xend = x_max,
            y = divider_y, yend = divider_y,
            linewidth = 0.3, color = "black"
          )
      }
    }
    
    # Position labels in the right margin (outside plot area)
    # Use Inf to place at far right edge of plot panel
    
    # Calculate x position for labels (just inside the right edge of plot)
    x_max <- max(plot_data$log2_fc_abs, na.rm = TRUE)
    label_x <- x_max * 1.15  # Position at 95% of max value
    
    # Add group labels inside plot area with line breaks for multi-word labels
    for (i in 1:nrow(group_info)) {
      # Add line break if label has multiple words
      label_text <- gsub(" ", "\n", group_info$main_group[i])
      
      p <- p +
        annotate("text",
          x = label_x, y = group_info$y_mid[i],
          label = label_text,
          hjust = 0.5, vjust = 0.5,
          size = 2.5 * text_scale, fontface = "bold.italic",
          family = base_family, color = "black"
        )
    }
    
    # Turn off clipping to show text in margins
    p <- p + coord_cartesian(clip = "off")
  }
  
  return(p)
}