volc_targ <- function(FC_list, ttest_res, key){
      # Configure FC_table variable to merge fold-change of each metabolite with t-test results
      FC_table <- data.frame(
      Metabolite = unlist(rownames(FC_list)),
      Value = unlist(FC_list$Value),
      Time = unlist(FC_list$Time),
      row.names = NULL
      )
      FC_table$Time <- as.character(FC_table$Time)
      ttest_res$Time <- as.character(ttest_res$Time)
      # Merge t-test results with fold-change data into the volcano_data variable for plotting
      volcano_data <- merge(FC_table, ttest_res, by = c("Metabolite", "Time"), all = TRUE)
      # View(volcano_data)
      volcano_data <- volcano_data %>%
        left_join(key, by = c("Metabolite" = "Name"))
      # View(volcano_data)
      # Add a column for color coding
      volcano_data <- volcano_data %>%
        mutate(Legend = case_when(
          p_value < 0.05 & Value > log2(1.5) ~ "Up in PGD",
          p_value < 0.05 & Value < -log2(1.5) ~ "Down in PGD",
          TRUE ~ "Not Significant"
        ))
      # Add a column for significance, used for labeling significant metabolites
      # volcano_data$Label <- ifelse(volcano_data$Legend != 'Not Significant', 'true', 'false')
      top10sig <- volcano_data %>%
        filter(p_value < 0.05) %>%
        arrange(Value) %>%
        slice_head(n = 10)
      bottom10sig <- volcano_data %>%
        filter(p_value < 0.05) %>%
        arrange(desc(Value)) %>%
        slice_head(n = 10)
      # Combine top and bottom metabolites for labeling
      top_bottom <- rbind(top10sig, bottom10sig)
      # Add a column for labeling metabolites
      volcano_data$Label <- ifelse(volcano_data$Metabolite %in% top_bottom$Metabolite, TRUE, FALSE)

      # Create the volcano plot
      volcano_plot <- ggplot(volcano_data, aes(x = Value, y = -log10(p_value), color = Legend)) +
        geom_point(size = 2, alpha = 0.7) +
        scale_color_manual(values = c('Not Significant' = "gray", 'Up in PGD' = '#800017', 'Down in PGD' = '#113d6a')) +
        theme_light(base_family = "Arial") +
        labs(x = "Log2 Fold Change",
            y = "-Log10(p-value)") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        scale_x_continuous(
          limits = c(-3, 3),
          breaks = seq(-3, 3, 1),
          minor_breaks = seq(-3, 3, 1)
        ) +
        scale_y_continuous(
          limits = c(-0.5, 4),
          breaks = seq(0, 4, 1),
          minor_breaks = seq(-0.5, 4, 0.5)
        ) +
        geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "black") +
        theme(
          legend.position = "right",
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 13, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          ) +
          guides(color = guide_legend(override.aes = list(shape = 16, size=4)))
        return(list(
            volcano_data = volcano_data,
            volcano_plot = volcano_plot))
    }