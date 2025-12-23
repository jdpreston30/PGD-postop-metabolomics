 TFT_combined_succinate <- TFT_combined_LIMMA |>
  select(Patient, Time, severe_PGD, HILIC_162.9978_91.9)

# Faceted line graph
succinate_plot <- ggplot(TFT_combined_succinate, aes(x = Time, y = HILIC_162.9978_91.9, group = Patient, color = severe_PGD)) +
  geom_line(alpha = 0.5, linewidth = 0.8) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = c("N" = "#03507D", "Y" = "#94001E"), guide = "none") +
  facet_wrap(~ severe_PGD, labeller = labeller(severe_PGD = c("N" = "No Severe PGD", "Y" = "Severe PGD"))) +
  labs(
    x = "Time (hours)",
    y = "Succinate",
    title = "Succinate Levels by Time and PGD Status"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
print_to_png(succinate_plot, "succinate_facet_plot.png", dpi = 600, width = 4, height = 3)

# Calculate averages for summary plot
TFT_succinate_summary <- TFT_combined_succinate |>
  group_by(severe_PGD, Time) |>
  summarise(
    mean_succinate = mean(HILIC_162.9978_91.9, na.rm = TRUE),
    se_succinate = sd(HILIC_162.9978_91.9, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Summary line graph with error bars
succinate_summary_plot <- ggplot(TFT_succinate_summary, aes(x = Time, y = mean_succinate, group = severe_PGD, color = severe_PGD)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_succinate - se_succinate, ymax = mean_succinate + se_succinate), 
                width = 0.1, linewidth = 0.8) +
  scale_color_manual(
    values = c("N" = "#03507D", "Y" = "#94001E"),
    labels = c("N" = "No Severe PGD", "Y" = "Severe PGD"),
    name = NULL
  ) +
  labs(
    x = "Time (hours)",
    y = "Succinate (mean ± SE)",
    title = "Average Succinate Levels by Time and PGD Status"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "top",
    legend.text = element_text(size = 10, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print_to_png(succinate_summary_plot, "succinate_summary_plot.png", dpi = 600, width = 4, height = 3)

# Calculate averages for summary plot (untransformed data)
TFT_succinate_summary_untransformed <- TFT_combined_succinate |>
  mutate(succinate_original = 2^HILIC_162.9978_91.9) |>
  group_by(severe_PGD, Time) |>
  summarise(
    mean_succinate = mean(succinate_original, na.rm = TRUE),
    se_succinate = sd(succinate_original, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Summary line graph with error bars (untransformed)
succinate_summary_plot_untransformed <- ggplot(TFT_succinate_summary_untransformed, 
                                                aes(x = Time, y = mean_succinate, group = severe_PGD, color = severe_PGD)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_succinate - se_succinate, ymax = mean_succinate + se_succinate), 
                width = 0.1, linewidth = 0.8) +
  scale_color_manual(
    values = c("N" = "#03507D", "Y" = "#94001E"),
    labels = c("N" = "No Severe PGD", "Y" = "Severe PGD"),
    name = NULL
  ) +
  labs(
    x = "Time (hours)",
    y = "Succinate (original scale, mean ± SE)",
    title = "Average Succinate Levels (Untransformed)"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "top",
    legend.text = element_text(size = 10, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print_to_png(succinate_summary_plot_untransformed, "succinate_summary_plot_untransformed.png", dpi = 600, width = 4, height = 3)

# Calculate averages in original scale (undo log2 transform) and prepare individual points
TFT_succinate_untransformed <- TFT_combined_succinate |>
  mutate(succinate_original = 2^HILIC_162.9978_91.9)

TFT_succinate_means <- TFT_succinate_untransformed |>
  group_by(severe_PGD, Time) |>
  summarise(
    mean_succinate = mean(succinate_original, na.rm = TRUE),
    .groups = "drop"
  )

# Faceted bar plot with original scale
succinate_bar_plot <- ggplot() +
  geom_col(data = TFT_succinate_means, 
           aes(x = Time, y = mean_succinate, fill = severe_PGD, color = severe_PGD),
           width = 0.7, linewidth = 0.8, alpha = 0.5) +
  geom_point(data = TFT_succinate_untransformed,
             aes(x = Time, y = succinate_original, color = severe_PGD),
             size = 2.5, alpha = 0.7, position = position_jitter(width = 0.1, height = 0, seed = 42)) +
  scale_fill_manual(values = c("N" = "#03507D", "Y" = "#94001E"), guide = "none") +
  scale_color_manual(values = c("N" = "#03507D", "Y" = "#94001E"), guide = "none") +
  facet_wrap(~ severe_PGD, labeller = labeller(severe_PGD = c("N" = "No Severe PGD", "Y" = "Severe PGD"))) +
  labs(
    x = "Time (hours)",
    y = "Succinate (original scale)",
    title = "Average Succinate Levels (Untransformed)"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print_to_png(succinate_bar_plot, "succinate_bar_plot.png", dpi = 600, width = 4, height = 3)
