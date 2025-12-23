plot_enrichment <- function(enrichment_data) {
  # Scale factor: 8/14 ≈ 0.571 to match other plots
  scale_factor <- 8/14
  
  # Rename specific pathway
  enrichment_data <- enrichment_data %>%
    mutate(pathway_name = if_else(
      as.character(pathway_name) == "Mono-Unsaturated Fatty Acid Beta-Oxidation",
      "MUFA Beta-Oxidation",
      as.character(pathway_name)
    )) %>%
    mutate(pathway_name = factor(pathway_name, levels = unique(pathway_name)))
  
  ggplot(
    enrichment_data,
    aes(x = 0.5, y = 0.5, size = enrichment_factor, color = p_fisher)
  ) +
  # One dummy row per facet -> avoid the warning
  geom_tile(
    data = data.frame(x = 0.5, y = 0.5),
    aes(x = x, y = y),
    width = 1, height = 1,
    fill = "white", colour = "grey80", linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  geom_point(
    alpha = 0.95, shape = 16, stroke = 0,
    na.rm = TRUE, show.legend = TRUE
  ) +
  facet_grid(
    rows = vars(pathway_name),
    cols = vars(Comparisons),
    switch = "y", drop = FALSE
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  
  # Keep limits ascending; reverse legend order via guide
  scale_size_continuous(
    range = c(5 * scale_factor, 20 * scale_factor),
    limits = c(0, 5),
    breaks = c(5, 3, 1),
    labels = c("5+", "3", "1"),
    name = "Enrichment factor",
    guide = guide_legend(reverse = TRUE) # …because we reverse the legend
  ) +
  # Keep p limits ascending; reverse colorbar via guide
  scale_color_gradient(
    low = "#0a2256", high = "#c3dbe9", # Dark (small p) -> light (large p)
    limits = c(0.01, 0.05),
    oob = scales::squish,
    name = "p-value\n",
    guide = guide_colorbar(
      reverse = TRUE, # 0.01 at top, 0.05 at bottom
      barheight = unit(5 * scale_factor, "cm"),
      barwidth = unit(0.9 * scale_factor, "cm")
    )
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing.x = unit(0, "pt"),
    panel.spacing.y = unit(0, "pt"),
    strip.placement = "outside",
    strip.text.x.top = element_text(
      angle = 0, vjust = 1,
      face = "bold", family = "Arial", size = 13 * scale_factor
    ),
    strip.text.y.left = element_text(
      angle = 0, hjust = 1,
      face = "bold", family = "Arial", size = 12 * scale_factor,
      margin = margin(r = 6 * scale_factor)
    ),
    legend.title = element_text(size = 8, face = "bold", family = "Arial"),
    legend.text = element_text(size = 8, face = "bold", color = "black", margin = margin(l = 4, r = 4)),
    plot.margin = margin(t = 20 * scale_factor, r = 40 * scale_factor, b = 10 * scale_factor, l = 120)
  ) +
  coord_cartesian(clip = "off")
}
