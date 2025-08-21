#+ 3.5: Run t-tests for pairwise variant comparisons
#- 3.5.1: FV-PTC vs PTC comparison
# _Extract variant assignments (excluding FTC)
variant_assignments_FVPTC_PTC <- UFT_metaboanalyst_log2 %>%
  dplyr::filter(Variant != "FTC") %>%
  dplyr::select(Patient_ID, Variant)
FVPTC_vs_PTC <- mummichog_ttests(
  data = UFT_metaboanalyst_log2 %>% dplyr::filter(Variant != "FTC"),
  group_assignments = variant_assignments_FVPTC_PTC,
  group_column = "Variant",
  output_filename = "FVPTC_vs_PTC.csv",
  group1_value = "PTC",
  group2_value = "FV-PTC"
)
#- 3.5.2: FTC vs PTC comparison
# _Extract variant assignments (excluding FV-PTC)
variant_assignments_FTC_PTC <- UFT_metaboanalyst_log2 %>%
  dplyr::filter(Variant != "FV-PTC") %>%
  dplyr::select(Patient_ID, Variant)
FTC_vs_PTC <- mummichog_ttests(
  data = UFT_metaboanalyst_log2 %>% dplyr::filter(Variant != "FV-PTC"),
  group_assignments = variant_assignments_FTC_PTC,
  group_column = "Variant",
  output_filename = "FTC_vs_PTC.csv",
  group1_value = "PTC",
  group2_value = "FTC"
)
#- 3.5.3: FTC vs FV-PTC comparison
# _Extract variant assignments (excluding PTC)
variant_assignments_FTC_FVPTC <- UFT_metaboanalyst_log2 %>%
  dplyr::filter(Variant != "PTC") %>%
  dplyr::select(Patient_ID, Variant)
FTC_vs_FVPTC <- mummichog_ttests(
  data = UFT_metaboanalyst_log2 %>% dplyr::filter(Variant != "PTC"),
  group_assignments = variant_assignments_FTC_FVPTC,
  group_column = "Variant",
  output_filename = "FTC_vs_FVPTC.csv",
  group1_value = "FV-PTC",
  group2_value = "FTC"
)
#+ 3.6: Pathway Enrichment Plot
#- 3.6.1: Import results from mummichog
FVPTC_PTC <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FVPTC_PTC") %>%
  mutate(Comparisons = "FVPTC_PTC")
FTC_PTC <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FTC_PTC") %>%
  mutate(Comparisons = "FTC_PTC")
FVPTC_FTC <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FVPTC_FTC") %>%
  mutate(Comparisons = "FVPTC_FTC")
#- 3.6.2: Bind rows then filter to important variables
variant_enrichment <- bind_rows(FVPTC_PTC, FTC_PTC, FVPTC_FTC) %>%
  mutate(enrichment_factor = Hits.sig / Expected) %>%
  select(Comparisons, pathway_name, p_fisher, enrichment_factor) %>%
  filter(p_fisher < 0.05) %>%
  mutate(
    # Label comparisons nicely
    Comparisons = dplyr::case_when(
      Comparisons == "FVPTC_PTC" ~ "PTC vs. FV-PTC",
      Comparisons == "FTC_PTC" ~ "FTC vs. PTC",
      Comparisons == "FVPTC_FTC" ~ "FTC vs. FV-PTC",
      TRUE ~ Comparisons
    ),
    Comparisons = factor(Comparisons, levels = c("FTC vs. FV-PTC", "FTC vs. PTC", "PTC vs. FV-PTC")),
    pathway_name = forcats::fct_reorder(pathway_name, enrichment_factor, .fun = max)
  ) %>%
  tidyr::complete(pathway_name, Comparisons) %>%
  mutate(enrichment_factor = pmin(enrichment_factor, 7)) %>%
  mutate(
    # First handle Beta/β
    pathway_name = stringr::str_replace_all(
      pathway_name,
      stringr::regex("\\bBeta[- ]Alanine\\b", ignore_case = TRUE),
      "β-Alanine"
    ),
    # Capitalization fixes
    pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bmetabolism\\b", TRUE), "Metabolism"),
    pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bleucine\\b", TRUE), "Leucine"),
    pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bisoleucine\\b", TRUE), "Isoleucine"),
    pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bdegradation\\b", TRUE), "Degradation"),
    pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bmannose\\b", TRUE), "Mannose"),
    pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bacid\\b", TRUE), "Acid"),
    pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\boxidation\\b", TRUE), "Oxidation"),
    pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bperoxisome\\b", TRUE), "Peroxisome")
  ) %>%
  mutate(
    pathway_name = stringr::str_replace_all(
      pathway_name,
      stringr::regex("\\bretinol\\b", ignore_case = TRUE),
      "Retinol"
    )
  ) %>%
  mutate(
    pathway_name = factor(
      pathway_name,
      levels = filter(., Comparisons == "FTC vs. FV-PTC") %>%
        arrange(desc(enrichment_factor)) %>%
        pull(pathway_name) %>%
        unique()
    )
  )
#- 3.6.3: Plot
variant_enrichment_plot <- ggplot(
  variant_enrichment,
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
  coord_fixed(clip = "off") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +

  # Keep limits ascending; reverse legend order via guide
  scale_size_continuous(
    range = c(5, 20),
    limits = c(0, 7),
    breaks = c(7, 5, 3, 1), # Labels show in this order…
    name = "Enrichment factor",
    guide = guide_legend(reverse = TRUE) # …because we reverse the legend
  ) +

  # Keep p limits ascending; reverse colorbar via guide
  scale_color_gradient(
    low = "#0a2256", high = "#c3dbe9", # Dark (small p) -> light (large p)
    limits = c(0.01, 0.05),
    oob = scales::squish,
    name = NULL,
    guide = guide_colorbar(
      reverse = TRUE, # 0.01 at top, 0.05 at bottom
      barheight = unit(5, "cm"),
      barwidth = unit(0.9, "cm")
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
      face = "bold", family = "Arial", size = 6
    ),
    strip.text.y.left = element_text(
      angle = 0, hjust = 1,
      face = "bold", family = "Arial", size = 14,
      margin = margin(r = 6)
    ),
    legend.title = element_text(size = 14, face = "bold", family = "Arial"),
    legend.text = element_text(size = 14, family = "Arial"),
    plot.margin = margin(t = 20, r = 40, b = 10, l = 40)
  ) +
  coord_cartesian(clip = "off")
#- 3.6.4: Export Plot
panel_size <- 0.3 # tweak this until happy
n_rows <- length(unique(variant_enrichment$pathway_name))
n_cols <- length(unique(variant_enrichment$Comparisons))
# Add space for legends, strip labels, margins
extra_width <- 3 # inches for legends on the right
extra_height <- 1 # inches for titles/margins
total_width <- n_cols * panel_size + extra_width
total_height <- n_rows * panel_size + extra_height

ggsave(
  "variant_enrichment_plot.svg",
  variant_enrichment_plot,
  width = total_width,
  height = total_height,
  units = "in"
)
