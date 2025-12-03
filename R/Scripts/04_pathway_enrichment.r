#* 4: Pathway Enrichment Analysis
#+ 4.1: Run Mummichog on three LIMMA contrasts
#- 4.1.1: Run Mummichog (MFN) - Group/PGD effect
mummichog_pgd_mfn <- run_mummichog_analysis(
  ttest_results = limma_untarg$mummichog$group_main,
  output_dir = "Outputs/Enrichment/PGD",
  database = "hsa_mfn",
  instrumentOpt = 5.0,
  msModeOpt = "mixed",
  force_primary_ion = "yes"
)
#- 4.1.2: Run Mummichog (MFN) - Time effect
mummichog_time_mfn <- run_mummichog_analysis(
  ttest_results = limma_untarg$mummichog$time_main,
  output_dir = "Outputs/Enrichment/Time",
  database = "hsa_mfn",
  instrumentOpt = 5.0,
  msModeOpt = "mixed",
  force_primary_ion = "yes"
)
#- 4.1.3: Run Mummichog (MFN) - Interaction effect
mummichog_interaction_mfn <- run_mummichog_analysis(
  ttest_results = limma_untarg$mummichog$interaction,
  output_dir = "Outputs/Enrichment/Interaction",
  database = "hsa_mfn",
  instrumentOpt = 5.0,
  msModeOpt = "mixed",
  force_primary_ion = "yes"
)
#+ 4.2: Extract pathway results from CSV outputs
#- 4.2.1: Extract PGD pathways from CSV
PGD_MFN <- read_csv(
  list.files("Outputs/Enrichment/PGD", 
             pattern = "mummichog_pathway_enrichment_mummichog.csv", 
             full.names = TRUE, 
             recursive = TRUE)[1],
  show_col_types = FALSE
) |>
  select(pathway_name = 1, `Pathway total`, Hits.total, Hits.sig, Expected, `P(Fisher)`) |>
  mutate(Comparisons = "PGD")
#- 4.2.2: Extract Time pathways from CSV
time_MFN <- read_csv(
  list.files("Outputs/Enrichment/Time", 
             pattern = "mummichog_pathway_enrichment_mummichog.csv", 
             full.names = TRUE, 
             recursive = TRUE)[1],
  show_col_types = FALSE
) |>
  select(pathway_name = 1, `Pathway total`, Hits.total, Hits.sig, Expected, `P(Fisher)`) |>
  mutate(Comparisons = "Time")
#- 4.2.3: Extract Interaction pathways from CSV
interaction_MFN <- read_csv(
  list.files("Outputs/Enrichment/Interaction", 
             pattern = "mummichog_pathway_enrichment_mummichog.csv", 
             full.names = TRUE, 
             recursive = TRUE)[1],
  show_col_types = FALSE
) |>
  select(pathway_name = 1, `Pathway total`, Hits.total, Hits.sig, Expected, `P(Fisher)`) |>
  mutate(Comparisons = "Interaction")
#+ 4.3: Prep data for enrichment plots
MFN_enrichment <- bind_rows(PGD_MFN, time_MFN, interaction_MFN) |>
  rename(p_fisher = "P(Fisher)") |>
  mutate(enrichment_factor = Hits.sig / Expected) |>
  select(Comparisons, pathway_name, p_fisher, enrichment_factor, Hits.sig, Expected) |>
  rowwise() |>
  mutate(
    pathway_name = {
      words <- str_split(pathway_name, "\\s+")[[1]]
      new_name <- paste(sapply(words, function(w) {
        w_clean <- str_to_title(w)
        if (tolower(w) %in% c("and", "from")) {
          return(tolower(w))
        }
        if (tolower(w) == "epa") {
          return("EPA")
        }
        return(w_clean)
      }), collapse = " ")
      new_name
    }
  ) |>
  ungroup() |>
  filter(p_fisher <= 0.05) |>
  mutate(
    Comparisons = factor(Comparisons, levels = c("PGD", "Time", "Interaction")),
    pathway_name = forcats::fct_reorder(pathway_name, enrichment_factor, .fun = max)
  ) |>
    complete(pathway_name, Comparisons) |>
    mutate(enrichment_factor = pmin(enrichment_factor, 5)) |>
    mutate(
      pathway_name = factor(
        pathway_name,
        levels = filter(., Comparisons == "PGD") |>
          arrange(desc(enrichment_factor)) |>
          pull(pathway_name) |>
          unique()
      )
    ) |>
    filter(!is.na(p_fisher)) |>
    {
      df <- .
      # order by Combined, then append the rest (so nothing becomes NA)
      ordered <- df |>
        dplyr::filter(Comparisons == "PGD") |>
        dplyr::arrange(dplyr::desc(enrichment_factor)) |>
        dplyr::pull(pathway_name) |>
        as.character() |>
        unique()

      all_names <- unique(as.character(df$pathway_name))
      levels_all <- c(ordered, setdiff(all_names, ordered))

      df |>
        dplyr::mutate(pathway_name = factor(as.character(pathway_name),
          levels = levels_all
        ))
    }
#+ 4.4: Plot
#- 4.4.2: Plot
MFN_plot <- ggplot(
  MFN_enrichment,
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
    face = "bold", family = "Arial", size = 13
  ),
  strip.text.y.left = element_text(
    angle = 0, hjust = 1,
    face = "bold", family = "Arial", size = 13,
    margin = margin(r = 6)
  ),
  legend.title = element_text(size = 14, face = "bold", family = "Arial"),
  legend.text = element_text(size = 14, family = "Arial"),
  plot.margin = margin(t = 20, r = 40, b = 10, l = 40)
) +
coord_cartesian(clip = "off")
#- 4.4.3: Export Plot
{
  panel_size <- 0.3 # tweak this until happy
  n_rows <- length(unique(MFN_enrichment$pathway_name))
  n_cols <- length(unique(MFN_enrichment$Comparisons))

  # Add space for legends, strip labels, margins
  extra_width <- 9 # inches for legends on the right
  extra_height <- 6 # inches for titles/margins
  total_width <- n_cols * panel_size + extra_width
  total_height <- n_rows * panel_size + extra_height

  ggsave(
    "Outputs/Enrichment/MFN_enrichment_plot.png",
    MFN_plot,
    width = total_width,
    height = total_height,
    units = "in",
    dpi = 1200,
    bg = "white"
  )
}
