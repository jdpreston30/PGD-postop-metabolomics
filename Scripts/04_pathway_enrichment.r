#+ 3.1: Run t-tests for pairwise variant comparisons
  #- 3.5.0: Set up data for comparisons
    combined_UFT_mummichog <- combined_UFT %>%
      rename(PGD = Clinical_PGD, Patient_ID = Patient)
  #- 3.5.1: 12h PGD comparison
    PGD_assign_12 <- combined_UFT_mummichog %>%
      filter(Time != "24") %>%
      select(Patient_ID, PGD)
    PGD_12h_ttests <- mummichog_ttests(
      data = combined_UFT_mummichog %>% filter(Time != "24"),
      group_assignments = PGD_assign_12,
      group_column = "PGD",
      output_filename = "PGD_12h_ttests.csv",
      group1_value = "N",
      group2_value = "Y"
    )
  #- 3.5.2: 24h Comparison
    PGD_assign_24 <- combined_UFT_mummichog %>%
      filter(Time == "24") %>%
      select(Patient_ID, PGD)
    PGD_24h_ttests <- mummichog_ttests(
      data = combined_UFT_mummichog %>% filter(Time == "24"),
      group_assignments = PGD_assign_24,
      group_column = "PGD",
      output_filename = "PGD_24h_ttests.csv",
      group1_value = "N",
      group2_value = "Y"
    )
  #- 3.5.3: Full Comparison
    PGD_assign_full <- combined_UFT_mummichog %>%
      select(Patient_ID, PGD) %>%
      unique()
    PGD_full_ttests <- mummichog_ttests(
      data = combined_UFT_mummichog,
      group_assignments = PGD_assign_full,
      group_column = "PGD",
      output_filename = "PGD_full_ttests.csv",
      group1_value = "N",
      group2_value = "Y"
    )
#+ 3.1: Pathway Enrichment Plot Import
  #- 3.6.1: Import results from mummichog (MFN)
    PGD_12h_MFN <- read_xlsx("Outputs/Mummichog Outputs/mummichog_outputs.xlsx", sheet = "12h") %>% mutate(Comparisons = "12h")
    PGD_24h_MFN <- read_xlsx("Outputs/Mummichog Outputs/mummichog_outputs.xlsx", sheet = "24h") %>% mutate(Comparisons = "24h")
    PGD_all_MFN <- read_xlsx("Outputs/Mummichog Outputs/mummichog_outputs.xlsx", sheet = "All") %>% mutate(Comparisons = "Combined")
  #- 3.6.1: Import results from mummichog (KEGG)
    PGD_12h_KEGG <- read_xlsx("Outputs/Mummichog Outputs/mummichog_outputs.xlsx", sheet = "12h_KEGG") %>% mutate(Comparisons = "12h")
    PGD_24h_KEGG <- read_xlsx("Outputs/Mummichog Outputs/mummichog_outputs.xlsx", sheet = "24h_KEGG") %>% mutate(Comparisons = "24h")
    PGD_all_KEGG <- read_xlsx("Outputs/Mummichog Outputs/mummichog_outputs.xlsx", sheet = "All_KEGG") %>% mutate(Comparisons = "Combined")
#+ 3.3: Prep data for enrichment plots
#- 3.3.0: Define abbreviations to rename amino acids
  amino_map <- c(
    "Valine" = "Val",
    "Leucine" = "Leu",
    "Isoleucine" = "Ile",
    "Arginine" = "Arg",
    "Proline" = "Pro",
    "Glutamate" = "Glu",
    "Glutathione" = "GSH",
    "Tyrosine" = "Tyr",
    "Cysteine" = "Cys",
    "Valine" = "Val"
  )
#- 3.3.1: Bind rows then filter to important variables (MFN)
  KEGG_enrichment <- bind_rows(PGD_12h_KEGG, PGD_24h_KEGG, PGD_all_KEGG) %>%
    rename(p_fisher = "P(Fisher)") %>%
      mutate(enrichment_factor = Hits.sig / Expected) %>%
      select(Comparisons, pathway_name, p_fisher, enrichment_factor) %>%
      mutate(
        pathway_name = str_split(pathway_name, "\\s+") %>%
          lapply(function(words) {
            sapply(words, function(w) {
              # normalize case first
              w_clean <- str_to_title(w)
              # rules
              if (tolower(w) %in% c("and", "from")) {
                return(tolower(w))
              }
              if (str_to_lower(w) == "epa") {
                return("EPA")
              }
              if (w_clean %in% names(amino_map)) {
                return(amino_map[[w_clean]])
              }
              return(w_clean)
            }) %>%
              paste(collapse = " ")
          }) %>%
          unlist()
      ) %>%
      filter(p_fisher < 0.05) %>%
      mutate(
        Comparisons = factor(Comparisons, levels = c("12h", "24h", "Combined")),
        pathway_name = forcats::fct_reorder(pathway_name, enrichment_factor, .fun = max)
      ) %>%
      complete(pathway_name, Comparisons) %>%
      mutate(enrichment_factor = pmin(enrichment_factor, 7)) %>%
      mutate(
        pathway_name = factor(
          pathway_name,
          levels = filter(., Comparisons == "Combined") %>%
            arrange(desc(enrichment_factor)) %>%
            pull(pathway_name) %>%
            unique()
        )
      ) %>%
      filter(!is.na(p_fisher))
#- 3.3.2: BInd rows then filter to important variables (MFN)
  MFN_enrichment <- bind_rows(PGD_12h_MFN, PGD_24h_MFN, PGD_all_MFN) %>%
    rename(p_fisher = "P(Fisher)") %>%
    mutate(enrichment_factor = Hits.sig / Expected) %>%
    select(Comparisons, pathway_name, p_fisher, enrichment_factor) %>%
    mutate(
      pathway_name = str_split(pathway_name, "\\s+") %>%
        lapply(function(words) {
          sapply(words, function(w) {
            # normalize case first
            w_clean <- str_to_title(w)
            # rules
            if (tolower(w) %in% c("and", "from")) {
              return(tolower(w))
            }
            if (str_to_lower(w) == "epa") {
              return("EPA")
            }
            if (w_clean %in% names(amino_map)) {
              return(amino_map[[w_clean]])
            }
            return(w_clean)
          }) %>%
            paste(collapse = " ")
        }) %>%
        unlist()
    ) %>%
    filter(p_fisher < 0.05) %>%
    mutate(
      Comparisons = factor(Comparisons, levels = c("12h", "24h", "Combined")),
      pathway_name = forcats::fct_reorder(pathway_name, enrichment_factor, .fun = max)
    ) %>%
    complete(pathway_name, Comparisons) %>%
    mutate(enrichment_factor = pmin(enrichment_factor, 7)) %>%
    mutate(
      pathway_name = factor(
        pathway_name,
        levels = filter(., Comparisons == "Combined") %>%
          arrange(desc(enrichment_factor)) %>%
          pull(pathway_name) %>%
          unique()
      )
    ) %>%
    filter(!is.na(p_fisher)) %>%
    mutate(
      pathway_name = if_else(
        pathway_name == "Putative Anti-Inflammatory Metabolites Formation from EPA",
        "Anti-Inflammatory Metab. Formation from EPA",
        pathway_name
      ),
      pathway_name = if_else(
        pathway_name == "Ascorbate (Vitamin C) and Aldarate Metabolism",
        "Ascorbate and Aldarate Metabolism",
        pathway_name
      ),
      pathway_name = if_else(
        pathway_name == "Vitamin D3 (Cholecalciferol) Metabolism",
        "Cholecalciferol Metabolism",
        pathway_name
      )
    ) %>%
    {
      df <- .
      # order by Combined, then append the rest (so nothing becomes NA)
      ordered <- df %>%
        dplyr::filter(Comparisons == "Combined") %>%
        dplyr::arrange(dplyr::desc(enrichment_factor)) %>%
        dplyr::pull(pathway_name) %>%
        as.character() %>%
        unique()

      all_names <- unique(as.character(df$pathway_name))
      levels_all <- c(ordered, setdiff(all_names, ordered))

      df %>%
        dplyr::mutate(pathway_name = factor(as.character(pathway_name),
          levels = levels_all
        ))
    }
#+ 3.4: Plot
  #- 3.4.0: Set conflicts
    conflicts_prefer(ggplot2::margin)
  #- 3.6.3: Plot
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
  #- 3.6.4: Export Plot
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
    "MFN_enrichment_plot.png",
    MFN_plot,
    width = total_width,
    height = total_height,
    units = "in",
    dpi = 600
  )
}
