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
) %>%
  select(pathway_name = 1, `Pathway total`, Hits.total, Hits.sig, Expected, `P(Fisher)`) %>%
  mutate(Comparisons = "PGD")
#- 4.2.2: Extract Time pathways from CSV
time_MFN <- read_csv(
  list.files("Outputs/Enrichment/Time", 
             pattern = "mummichog_pathway_enrichment_mummichog.csv", 
             full.names = TRUE, 
             recursive = TRUE)[1],
  show_col_types = FALSE
) %>%
  select(pathway_name = 1, `Pathway total`, Hits.total, Hits.sig, Expected, `P(Fisher)`) %>%
  mutate(Comparisons = "Time")
#- 4.2.3: Extract Interaction pathways from CSV
interaction_MFN <- read_csv(
  list.files("Outputs/Enrichment/Interaction", 
             pattern = "mummichog_pathway_enrichment_mummichog.csv", 
             full.names = TRUE, 
             recursive = TRUE)[1],
  show_col_types = FALSE
) %>%
  select(pathway_name = 1, `Pathway total`, Hits.total, Hits.sig, Expected, `P(Fisher)`) %>%
  mutate(Comparisons = "Interaction")
#+ 4.3: Prep data for enrichment plots
MFN_enrichment <- bind_rows(PGD_MFN, time_MFN, interaction_MFN) %>%
  rename(p_fisher = "P(Fisher)") %>%
  mutate(enrichment_factor = Hits.sig / Expected) %>%
  select(Comparisons, pathway_name, p_fisher, enrichment_factor, Hits.sig, Expected) %>%
  rowwise() %>%
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
  ) %>%
  ungroup() %>%
  filter(p_fisher <= 0.05) %>%
  mutate(
    Comparisons = factor(Comparisons, levels = c("PGD", "Time", "Interaction")),
    pathway_name = forcats::fct_reorder(pathway_name, enrichment_factor, .fun = max)
  ) %>%
    complete(pathway_name, Comparisons) %>%
    mutate(enrichment_factor = pmin(enrichment_factor, 5)) %>%
    mutate(
      pathway_name = factor(
        pathway_name,
        levels = filter(., Comparisons == "PGD") %>%
          arrange(desc(enrichment_factor)) %>%
          pull(pathway_name) %>%
          unique()
      )
    ) %>%
    filter(!is.na(p_fisher)) %>%
    {
      df <- .
      # order by Combined, then append the rest (so nothing becomes NA)
      ordered <- df %>%
        dplyr::filter(Comparisons == "PGD") %>%
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
#+ 4.4: Plot enrichment function
MFN_plot <- plot_enrichment(MFN_enrichment)
