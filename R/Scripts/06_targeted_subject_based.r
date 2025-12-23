#* 6: Subject-specific Analysis
#+ 6.1: Clean up interaction metabolites for graphing
#- 6.1.1: Get metadata to join
QC_intx <- read_xlsx(config$paths$manual_QC, sheet = "intx_curated")
intx_names_curated <- QC_intx |>
  select(feature, identified_name, display_name, main_group, subgroup)
#- 6.1.2: Leftjoin with limma results
intx_display_data <- intx_names_curated |>
  left_join(limma_targ$interaction |>
              rename(feature = Metabolite),
            by = "feature") |>
  select(feature, identified_name, display_name, main_group, subgroup, p_value = p.value, log2_fc = logFC)
#+ 6.2: Graph
div_bars_intx <- plot_diverging_bars(intx_display_data, 
  group_ordering = TRUE, 
  add_group_labels = TRUE,
  max_features = 50,  
  fc_threshold = 0,
  x_max = 5,
  label_pos = 1
)
ggsave(
  "test.png",
  div_bars_intx,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)
