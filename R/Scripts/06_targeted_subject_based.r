#* 6: Subject-specific Analysis
#+ 6.1: Clean up interaction metabolites for graphing
#- 6.1.1: Get metadata to join
intx_names_curated <- QC_intx |>
  select(feature, identified_name, display_name, main_group, subgroup)
#- 6.1.2: Leftjoin with limma results
intx_display_data <- intx_names_curated |>
  left_join(limma_targ$interaction |>
              rename(feature = Metabolite),
            by = "feature") |>
  select(feature, identified_name, display_name, main_group, subgroup, p_value = p.value, log2_fc = logFC)
#+ 6.2: Graph Diverging Bars of Interaction Metabolites
div_bars_intx <- plot_diverging_bars(intx_display_data, 
  group_ordering = TRUE, 
  add_group_labels = TRUE,
  max_features = 50,  
  fc_threshold = 0,
  x_max = 5.05,
  lower_expand = 0.0001,
  label_pos = 0.9,
  legend_labels = "intx"
)
#+ 6.3: Prepare data for Succcinate and AMP Plot
#- 6.3.1: Succinate Data
TFT_combined_succinate <- TFT_combined_LIMMA |>
select(Patient, Time, severe_PGD, HILIC_162.9978_91.9) |>
group_by(severe_PGD, Time) |>
summarise(
  mean_succinate = mean(HILIC_162.9978_91.9, na.rm = TRUE),
  se_succinate = sd(HILIC_162.9978_91.9, na.rm = TRUE) / sqrt(n()),
  .groups = "drop"
)
#- 6.3.2: AMP Data
TFT_combined_AMP <- TFT_combined_LIMMA |>
select(Patient, Time, severe_PGD, HILIC_348.0705_262.3) |>
group_by(severe_PGD, Time) |>
summarise(
  mean_AMP = mean(HILIC_348.0705_262.3, na.rm = TRUE),
  se_AMP = sd(HILIC_348.0705_262.3, na.rm = TRUE) / sqrt(n()),
  .groups = "drop"
)
#+ 6.4: Graph Succinate Summary Plot
#- 6.4.1: Succinate
succinate_summary_plot <- plot_line_summary(
  data = TFT_combined_succinate,
  time_var = "Time",
  mean_var = "mean_succinate",
  error_var = "se_succinate",
  group_var = "severe_PGD",
  legend_position = "BR",
  title = "Succinate"
)
#- 6.4.2: AMP
AMP_summary_plot <- plot_line_summary(
  data = TFT_combined_AMP,
  time_var = "Time",
  mean_var = "mean_AMP",
  error_var = "se_AMP",
  group_var = "severe_PGD",
  legend_position = "TR",
  title = "AMP"
)
