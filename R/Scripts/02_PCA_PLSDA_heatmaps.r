#*2: PCA, PLSDA, and Heatmaps
#+ 2.1: Run Exploratory PLSDA and PCA
#- 2.1.1: Define datasets
UFT_12h <- UFT |> 
  filter(Time == 12) |>
  select(Sample_ID, Time, severe_PGD, all_of(untargeted_features))
UFT_24h <- UFT |>
  filter(Time == 24) |>
  select(Sample_ID, Time, severe_PGD, all_of(untargeted_features))
UFT_12and24 <- UFT |>
  select(Sample_ID, Time, severe_PGD, all_of(untargeted_features))
#- 2.1.2: Define colors
cluster_colors <- c(
  "Y" = "#94001E",
  "N" = "#03507D"
)
#- 2.1.3: Set Plot Specs
plot_specs <- tribble(
  ~data, ~comp_x, ~comp_y, ~method, ~outpath,
  UFT_12h, 1, 2, "PLSDA", "Outputs/PLSDA/PLSDA_12h_ComparePGD.png",
  UFT_24h, 1, 2, "PLSDA", "Outputs/PLSDA/PLSDA_24h_ComparePGD.png",
  UFT_12and24, 1, 2, "PLSDA", "Outputs/PLSDA/PLSDA_combinedTime_ComparePGD.png",
  UFT_12h, 1, 2, "PCA", "Outputs/PCA/PCA_12h_ComparePGD.png",
  UFT_24h, 1, 2, "PCA", "Outputs/PCA/PCA_24h_ComparePGD.png",
  UFT_12and24, 1, 2, "PCA", "Outputs/PCA/PCA_combinedTime_ComparePGD.png"
)
#- 2.1.4: Load or generate permutation test results
{
  # Check if all cached RDS files exist
  cached_files <- c(
    config$paths$permutation$plsda_12h,
    config$paths$permutation$plsda_24h,
    config$paths$permutation$plsda_combined
  )
  if (all(file.exists(cached_files))) {
    message("Loading cached permutation test results...")
    plsda_12h <- readRDS(config$paths$permutation$plsda_12h)
    plsda_24h <- readRDS(config$paths$permutation$plsda_24h)
    plsda_combined <- readRDS(config$paths$permutation$plsda_combined)
  } else {
    plsda_12h <- make_PCA(UFT_12h, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PLSDA", run_permutation = TRUE)
    plsda_24h <- make_PCA(UFT_24h, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PLSDA", run_permutation = TRUE)
    plsda_combined <- make_PCA(UFT_12and24, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PLSDA", run_permutation = TRUE)
    # Save combined results as list (12h and 24h already saved after generation)
    if (!file.exists(config$paths$permutation$plsda_combined)) {
      saveRDS(plsda_combined, config$paths$permutation$plsda_combined)
    }
  }
}
#- 2.1.5: Generate plots only (without permutation testing)
plsda_12h$plot <- make_PCA(UFT_12h, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PLSDA", run_permutation = FALSE)$plot
plsda_24h$plot <- make_PCA(UFT_24h, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PLSDA", run_permutation = FALSE)$plot
plsda_combined$plot <- make_PCA(UFT_12and24, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PLSDA", run_permutation = FALSE)$plot
# PCAs for supplemental
pca_12h <- make_PCA(UFT_12h, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PCA")
pca_24h <- make_PCA(UFT_24h, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PCA")
pca_combined <- make_PCA(UFT_12and24, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PCA")
