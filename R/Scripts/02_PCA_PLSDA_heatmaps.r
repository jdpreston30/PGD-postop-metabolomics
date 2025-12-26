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
#- 2.1.4: Generate all PLSDAs and PCAs
#- 2.1.4: Generate all PLSDAs and PCAs
# Check for cached permutation results in OneDrive
permutation_dir <- config$paths$permutation_results
if (file.exists(file.path(permutation_dir, "plsda_12h.rds")) && 
    file.exists(file.path(permutation_dir, "plsda_24h.rds")) && 
    file.exists(file.path(permutation_dir, "plsda_combined.rds"))) {
  message("Loading cached PLS-DA permutation results from OneDrive...")
  plsda_12h <- readRDS(file.path(permutation_dir, "plsda_12h.rds"))
  plsda_24h <- readRDS(file.path(permutation_dir, "plsda_24h.rds"))
  plsda_combined <- readRDS(file.path(permutation_dir, "plsda_combined.rds"))
} else {
  message("Running PLS-DA with permutation testing (this will take time)...")
  plsda_12h <- make_PCA(UFT_12h, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PLSDA", run_permutation = TRUE)
  plsda_24h <- make_PCA(UFT_24h, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PLSDA", run_permutation = TRUE)
  plsda_combined <- make_PCA(UFT_12and24, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PLSDA", run_permutation = TRUE)
}
saveRDS(plsda_combined, file.path(permutation_dir, "plsda_combined.rds"))
# PCAs for supplemental
pca_12h <- make_PCA(UFT_12h, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PCA")
pca_24h <- make_PCA(UFT_24h, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PCA")
pca_combined <- make_PCA(UFT_12and24, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PCA")
