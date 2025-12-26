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
# Load individual RDS components from local Outputs/Permutation folder
message("Loading cached PLS-DA permutation results...")

# Load 12h components
if (all(file.exists(c(
  "Outputs/Permutation/plsda_12h_plot.rds",
  "Outputs/Permutation/plsda_12h_scores.rds",
  "Outputs/Permutation/plsda_12h_cv.rds",
  "Outputs/Permutation/plsda_12h_permutation.rds"
)))) {
  plsda_12h <- list(
    plot = readRDS("Outputs/Permutation/plsda_12h_plot.rds"),
    model = NULL,
    scores = readRDS("Outputs/Permutation/plsda_12h_scores.rds"),
    scores_df = readRDS("Outputs/Permutation/plsda_12h_scores_df.rds"),
    explained = readRDS("Outputs/Permutation/plsda_12h_explained.rds"),
    Y = readRDS("Outputs/Permutation/plsda_12h_Y.rds"),
    cv = readRDS("Outputs/Permutation/plsda_12h_cv.rds"),
    permutation = readRDS("Outputs/Permutation/plsda_12h_permutation.rds")
  )
} else {
  stop("Missing 12h permutation results. Run permutation test first.")
}

# Load 24h components
if (all(file.exists(c(
  "Outputs/Permutation/plsda_24h_plot.rds",
  "Outputs/Permutation/plsda_24h_scores.rds",
  "Outputs/Permutation/plsda_24h_cv.rds",
  "Outputs/Permutation/plsda_24h_permutation.rds"
)))) {
  plsda_24h <- list(
    plot = readRDS("Outputs/Permutation/plsda_24h_plot.rds"),
    model = NULL,
    scores = readRDS("Outputs/Permutation/plsda_24h_scores.rds"),
    scores_df = readRDS("Outputs/Permutation/plsda_24h_scores_df.rds"),
    explained = readRDS("Outputs/Permutation/plsda_24h_explained.rds"),
    Y = readRDS("Outputs/Permutation/plsda_24h_Y.rds"),
    cv = readRDS("Outputs/Permutation/plsda_24h_cv.rds"),
    permutation = readRDS("Outputs/Permutation/plsda_24h_permutation.rds")
  )
} else {
  stop("Missing 24h permutation results. Run permutation test first.")
}

# For combined, run the permutation test (not cached)
message("Running PLS-DA for combined 12h + 24h timepoint...")
plsda_combined <- make_PCA(UFT_12and24, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PLSDA", run_permutation = TRUE)
# PCAs for supplemental
pca_12h <- make_PCA(UFT_12h, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PCA")
pca_24h <- make_PCA(UFT_24h, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PCA")
pca_combined <- make_PCA(UFT_12and24, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PCA")
#++++ Manual
