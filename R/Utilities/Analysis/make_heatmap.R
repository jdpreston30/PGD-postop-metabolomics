#' Create heatmap with optional feature selection (ANOVA / variance / MAD)
#' and optional T_stage annotation. Returns plot object for patchwork.
#'
#' @param data                 Data frame with Sample_ID, Variant, optional T_stage, then feature columns
#' @param top_features         NULL (default: show all features). If numeric >0, keep top N by `feature_selector`.
#' @param feature_selector     One of c("none","anova","variance","mad"). Default "none".
#' @param p_anova       ANOVA comparison. Choose either 'pgd', 'time', or 'interaction'.
#' @param n_clades             Number of clades to extract from sample HCA (default 2)
#' @param file_path       File path for saving heatmap
#' @param file_name       File name for saving heatmap
#' @return List with plot object, M, Mz, hc_cols, clade_df, clade_lists, ann_col, ann_colors, etc.
#' @export
make_heatmap <- function(
    data,
    top_features = NULL,
    feature_selector = c("none", "anova", "variance", "mad", "ttest"),
    p_anova = NULL,
    n_clades = 2,
    file_path,
    file_name) {
  feature_selector <- match.arg(feature_selector)

  # ---- Checks ----
  stopifnot(all(c("Sample_ID", "Clinical_PGD") %in% names(data)))
  has_time <- "Time" %in% names(data)

  # Keep ID/Variant/(optional)T_stage up front
  dat <- dplyr::select(
    data,
    dplyr::any_of(c("Sample_ID", "Clinical_PGD", if (has_time) "Time")),
    dplyr::everything()
  )

  # Coerce Clinical_PGD to factor in desired order
  dat$Clinical_PGD <- factor(dat$Clinical_PGD, levels = unique(dat$Clinical_PGD))

  # Build matrix: samples x features (drop ID/Clinical_PGD/Time)
  drop_cols <- c("Sample_ID", "Clinical_PGD", "Time")
  X <- as.matrix(dplyr::select(dat, -dplyr::all_of(drop_cols)))
  stopifnot(all(vapply(as.data.frame(X), is.numeric, TRUE)))

  # Sample IDs (use make.names for uniqueness + valid rownames)
  sample_ids <- make.names(dat$Sample_ID, unique = TRUE)
  rownames(X) <- sample_ids

  # Group factor for ANOVA ranking
  group_pgd <- dat$Clinical_PGD
  group_time <- dat$Time

  # ---- Optional feature ranking/selection ----
  if (!is.null(top_features) && is.numeric(top_features) && top_features > 0 &&
    feature_selector != "none") {
    top_n <- min(top_features, ncol(X))

    rank_idx <- switch(feature_selector,
      "anova" = {
        pvals <- apply(X, 2, function(x) {
        fit <- aov(x ~ group_pgd * group_time)
        switch (p_anova, 
        "pgd" = {
          summary(fit)[[1]][["Pr(>F)"]][1]
        },
        "time" = {
          summary(fit)[[1]][["Pr(>F)"]][2]
        },
        "interaction" = {
          summary(fit)[[1]][["Pr(>F)"]][3]
        })
      })
      },
      "variance" = {
        v <- apply(X, 2, stats::var, na.rm = TRUE)
        order(v, decreasing = FALSE, na.last = NA)
      },
      "mad" = {
        m <- apply(X, 2, stats::mad, na.rm = TRUE)
        order(m, decreasing = FALSE, na.last = NA)
      },
      "ttest" = {
        pvals <- apply(X, 2, function(x) {
          # Run t-test between two groups
          ttest_res <- t.test(x ~ group_pgd)
          ttest_res$p.value
        })
        order(pvals, decreasing = FALSE, na.last = NA)
      }
    )

    X <- X[, head(rank_idx, top_n), drop = FALSE]
  }

  # Drop zero-variance columns (after selection)
  nzv <- apply(X, 2, sd, na.rm = TRUE) > 0
  nzv[is.na(nzv)] <- FALSE
  if (!all(nzv)) X <- X[, nzv, drop = FALSE]
  if (!ncol(X)) stop("No features remain after selection/variance filtering.")

  # Heatmap matrix: features x samples
  M <- t(X)

  # ---- Sample clustering & clades ----
  Mz <- t(scale(t(M), center = TRUE, scale = TRUE))
  Mz[is.na(Mz)] <- 0
  d_cols <- dist(t(Mz), method = "euclidean")
  hc_cols <- hclust(d_cols, method = "complete")

  # Map back to original Sample_ID
  id_map <- setNames(dat$Sample_ID, sample_ids)
  clades_raw <- stats::cutree(hc_cols, k = n_clades)
  
  # Assign clusters based on dendrogram order: leftmost = Cluster 1, rightmost = Cluster 2
  # Get the order of samples from the dendrogram
  ordered_samples <- names(clades_raw)[hc_cols$order]
  
  # Find which raw cluster appears first (leftmost) in the dendrogram order
  first_cluster_raw <- clades_raw[ordered_samples[1]]
  
  # Assign final cluster numbers: leftmost cluster becomes Cluster 1, rightmost becomes Cluster 2
  clades <- ifelse(clades_raw == first_cluster_raw, 1, 2)
  
  ids_ordered_clean <- names(clades)[hc_cols$order]
  ids_ordered_orig <- unname(id_map[ids_ordered_clean])

  cluster_df <- tibble::tibble(
    Sample_ID = ids_ordered_orig,
    Cluster    = unname(clades[ids_ordered_clean])
  )

  cluster_lists <- lapply(seq_len(n_clades), function(i) {
    cluster_df %>%
      dplyr::filter(Cluster == i) %>%
      dplyr::pull(Sample_ID)
  })
  names(cluster_lists) <- paste0("cluster", seq_len(n_clades), "_ids")

  # ---- Column annotation (aligned to columns of M) ----
  # Build annotation in REVERSE order since pheatmap displays first column at bottom
  # Start with variant annotation (will be on bottom - first column)
  ann_col <- data.frame(Clinical_PGD = dat$Clinical_PGD, row.names = sample_ids)
  
  # Add cluster annotation last (will be on top - last column)
  # Create a mapping from sample_ids to final cluster assignments
  cluster_mapping <- setNames(cluster_df$Cluster, cluster_df$Sample_ID)
  original_ids <- unname(id_map[sample_ids])  # Convert sample_ids back to original Sample_IDs
  cluster_labels <- paste0("Cluster ", cluster_mapping[original_ids])
  names(cluster_labels) <- sample_ids
  
  # Use standard factor levels - legend order controlled by color order
  ann_col$Clinical_PGD <- dat$Clinical_PGD

  # Reorder rows of ann_col to match M's columns
  ann_col <- ann_col[colnames(M), , drop = FALSE]

  # ---- Annotation color lists ----
  # Build colors in same order as ann_col data frame (pheatmap displays in reverse)
  ann_colors <- list()
  
  # Add variant colors first (matches first column - will display at bottom)
  ann_colors$Clinical_PGD <- c(
      "Y" = "#94001E",
      "N" = "#03507D"
    )
  # ---- Heatmap (for screen) ----
  heatmap_plot <- pheatmap::pheatmap(
    M,
    scale = "row",
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(255),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    fontsize = 10,
    na_col = "#DDDDDD",
    legend_labels = "Z-Score"
  )

  # Create heatmap plot object for patchwork
  heatmap_plot <- pheatmap::pheatmap(
    M,
    scale = "row",
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(255),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    fontsize = 8,
    na_col = "#DDDDDD",
    silent = TRUE,  # Prevents auto-display
    legend_labels = "Z-Score"
  )

  # Create heatmap plot to save
  png(
      filename = paste0(file_path, file_name, ".png"),
      width = 4,
      height = 4,
      units = "in",
      res = 600
    )
  print(
    pheatmap::pheatmap(
    M,
    scale = "row",
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(255),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    fontsize = 8,
    na_col = "#DDDDDD",
    legend_labels = "Z-Score"
  ))
  dev.off()

  list(
    M = M,
    Mz = Mz,
    hc_cols = hc_cols,
    sample_ids = sample_ids,
    group_pgd = group_pgd,
    group_time = group_time,
    feature_selector = feature_selector,
    top_features = top_features,
    ann_col = ann_col,
    ann_colors = ann_colors,
    cluster_df = cluster_df,  # Updated from clade_df
    clusters = clades,        # Updated from clades
    cluster_lists = cluster_lists,  # Updated from clade_lists
    heatmap_plot = heatmap_plot  # Plot object for patchwork
  )
}
