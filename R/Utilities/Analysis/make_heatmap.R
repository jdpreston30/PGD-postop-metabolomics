#' Create heatmap with optional feature selection (ANOVA / variance / MAD)
#' and optional T_stage annotation. Returns plot object for patchwork.
#'
#' @param data                 Data frame with Sample_ID, Variant, optional T_stage, then feature columns
#' @param variant_colors       Named color vector for Variant (names = levels)
#' @param top_features         NULL (default: show all features). If numeric >0, keep top N by `feature_selector`.
#' @param feature_selector     One of c("none","anova","variance","mad"). Default "none".
#' @param variant_levels       Factor order for Variant (default c("PTC","FV-PTC","FTC"))
#' @param n_clades             Number of clades to extract from sample HCA (default 2)
#' @param annotate_t_stage     Logical; if TRUE and T_stage exists, add as column annotation
#' @param T_stage_colors       Optional named color vector for T_stage; include "T1","T2","T3","T4". "Unknown" added if missing.
#' @param cluster_colors       Named color vector for cluster annotation (e.g., c("Cluster 1" = "#color1", "Cluster 2" = "#color2"))
#' @return List with plot object, M, Mz, hc_cols, clade_df, clade_lists, ann_col, ann_colors, etc.
#' @export
make_heatmap <- function(
    data,
    variant_colors = c("PTC" = "#DF8D0A", "FV-PTC" = "#23744E", "FTC" = "#194992"),
    top_features = NULL,
    feature_selector = c("none", "anova", "variance", "mad", "ttest"),
    variant_levels = c("PTC", "FV-PTC", "FTC"),
    n_clades = 2,
    annotate_t_stage = FALSE,
    T_stage_colors = c("T1-T2" = "#dfba37", "T3-T4" = "#72061c"),
    p_anova = c("pgd", "time", "interaction")) {
  feature_selector <- match.arg(feature_selector)

  # ---- Checks ----
  stopifnot(all(c("Sample_ID", "Clinical_PGD") %in% names(data)))
  has_t <- "T_stage" %in% names(data)
  has_time <- "Time" %in% names(data)

  # Keep ID/Variant/(optional)T_stage up front
  dat <- dplyr::select(
    data,
    dplyr::any_of(c("Sample_ID", "Clinical_PGD", if (has_t) "T_stage")),
    dplyr::everything()
  )

  # Coerce Variant to factor in desired order
  dat$Clinical_PGD <- factor(dat$Clinical_PGD, levels = unique(dat$Clinical_PGD))

  # Build matrix: samples x features (drop ID/Variant/T_stage)
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
  ann_col <- data.frame(Variant = dat$Clinical_PGD, row.names = sample_ids)

  # Add T-stage annotation if requested (will be in middle - second column)
  if (annotate_t_stage && has_t) {
    # Turn T_stage into character - no need to handle NA since dataset is complete
    t_raw <- as.character(dat$T_stage)

    # Handle both individual T stages (T1, T2, T3, T4) and binned stages (T1-T2, T3-T4)
    # If values already start with "T" (e.g., "T3" or "T1-T2"), keep them;
    # if they look numeric (e.g., "3"), prepend "T".
    needs_T <- grepl("^[0-9]+$", t_raw)
    t_clean <- ifelse(needs_T, paste0("T", t_raw), t_raw)

    # Factor with canonical order (no Unknown needed since dataset is complete)
    # First check what format we have in the data
    unique_vals <- unique(t_clean)
    if (any(grepl("-", unique_vals))) {
      # We have binned values like "T1-T2", "T3-T4"
      t_lvls <- c("T1-T2", "T3-T4")
    } else {
      # We have individual values like "T1", "T2", "T3", "T4"
      t_lvls <- c("T1", "T2", "T3", "T4")
    }
    
    ann_col$`T Stage` <- factor(t_clean, levels = t_lvls)  # Use "T Stage" as column name for display
  }
  
  # Add cluster annotation last (will be on top - last column)
  # Create a mapping from sample_ids to final cluster assignments
  cluster_mapping <- setNames(cluster_df$Cluster, cluster_df$Sample_ID)
  original_ids <- unname(id_map[sample_ids])  # Convert sample_ids back to original Sample_IDs
  cluster_labels <- paste0("Cluster ", cluster_mapping[original_ids])
  names(cluster_labels) <- sample_ids
  
  # Use standard factor levels - legend order controlled by color order
  ann_col$Cluster <- factor(cluster_labels[sample_ids], levels = c("Cluster 1", "Cluster 2"))
  ann_col$Clinical_PGD <- dat$Clinical_PGD

  # Reorder rows of ann_col to match M's columns
  ann_col <- ann_col[colnames(M), , drop = FALSE]

  # ---- Annotation color lists ----
  # Build colors in same order as ann_col data frame (pheatmap displays in reverse)
  ann_colors <- list()
  
  # Add variant colors first (matches first column - will display at bottom)
  ann_colors$Clinical_PGD <- cluster_colors

  if (annotate_t_stage && has_t) {
    # Seed T-stage colors if not provided
    if (is.null(T_stage_colors)) T_stage_colors <- c()
    
    
    # Get the levels that actually appear in the data
    needed <- levels(ann_col$`T Stage`)  # Use "T Stage" column name
    
    # Find any missing colors and assign defaults
    missing <- setdiff(needed, names(T_stage_colors))
    if (length(missing)) {
      # Check if we have binned values and provide appropriate defaults
      if (any(grepl("-", missing))) {
        # Default colors for binned values
        default_binned <- c("T1-T2" = "#4575b4", "T3-T4" = "#d73027")
        fill_cols <- default_binned[missing]
        # For any still missing, use grays
        still_missing <- missing[is.na(fill_cols)]
        if (length(still_missing)) {
          gray_cols <- grDevices::gray.colors(length(still_missing), start = 0.3, end = 0.7)
          names(gray_cols) <- still_missing
          fill_cols[still_missing] <- gray_cols
        }
      } else {
        # Default for individual T stages
        fill_cols <- grDevices::gray.colors(length(missing), start = 0.3, end = 0.7)
        names(fill_cols) <- missing
      }
      T_stage_colors <- c(T_stage_colors, fill_cols)
    }
    
    ann_colors$`T Stage` <- T_stage_colors[needed]  # Use "T Stage" as key
  }
  
  # Add cluster colors last (matches last column - will display at top)
  # Order colors so Cluster 1 appears on top in legend
  if (!is.null(cluster_colors)) {
    # Keep the standard order of cluster colors
    ann_colors$Cluster <- cluster_colors[c("Cluster 1", "Cluster 2")]
  } else {
    # Default cluster colors if not provided
    default_cluster_colors <- c("Cluster 1" = "#94001E", "Cluster 2" = "#03507D")
    ann_colors$Cluster <- default_cluster_colors
  }
  
  # Add variant colors last (bottom annotation)
  ann_colors$Clinical_PGD <- cluster_colors

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
