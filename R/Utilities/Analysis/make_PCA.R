#' Create PCA plot with ellipses
#'
#' @param data Data frame with Patient_ID, Clinical_PGD, and feature columns
#' @param plot_title Optional title for the plot (default: "")
#' @param ellipse_colors Named vector of colors for each Clinical_PGD/group
#' @param point_size Size of the points (default: 3 for standalone, 0.5 for multi-panel)
#' @param comp_x Which principal component to plot on x-axis (default: 1)
#' @param comp_y Which principal component to plot on y-axis (default: 2)
#' @return List containing the plot, PCA object, scores, scores_df, and explained variance
#' @export
make_PCA <- function(data, plot_title = "", 
                        ellipse_colors = c("Y" = "#94001E", "N" = "#03507D"), point_size = 3, comp_x = 1, comp_y = 2, group_var, method) {
      #_Data preparation
      df <- as.data.frame(data)
      cls_col <- if (group_var %in% names(df)) group_var else names(df)[2]
      X <- df[, -c(1, 2, 3), drop = FALSE]
      #_Coerce to numeric safely
      X[] <- lapply(X, function(v) suppressWarnings(as.numeric(v)))
      #_Handle NAs with median imputation
      if (anyNA(X)) {
        X[] <- lapply(X, function(v) {
          v[is.na(v)] <- stats::median(v, na.rm = TRUE)
          v
        })
      }
      Y <- factor(df[[cls_col]])
      if (method == "PCA"){
      #_Perform PCA
      pca <- stats::prcomp(X, center = TRUE, scale. = TRUE)
      
      #_Validate component indices
      max_comp <- min(ncol(X), nrow(X) - 1)
      if (comp_x > max_comp || comp_y > max_comp) {
        stop(paste("Requested components exceed available components. Max components:", max_comp))
      }
      
      scores <- pca$x[, c(comp_x, comp_y), drop = FALSE]
      explained <- round((pca$sdev^2 / sum(pca$sdev^2))[c(comp_x, comp_y)] * 100)
      model_obj <- pca
      } else if (method == "PLSDA") {
        if (!requireNamespace("mixOmics", quietly = TRUE)) {
          stop("Please install the 'mixOmics' package to use PLS-DA")
        }
        
        plsda_obj <- mixOmics::plsda(X, Y, ncomp = max(comp_x, comp_y))
        scores <- plsda_obj$variates$X[, c(comp_x, comp_y), drop = FALSE]
        explained <- round(plsda_obj$prop_expl_var$X[c(comp_x, comp_y)] * 100)
        
        model_obj <- plsda_obj
      }

      #_Prepare plot data
      scores_df <- data.frame(
        Comp1 = scores[, 1],
        Comp2 = scores[, 2],
        Class = Y
      )
      
      # Identify NA values and create separate datasets
      na_mask <- is.na(scores_df$Class)
      scores_df_complete <- scores_df[!na_mask, , drop = FALSE]
      scores_df_na <- scores_df[na_mask, , drop = FALSE]
      
      #_Create PCA plot
      axis_prefix <- ifelse(method == "PCA", "PC", "LV")
      pca_plot <- ggplot2::ggplot() +
        # Plot complete cases with colors and ellipses
        {if(nrow(scores_df_complete) > 0) {
          list(
            ggplot2::geom_point(
              data = scores_df_complete,
              ggplot2::aes(x = Comp1, y = Comp2, group = Class, fill = Class),
              size = point_size, shape = 21
            ),
            ggplot2::stat_ellipse(
              data = scores_df_complete,
              ggplot2::aes(x = Comp1, y = Comp2, group = Class, fill = Class),
              geom = "polygon", alpha = 0.3, color = NA
            )
          )
        }} +
        # Plot NA values as open circles without color
        {if(nrow(scores_df_na) > 0) {
          ggplot2::geom_point(
            data = scores_df_na,
            ggplot2::aes(x = Comp1, y = Comp2),
            size = point_size, shape = 1, color = "black", fill = NA
          )
        }} +
        ggplot2::scale_color_manual(values = ellipse_colors, name = "Clinical PGD", drop = TRUE, na.translate = FALSE) +
        ggplot2::scale_fill_manual(values = ellipse_colors, name = "Clinical PGD", drop = TRUE, na.translate = FALSE) +
        ggplot2::theme_minimal(base_family = "Arial") +
        ggplot2::labs(
          x = paste0(axis_prefix, comp_x, " (", explained[1], "%)"),
          y = paste0(axis_prefix, comp_y, " (", explained[2], "%)")
        ) +
        ggplot2::theme(
          axis.title = ggplot2::element_text(size = 25, face = "bold"),
          axis.text = ggplot2::element_text(size = 22, face = "bold", color = "black"),
          panel.grid.major = ggplot2::element_line(color = "gray80", linewidth = 0.8, linetype = "solid"),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 3.2),
          panel.background = ggplot2::element_blank(),
          aspect.ratio = 1,
          plot.margin = grid::unit(c(2, 8, 8, 8), "pt"),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.box.margin = ggplot2::margin(0, 0, 0, 0),
          legend.margin = ggplot2::margin(t = 0, r = 0, b = -3, l = 0),

          # 🔽 controls symbol size
          legend.key.width = grid::unit(0.35, "cm"),
          legend.key.height = grid::unit(0.35, "cm"),
          legend.key.size = grid::unit(0.35, "cm"),

          # 🔽 distance between symbol and text
          legend.spacing.x = grid::unit(0, "cm"),

          # 🔽 fine-tunes text alignment (0 = left, 1 = right; <0 moves text closer to key)
          legend.text.align = -10,
          legend.text = element_text(size = 12),
          legend.title = ggplot2::element_text(size = 14),
          axis.title.y = ggplot2::element_text(
            margin = ggplot2::margin(r = 0), hjust = 0.5
          )
        ) + 
        ggplot2::guides(
          color = guide_legend(title = "Clinical PGD")
        )
      
      #_Return useful objects for further analysis
      return(list(
        plot = pca_plot,
        model = model_obj,
        scores = scores,
        scores_df = scores_df,
        explained = explained,
        Y = Y
      ))
}
