relabel_pgd <- function(plot) {
  # Extract and modify existing scales to preserve colors
  plot_scales <- plot$scales$scales
  
  for (i in seq_along(plot_scales)) {
    scale <- plot_scales[[i]]
    
    # Check if this scale controls color or fill aesthetics
    if (any(c("colour", "fill") %in% scale$aesthetics)) {
      # Create label transformation function
      original_labels <- scale$get_labels
      
      scale$get_labels <- function(breaks) {
        labels <- if (!is.null(original_labels)) {
          original_labels(breaks)
        } else {
          breaks
        }
        gsub("Severe PGD", "sPGD", gsub("No Severe PGD", "No sPGD", labels))
      }
      
      plot_scales[[i]] <- scale
    }
  }
  
  plot$scales$scales <- plot_scales
  return(plot)
}
