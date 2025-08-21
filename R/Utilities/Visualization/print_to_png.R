#' Print plot to PNG with auto-refresh for macOS Preview
#'
#' @param plot The plot object to print
#' @param filename Name of the PNG file (with or without .png extension)
#' @param width Width in inches (default: 8.5)
#' @param height Height in inches (default: 11)
#' @param dpi Resolution in DPI (default: 300 for high quality)
#' @param output_dir Directory to save the PNG (default: "Outputs")
#' @param auto_open Whether to automatically open in Preview on first run (default: TRUE)
#' @return Invisible path to the created PNG file
#' @export
print_to_png <- function(plot, filename, width = 8.5, height = 11, dpi = 300, 
                         output_dir = "Outputs", auto_open = TRUE) {
  # Ensure filename has .png extension
  if (!grepl("\\.png$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".png")
  }
  
  # Create full path
  filepath <- file.path(output_dir, filename)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Check if file already exists (for auto-open logic)
  file_exists <- file.exists(filepath)
  
  # Save the plot as PNG
  ggplot2::ggsave(
    filename = filepath,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    device = "png"
  )
  
  # Auto-open in Preview only on first run (or if specified)
  if (auto_open && !file_exists) {
    system(paste("open", shQuote(filepath)))
    cat("PNG saved and opened in Preview:", filepath, "\n")
    cat("Preview will auto-refresh when you re-run this function!\n")
  } else {
    cat("PNG updated:", filepath, "\n")
  }
  
  # Return path invisibly
  invisible(filepath)
}
