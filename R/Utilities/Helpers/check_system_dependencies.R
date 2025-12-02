#' Check for required system dependencies
#'
#' This function verifies that necessary system tools and libraries are
#' available for the analysis pipeline.
#'
#' @return Invisible NULL
#' @export
check_system_dependencies <- function() {
  cat("🔍 Checking system dependencies...\n")
  
  # Check for pandoc (required for rmarkdown)
  pandoc_available <- nzchar(Sys.which("pandoc"))
  if (!pandoc_available) {
    warning("⚠️  pandoc not found. Install from https://pandoc.org/ or use RStudio")
  } else {
    cat("✅ pandoc found\n")
  }
  
  # Check for ImageMagick (required for magick package)
  magick_available <- nzchar(Sys.which("convert"))
  if (!magick_available) {
    warning("⚠️  ImageMagick not found. Install with: brew install imagemagick")
  } else {
    cat("✅ ImageMagick found\n")
  }
  
  # Check R version
  r_version <- as.numeric_version(paste(R.version$major, R.version$minor, sep = "."))
  if (r_version < "4.5.1") {
    warning("⚠️  R version ", r_version, " detected. Recommended: >= 4.5.1")
  } else {
    cat("✅ R version", as.character(r_version), "\n")
  }
  
  cat("✅ System dependency check complete\n")
  invisible(NULL)
}
