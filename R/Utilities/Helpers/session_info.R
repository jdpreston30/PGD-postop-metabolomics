#' Session Information
#' 
#' Captures the exact R and package versions used for this analysis.
#' This function generates and saves comprehensive session information
#' for reproducibility documentation.
#' 
#' @param session_file Path to save session info. Default: "session_info.txt" in project root
#' @param verbose Logical. Print session info to console? Default TRUE
#' @return Invisibly returns the session info output
#' @export
session_info <- function(session_file = here::here("session_info.txt"), verbose = TRUE) {
  
  if (verbose) {
    cat("=== R Session Information ===\n")
    cat("Analysis conducted on:", format(Sys.time()), "\n\n")
  }
  
  # Capture R version and platform info
  session_output <- capture.output(sessionInfo())
  
  if (verbose) {
    print(sessionInfo())
  }

  # Check and capture system dependency versions
  if (verbose) {
    cat("\n=== System Dependencies ===\n")
  }
  
  system_deps_info <- character()
  
  # Ghostscript
  gs_version <- tryCatch({
    system("gs --version", intern = TRUE)
  }, error = function(e) "Not found")
  if (verbose) cat("Ghostscript:", gs_version, "\n")
  system_deps_info <- c(system_deps_info, paste("Ghostscript:", gs_version))
  
  # ImageMagick
  im_version <- tryCatch({
    system("convert --version | head -n1", intern = TRUE)
  }, error = function(e) "Not found")
  if (verbose) cat("ImageMagick:", im_version, "\n")
  system_deps_info <- c(system_deps_info, paste("ImageMagick:", im_version))
  
  # Pandoc
  pandoc_version <- tryCatch({
    system("pandoc --version | head -n1", intern = TRUE)
  }, error = function(e) "Not found")
  if (verbose) cat("Pandoc:", pandoc_version, "\n")
  system_deps_info <- c(system_deps_info, paste("Pandoc:", pandoc_version))
  
  # LaTeX (pdflatex)
  latex_version <- tryCatch({
    system("pdflatex --version | head -n1", intern = TRUE)
  }, error = function(e) "Not found")
  if (verbose) cat("LaTeX:", latex_version, "\n")
  system_deps_info <- c(system_deps_info, paste("LaTeX:", latex_version))
  
  # Prepare complete session output for file
  complete_output <- c(
    paste("Analysis run on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    session_output,
    "",
    "=== System Dependencies ===",
    system_deps_info
  )
  
  # Save to file
  writeLines(complete_output, session_file)
  if (verbose) {
    cat("\n📄 Session info saved to:", session_file, "\n")
  }
  
  # Print package versions in a clean format
  if (verbose) {
    cat("\n=== Key Package Versions ===\n")
  }

  # Read packages from DESCRIPTION file
  desc_file <- here::here("DESCRIPTION")
  desc_lines <- readLines(desc_file)
  imports_start <- which(grepl("^Imports:", desc_lines))
  next_field <- which(grepl("^[A-Z]", desc_lines[(imports_start + 1):length(desc_lines)]))
  if (length(next_field) > 0) {
    imports_end <- imports_start + next_field[1] - 1
  } else {
    imports_end <- length(desc_lines)
  }
  imports_lines <- desc_lines[imports_start:imports_end]
  imports_text <- paste(imports_lines, collapse = " ")
  imports_text <- gsub("Imports:", "", imports_text)
  imports_text <- gsub("\\s+", " ", imports_text)
  packages <- strsplit(imports_text, ",")[[1]]
  required_packages <- trimws(packages)
  required_packages <- required_packages[required_packages != ""]
  
  if (verbose) {
    for (pkg in required_packages) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        version <- as.character(packageVersion(pkg))
        cat(sprintf("%-15s %s\n", pkg, version))
      }
    }
  }
  
  invisible(complete_output)
}
