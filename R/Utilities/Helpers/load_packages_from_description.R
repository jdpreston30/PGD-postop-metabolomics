#' Load all packages specified in DESCRIPTION file
#'
#' This function reads the DESCRIPTION file and loads all packages listed in
#' the Imports and Bioconductor fields. It provides progress feedback and
#' handles errors gracefully.
#'
#' @return Invisible NULL
#' @export
load_packages_from_description <- function() {
  desc_file <- "DESCRIPTION"
  
  if (!file.exists(desc_file)) {
    stop("DESCRIPTION file not found. Please ensure you're in the project root directory.")
  }
  
  # Read DESCRIPTION file
  desc_lines <- readLines(desc_file)
  
  # Extract Imports section
  imports_start <- which(grepl("^Imports:", desc_lines))
  if (length(imports_start) == 0) {
    stop("No Imports section found in DESCRIPTION file.")
  }
  
  # Find where Imports section ends
  next_field <- which(grepl("^[A-Z]", desc_lines[(imports_start + 1):length(desc_lines)]))
  if (length(next_field) > 0) {
    imports_end <- imports_start + next_field[1] - 1
  } else {
    imports_end <- length(desc_lines)
  }
  
  # Extract package names from Imports
  imports_lines <- desc_lines[imports_start:imports_end]
  imports_text <- paste(imports_lines, collapse = " ")
  imports_text <- gsub("^Imports:", "", imports_text)
  imports_text <- gsub(",", " ", imports_text)
  imports_text <- gsub("\\(.*?\\)", "", imports_text)  # Remove version specifications
  packages <- trimws(unlist(strsplit(imports_text, " ")))
  packages <- packages[packages != ""]
  
  # Extract Bioconductor packages
  bioc_start <- which(grepl("^Bioconductor:", desc_lines))
  bioc_packages <- character(0)
  if (length(bioc_start) > 0) {
    next_field <- which(grepl("^[A-Z]", desc_lines[(bioc_start + 1):length(desc_lines)]))
    if (length(next_field) > 0) {
      bioc_end <- bioc_start + next_field[1] - 1
    } else {
      bioc_end <- length(desc_lines)
    }
    bioc_lines <- desc_lines[bioc_start:bioc_end]
    bioc_text <- paste(bioc_lines, collapse = " ")
    bioc_text <- gsub("^Bioconductor:", "", bioc_text)
    bioc_text <- gsub(",", " ", bioc_text)
    bioc_packages <- trimws(unlist(strsplit(bioc_text, " ")))
    bioc_packages <- bioc_packages[bioc_packages != ""]
  }
  
  # Combine all packages
  all_packages <- c(packages, bioc_packages)
  all_packages <- unique(all_packages)
  
  cat("📚 Loading", length(all_packages), "packages from DESCRIPTION...\n")
  
  # Load packages
  for (pkg in all_packages) {
    suppressPackageStartupMessages({
      library(pkg, character.only = TRUE)
    })
  }
  
  cat("✅ All packages loaded successfully!\n")
  invisible(NULL)
}
