#+ 0a.1: Read required packages from DESCRIPTION file
desc_file <- "DESCRIPTION"
if (!file.exists(desc_file)) {
  stop("DESCRIPTION file not found. Please ensure you're in the project root directory.")
}
#- 0a.1.1: Read DESCRIPTION file
desc_lines <- readLines(desc_file)
#- 0a.1.2: Extract Imports section
imports_start <- which(grepl("^Imports:", desc_lines))
if (length(imports_start) == 0) {
  stop("No Imports section found in DESCRIPTION file.")
}
#- 0a.1.3: Find where Imports section ends (next field or end of file)
next_field <- which(grepl("^[A-Z]", desc_lines[(imports_start + 1):length(desc_lines)])) # nolint
if (length(next_field) > 0) {
  imports_end <- imports_start + next_field[1] - 1
} else {
  imports_end <- length(desc_lines)
}
#- 0a.1.4: Extract package names
imports_lines <- desc_lines[imports_start:imports_end]
imports_text <- paste(imports_lines, collapse = " ")
#- 0a.1.5: Clean up and extract package names
imports_text <- gsub("Imports:", "", imports_text)
imports_text <- gsub("\\s+", " ", imports_text)
packages <- strsplit(imports_text, ",")[[1]]
required_packages <- trimws(packages)
required_packages <- required_packages[required_packages != ""]
cat("📋 Found", length(required_packages), "packages in DESCRIPTION file\n")
#+ 0a.2: Install missing packages
#- 0a.2.1: Check for missing packages
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  cat("📦 Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  #- 0a.2.2: Install CRAN packages
  cran_packages <- setdiff(missing_packages, c("mixOmics", "KEGGREST"))
  if (length(cran_packages) > 0) {
    install.packages(cran_packages, repos = "https://cran.rstudio.com/")
  }
  #- 0a.2.3: Install Bioconductor packages
  bioc_packages <- intersect(missing_packages, c("mixOmics", "KEGGREST"))
  if (length(bioc_packages) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(bioc_packages)
  }
}
#- 0a.2.4: Install GitHub packages
if (!require(TernTablesR, quietly = TRUE)) {
  remotes::install_github("jdpreston30/TernTablesR")
  library(TernTablesR)
}

#- 0a.2.5: Install MetaboAnalystR and dependencies (if needed)
if (!require(MetaboAnalystR, quietly = TRUE)) {
  cat("🧪 Installing MetaboAnalystR and dependencies...\n")

  # Install BiocManager if not available
  if (!require(BiocManager, quietly = TRUE)) {
    install.packages("BiocManager")
    library(BiocManager)
  }

  # Install Bioconductor dependencies
  BiocManager::install(c(
    "RBGL",
    "Rgraphviz",
    "fgsea",
    "globaltest",
    "GlobalAncova",
    "Rsamtools",
    "edgeR",
    "siggenes",
    "BiocParallel",
    "MSnbase",
    "xcms",
    "CAMERA",
    "multtest"
  ), update = FALSE)

  # Install CRAN dependencies
  install.packages(c(
    "igraph",
    "RColorBrewer",
    "xtable",
    "som",
    "RJSONIO",
    "gplots",
    "e1071",
    "caTools",
    "randomForest",
    "Cairo",
    "pls",
    "caret",
    "lattice",
    "scatterplot3d",
    "devtools",
    "qqconf",
    "metap"
  ))

  # Install MetaboAnalystR from GitHub
  remotes::install_github("xia-lab/MetaboAnalystR", dependencies = FALSE)
  library(MetaboAnalystR)
  cat("✅ MetaboAnalystR installation complete!\n")
}

#+ 0a.3: Load all required packages
cat("📚 Loading required packages...\n")
invisible(sapply(required_packages, library, character.only = TRUE, quietly = TRUE))
#+ 0a.4: Environment setup complete
cat("✅ Environment setup complete! All required packages loaded.\n")

  #+ 0.1: Dependencies
    #- 0.1.1: Install all packages
      {if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(c(
        "dplyr", "tidyr", "readr", "stringr", "purrr", "here",
        "ggplot2", "tibble", "forcats", "scales", "limma"
      ))}
    #- 0.1.2: Load libraries
        library(dplyr)
        library(tidyr)
        library(readr)
        library(stringr)
        library(limma)
        library(purrr)
        library(here)
        library(ggplot2)
        library(tibble)
        library(forcats)
        library(scales)
        library(limma)
  #+ 0.2: Call all utility and modeling functions
    purrr::walk(
      list.files(
        here::here("R", "Utilities"),
        pattern = "\\.[rR]$",
        full.names = TRUE,
        recursive = TRUE
      ),
      source
    )
  #+ 0.3: Set conflicts
    conflicts_prefer(ggplot2::margin)
  #+ 0.4: Set seed for reproducibility
    set.seed(2025)
  #+ 0.5: Set Vroom connection size
    Sys.setenv("VROOM_CONNECTION_SIZE" = 2^22)
  #+ Set Vroom connection size
    Sys.setenv("VROOM_CONNECTION_SIZE" = 2^22)

