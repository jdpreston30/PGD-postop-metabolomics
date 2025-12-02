#* 0a: Environment Setup
#' Package management is handled by renv for reproducibility.
#' All required packages and their exact versions are specified in renv.lock.
#' 
#' This script automatically checks for missing packages and runs renv::restore()
#' if needed, making the pipeline self-sufficient for first-time setup.
#' 
#' The .Rprofile automatically activates the renv environment when R starts.

# Suppress package update checks
Sys.setenv(R_BIOC_VERSION_CHECK_SKIP = "TRUE")

#+ 0a.1: Verify renv is active
cat("📦 Package environment managed by renv\n")
if (!("renv" %in% loadedNamespaces())) {
  warning("⚠️  renv is not active. Attempting to activate...")
  if (file.exists("renv/activate.R")) {
    source("renv/activate.R")
  } else {
    warning("⚠️  renv not initialized. Run renv::init() to set up.\n")
    cat("   For now, continuing without renv...\n")
  }
}

#+ 0a.2: Check if packages need to be installed
core_packages <- c("dplyr", "tidyr", "readr", "stringr", "purrr", "here", 
                   "ggplot2", "tibble", "forcats", "scales", "limma", "conflicted")
missing_core <- core_packages[!sapply(core_packages, requireNamespace, quietly = TRUE)]

if (length(missing_core) > 0) {
  cat("⚠️  Core packages missing:", paste(missing_core, collapse = ", "), "\n")
  
  # Check if renv.lock exists
  if (file.exists("renv.lock")) {
    cat("🔄 Running renv::restore() to install packages...\n")
    cat("   (This may take 10-20 minutes on first run)\n\n")
    
    # Run renv::restore() automatically
    tryCatch({
      renv::restore(prompt = FALSE)  # No prompt, automatic yes
      cat("\n✅ Package installation complete!\n")
    }, error = function(e) {
      stop("❌ Failed to restore packages: ", e$message, 
           "\n   Please run renv::restore() manually and check for errors.")
    })
    
    # Verify installation succeeded
    still_missing <- core_packages[!sapply(core_packages, requireNamespace, quietly = TRUE)]
    if (length(still_missing) > 0) {
      stop("❌ Packages still missing after restore: ", paste(still_missing, collapse = ", "),
           "\n   Please check renv::status() for details.")
    }
  } else {
    cat("⚠️  renv.lock not found. Installing packages manually...\n")
    
    # Fallback: Install packages manually (preserving your original logic)
    #- 0a.2.1: Install CRAN packages
    cran_packages <- setdiff(missing_core, c("limma", "mixOmics"))
    if (length(cran_packages) > 0) {
      install.packages(cran_packages, repos = "https://cran.rstudio.com/")
    }
    
    #- 0a.2.2: Install Bioconductor packages
    bioc_packages <- intersect(missing_core, c("limma", "mixOmics"))
    if (length(bioc_packages) > 0) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(bioc_packages)
    }
    
    #- 0a.2.3: Install GitHub packages
    if (!requireNamespace("TernTablesR", quietly = TRUE)) {
      if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
      remotes::install_github("jdpreston30/TernTablesR")
    }
    
    #- 0a.2.4: Install MetaboAnalystR and dependencies (if needed)
    if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
      cat("🧪 Installing MetaboAnalystR and dependencies...\n")
      
      # Install BiocManager if not available
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      
      # Install Bioconductor dependencies for MetaboAnalystR
      BiocManager::install(c(
        "RBGL", "Rgraphviz", "fgsea", "globaltest", "GlobalAncova",
        "Rsamtools", "edgeR", "siggenes", "BiocParallel", "MSnbase",
        "xcms", "CAMERA", "multtest"
      ), update = FALSE)
      
      # Install CRAN dependencies for MetaboAnalystR
      install.packages(c(
        "igraph", "RColorBrewer", "xtable", "som", "RJSONIO", "gplots",
        "e1071", "caTools", "randomForest", "Cairo", "pls", "caret",
        "lattice", "scatterplot3d", "devtools", "qqconf", "metap"
      ))
      
      # Install MetaboAnalystR from GitHub
      remotes::install_github("xia-lab/MetaboAnalystR", dependencies = FALSE)
      cat("✅ MetaboAnalystR installation complete!\n")
    }
  }
} else {
  cat("✅ renv environment verified. All core packages available.\n")
}

#+ 0a.3: Load conflicted and set preferences BEFORE loading other packages
library(conflicted)

# Set conflict preferences to prevent warnings during package loading
# Based on common conflicts in this project's package ecosystem
conflicts_prefer(ggplot2::margin)       # ggplot2 vs randomForest
conflicts_prefer(dplyr::filter)         # dplyr vs stats
conflicts_prefer(dplyr::select)         # dplyr vs MASS, raster
conflicts_prefer(dplyr::summarize)      # dplyr vs various
conflicts_prefer(dplyr::summarise)      # dplyr vs various
conflicts_prefer(dplyr::first)          # dplyr vs data.table
conflicts_prefer(dplyr::mutate)         # dplyr vs various
conflicts_prefer(dplyr::arrange)        # dplyr vs various
conflicts_prefer(dplyr::count)          # dplyr vs various
conflicts_prefer(dplyr::rename)         # dplyr vs various
conflicts_prefer(dplyr::combine)        # dplyr vs randomForest, gridExtra
conflicts_prefer(purrr::map)            # purrr vs mixOmics
conflicts_prefer(stats::chisq.test)     # stats vs janitor
conflicts_prefer(stats::fisher.test)    # stats vs janitor
conflicts_prefer(jsonlite::fromJSON)    # jsonlite vs RJSONIO
conflicts_prefer(readxl::read_xlsx)     # readxl vs officer
conflicts_prefer(raster::intersect)     # raster vs base, dplyr, lubridate
conflicts_prefer(igraph::compose)       # igraph vs flextable, purrr
conflicts_prefer(flextable::align)      # flextable vs xtable
conflicts_prefer(base::setdiff)         # base vs various
conflicts_prefer(base::intersect)       # base vs various
conflicts_prefer(base::union)           # base vs various, igraph
conflicts_prefer(base::as.factor)       # base vs various
conflicts_prefer(base::unique)          # base vs various
conflicts_prefer(base::as.data.frame)   # base vs various, tibble
conflicts_prefer(dplyr::desc)
conflicts_prefer(purrr::compact)
conflicts_prefer(ggplot2::margin)
conflicts_prefer(scales::alpha)
#+ 0a.4: Load all packages from DESCRIPTION file
cat("📚 Loading required packages...\n")
source("R/Utilities/Helpers/load_packages_from_description.R")
load_packages_from_description()

#+ 0a.5: Load GitHub Packages explicitly for renv detection
library(TernTablesR)
library(MetaboAnalystR)

#+ 0a.6: Check system dependencies
source("R/Utilities/Helpers/check_system_dependencies.R")
check_system_dependencies()

#+ 0a.7: Call all utility and modeling functions
purrr::walk(
  list.files(
    here::here("R", "Utilities"),
    pattern = "\\.[rR]$",
    full.names = TRUE,
    recursive = TRUE
  ),
  source
)

#+ 0a.8: Set seed for reproducibility
set.seed(2025)

#+ 0a.9: Set Vroom connection size
Sys.setenv("VROOM_CONNECTION_SIZE" = 2^22)

cat("✅ Environment setup complete! All required packages loaded.\n")

