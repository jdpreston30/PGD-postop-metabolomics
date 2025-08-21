#* 0: Dependencies and setting seeds
  #+ 0.1: Dependencies
    #- 0.1.1: Install all packages
      install.packages(c("dplyr", "tidyr", "readr", "stringr"))
      # devtools::install_github("jdpreston30/TernTablesR")
    #- 0.1.2: Load libraries
      library(dplyr)
      library(tidyr)
      library(readr)
      library(stringr)
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
