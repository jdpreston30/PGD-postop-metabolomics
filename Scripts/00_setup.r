#* 0: Dependencies and setting seeds
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

