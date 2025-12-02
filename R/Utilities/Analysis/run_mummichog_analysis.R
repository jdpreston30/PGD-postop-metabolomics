#' Run Mummichog Pathway Enrichment Analysis
#'
#' This function performs mummichog pathway enrichment analysis using MetaboAnalystR.
#' It takes tibble input, creates temporary CSV files, and runs the complete analysis
#' pipeline with configurable parameters.
#'
#' @param ttest_results A tibble containing t-test results with required columns
#' @param output_dir Character. Directory path where analysis results will be saved
#' @param database Character. Must be "hsa_kegg" or "hsa_mfn" for pathway database
#' @param instrumentOpt Numeric. Mass-spec instrument parameter (default: 5.0)
#' @param msModeOpt Character. Mass-spec mode (default: "mixed")
#' @param force_primary_ion Character. Primary ion filtering, "yes" or "no" (default: "yes")
#' @param rt_frac Numeric. RT fraction parameter (default: 0.02)
#' @param rt_tol Numeric. RT tolerance parameter (default: NA)
#'
#' @return A MetaboAnalyst mSet object containing analysis results
#'
#' @examples
#' \dontrun{
#' # MFN database analysis
#' result <- run_mummichog_analysis(
#'   ttest_results = my_ttest_data,
#'   output_dir = "outputs/analysis_mfn",
#'   database = "hsa_mfn"
#' )
#'
#' # KEGG database with custom parameters
#' result <- run_mummichog_analysis(
#'   ttest_results = my_ttest_data,
#'   output_dir = "outputs/analysis_kegg",
#'   database = "hsa_kegg",
#'   instrumentOpt = 3.0,
#'   msModeOpt = "positive"
#' )
#' }
#'
#' @export
run_mummichog_analysis <- function(
  ttest_results,
  output_dir,
  database, # Must be "hsa_kegg" or "hsa_mfn"
  instrumentOpt = 5.0, # Numeric - mass-spec instrument
  msModeOpt = "mixed", # Character - mass-spec mode
  force_primary_ion = "yes", # Character - primary ion filtering
  rt_frac = 0.02, # RT fraction
  rt_tol = NA # RT tolerance
) {
  # Validate required parameters
  if (missing(output_dir) || is.null(output_dir) || output_dir == "") {
    stop("output_dir parameter is required and cannot be empty")
  }
  if (missing(ttest_results) || is.null(ttest_results)) {
    stop("ttest_results parameter is required and cannot be NULL")
  }
  if (missing(database) || is.null(database)) {
    stop("database parameter is required and cannot be NULL")
  }
  
  library(MetaboAnalystR)

  # FIRST: Create output directory using absolute path
  cat("Original working directory:", getwd(), "\n")
  
  # Convert to absolute path if relative
  if (!startsWith(output_dir, "/") && !startsWith(output_dir, "~")) {
    output_dir <- file.path(getwd(), output_dir)
  }
  
  cat("Trying to create output directory:", output_dir, "\n")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  } else {
    cat("Output directory already exists:", output_dir, "\n")
  }

  # Normalize path AFTER creation so it exists
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  cat("✓ Directory confirmed to exist (normalized path):", output_dir, "\n")

  original_wd <- getwd()
  # Change to output directory so MetaboAnalystR writes files there
  setwd(output_dir)
  cat("Changed working directory to:", getwd(), "\n")
  cat("Output directory (absolute path):", output_dir, "\n")

  # Set up database caching to avoid re-downloading
  # Get the project root directory
  project_root <- original_wd # Use original working directory as project root
  while (!file.exists(file.path(project_root, "run.R")) && dirname(project_root) != project_root) {
    project_root <- dirname(project_root)
  }
  cache_dir <- file.path(project_root, "Databases")

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    cat("Created database cache directory:", cache_dir, "\n")
  }

  # Set MetaboAnalyst to use cache directory
  if (exists(".on.public.web")) {
    .on.public.web <- FALSE # Ensure local mode
  }

  # Validate database
  if (!database %in% c("hsa_kegg", "hsa_mfn")) {
    stop("Database must be 'hsa_kegg' or 'hsa_mfn'")
  }

  # Create temp CSV file in current directory (which is now output_dir)
  temp_csv <- "temp_input.csv"
  readr::write_csv(ttest_results, temp_csv)

  # Suppress common non-critical MetaboAnalystR warnings
  suppressWarnings({
    # EXACT WORKING SEQUENCE with configurable parameters
    mSet <- InitDataObjects("mass_all", "mummichog", FALSE, 150)

    # Check if mSet was properly initialized
    if (is.null(mSet) || !is.list(mSet)) {
      stop("Failed to initialize mSet object")
    }
    
    cat("✓ InitDataObjects completed\n")

    mSet <- SetPeakFormat(mSet, "rmp")
    if (!is.list(mSet)) stop("mSet corrupted after SetPeakFormat")
    cat("✓ SetPeakFormat completed\n")

    mSet <- UpdateInstrumentParameters(mSet, instrumentOpt, msModeOpt, force_primary_ion, rt_frac)
    if (!is.list(mSet)) stop("mSet corrupted after UpdateInstrumentParameters")
    cat("✓ UpdateInstrumentParameters completed\n")

    mSet <- Read.PeakListData(mSet, temp_csv)
    if (!is.list(mSet)) stop("mSet corrupted after Read.PeakListData")
    cat("✓ Read.PeakListData completed\n")

    mSet <- SanityCheckMummichogData(mSet)
    if (!is.list(mSet)) stop("mSet corrupted after SanityCheckMummichogData")
    cat("✓ SanityCheckMummichogData completed\n")

    mSet <- SetPeakEnrichMethod(mSet, "mum", "v2")
    if (!is.list(mSet)) stop("mSet corrupted after SetPeakEnrichMethod")
    cat("✓ SetPeakEnrichMethod completed\n")

    # Dynamic p-value threshold selection (top 10% of peaks)
    # Custom implementation of SetMummichogPvalFromPercent (fixed version)
    fraction <- 0.1 # Top 10% of peaks
    peakFormat <- mSet$paramSet$peakFormat

    if (peakFormat %in% c("rmp", "rmt")) {
      maxp <- 0
    } else {
      pvals <- c(0.25, 0.2, 0.15, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)
      ndat <- mSet$dataSet$mummi.proc
      n <- floor(fraction * length(ndat[, "p.value"]))
      actual_cutoff <- ndat[n + 1, 1] # This is the precise value you remember!
      if (!any(pvals <= actual_cutoff)) {
        maxp <- 0.00001
      } else {
        maxp <- max(pvals[pvals <= actual_cutoff])
      }
      # Store both values for reporting
      mSet$dataSet$actual_cutoff <- actual_cutoff
    }

    mSet$dataSet$cutoff <- maxp
    mSet <- SetMummichogPval(mSet, maxp) # Fixed: pass mSet instead of NA
    cat("✓ Dynamic p-value threshold calculated:", maxp, "(top 10% of peaks)\n")
    cat("✓ Actual 10th percentile p-value:", actual_cutoff, "\n")

    if (!is.list(mSet)) stop("mSet corrupted after SetMummichogPval")

    mSet <- PerformPSEA(mSet, database, "current", 3, 100)
    if (!is.list(mSet)) stop("mSet corrupted after PerformPSEA")
    cat("✓ PerformPSEA completed\n")

    # Debug: Check what pathway statistics are available
    cat("Available pathway fields:\n")
    if (!is.null(mSet$path.pval)) cat("  - path.pval: ", length(mSet$path.pval), " values\n")
    if (!is.null(mSet$path.fdr)) cat("  - path.fdr: ", length(mSet$path.fdr), " values\n")
    if (!is.null(mSet$pathway.fdr)) cat("  - pathway.fdr: ", length(mSet$pathway.fdr), " values\n")
    if (!is.null(mSet$dataSet$pathway.fdr)) cat("  - dataSet$pathway.fdr: ", length(mSet$dataSet$pathway.fdr), " values\n")

    mSet <- PlotPeaks2Paths(mSet, "peaks_to_paths_0_", "png", 150, width = NA)
    if (!is.list(mSet)) stop("mSet corrupted after PlotPeaks2Paths")
    cat("✓ PlotPeaks2Paths completed\n")

    # Step 10: Prepare Enrichment Network (equivalent to PrepareEnrichNet)
    build_network_path <- file.path(original_wd, "R/Utilities/Analysis/prepare_custom_enrichnet.R")
    if (file.exists(build_network_path)) {
      source(build_network_path)
      network_name <- paste0("enrichNet_", gsub("hsa_", "", database))
      cat("Creating network files in:", getwd(), "\n")
      mSet <- prepare_custom_enrichnet(mSet, network_name, "mixed")
      cat("✓ Network preparation completed\n")

      # Check if network files were created
      expected_network_files <- c(
        paste0(network_name, "_network_data.json"),
        paste0(network_name, ".sif")
      )
      for (nf in expected_network_files) {
        if (file.exists(nf)) {
          cat("✓ Created network file:", nf, "\n")
        } else {
          cat("✗ Network file not found:", nf, "\n")
        }
      }
    } else {
      cat("⚠ Network preparation skipped (prepare_custom_enrichnet.R not found at:", build_network_path, ")\n")
    }

    # DEBUG: Check what's in the current directory
    cat("=== DEBUGGING FILE CREATION ===\n")
    cat("Current working directory:", getwd(), "\n")
    created_files <- list.files(getwd(), full.names = FALSE)
    cat("Files in current directory:", paste(created_files, collapse = ", "), "\n")

    # Check if specific expected files exist
    expected_files <- c("scattermum.json", "peaks_to_paths_0_dpi150.png")
    for (f in expected_files) {
      if (file.exists(f)) {
        cat("✓ Found expected file:", f, "\n")
      } else {
        cat("✗ Missing expected file:", f, "\n")
      }
    }
    cat("=== END DEBUGGING ===\n")
  })

  # Cleanup - remove temp CSV but keep all analysis outputs
  if (file.exists(temp_csv)) {
    file.remove(temp_csv)
  }
  # Change back to original working directory
  setwd(original_wd)
  cat("✓ Changed back to original directory:", getwd(), "\n")

  # Extract analysis settings from mSet object
  peak_format <- if (!is.null(mSet$paramSet$peakFormat)) mSet$paramSet$peakFormat else "rmp"
  enrich_method <- if (!is.null(mSet$paramSet$anal.type)) mSet$paramSet$anal.type else "mummichog"
  version <- if (!is.null(mSet$paramSet$version)) mSet$paramSet$version else "v2"

  # Extract actual p-value thresholds from mSet object
  peak_filter_threshold <- if (!is.null(mSet$dataSet$cutoff)) {
    mSet$dataSet$cutoff
  } else {
    1.0 # fallback
  }

  # Extract the precise 10th percentile cutoff if available
  actual_cutoff <- if (!is.null(mSet$dataSet$actual_cutoff)) {
    mSet$dataSet$actual_cutoff
  } else {
    peak_filter_threshold # fallback to rounded value
  }

  # Extract pathway enrichment statistics
  pathway_pvals <- if (!is.null(mSet$path.pval)) {
    mSet$path.pval
  } else {
    NULL
  }

  # Check for FDR-adjusted pathway p-values (if they exist)
  pathway_fdr <- if (!is.null(mSet$path.fdr)) {
    mSet$path.fdr
  } else if (!is.null(mSet$pathway.fdr)) {
    mSet$pathway.fdr
  } else if (!is.null(mSet$dataSet$pathway.fdr)) {
    mSet$dataSet$pathway.fdr
  } else {
    NULL
  }

  # Count significant pathways (those with p < 0.05)
  sig_pathways <- if (!is.null(pathway_pvals)) {
    sum(pathway_pvals < 0.05, na.rm = TRUE)
  } else {
    0
  }

  total_pathways <- if (!is.null(pathway_pvals)) {
    length(pathway_pvals)
  } else {
    0
  }

  # Extract input data statistics
  input_pvals <- if (!is.null(mSet$pvals)) {
    mSet$pvals
  } else {
    NULL
  }

  peaks_analyzed <- if (!is.null(input_pvals)) {
    length(input_pvals)
  } else {
    nrow(ttest_results)
  }

  # P-value method description (always dynamic now)
  pval_method_desc <- "Top 10% of peaks (dynamic)"

  # Create the p-value setting description (always dynamic now)
  pval_setting_desc <- "- SetMummichogPvalFromPercent: 0.1 (top 10% of peaks)\n\n"

  # Create parameter summary markdown
  md_content <- paste0(
    "# Mummichog Analysis Parameters\n\n",
    "**Analysis Date:** ", Sys.Date(), "\n\n",
    "**Database:** ", database, "\n\n",
    "**MetaboAnalystR 'Set' Function Outputs:**\n",
    "- SetPeakFormat: ", peak_format, "\n",
    "- SetPeakEnrichMethod: ", enrich_method, " (", version, ")\n",
    pval_setting_desc,
    "**Instrument Parameters (UpdateInstrumentParameters):**\n",
    "- instrumentOpt: ", instrumentOpt, "\n",
    "- msModeOpt: ", msModeOpt, "\n",
    "- force_primary_ion: ", force_primary_ion, "\n",
    "- rt_frac: ", rt_frac, "\n\n",
    "**Analysis Parameters:**\n",
    "- Peak filtering method: ", pval_method_desc, "\n",
    "- Peak filtering threshold (rounded): ", peak_filter_threshold, "\n",
    "- Peak filtering threshold (precise): ", actual_cutoff, "\n",
    "- Peaks analyzed: ", peaks_analyzed, " out of ", nrow(ttest_results), "\n",
    "- Pathways analyzed: ", total_pathways, "\n",
    "- Significant pathways (p < 0.05): ", sig_pathways, "\n",
    if (!is.null(pathway_pvals)) {
      paste0("- Pathway p-values range: ", round(min(pathway_pvals, na.rm = TRUE), 6), " to ", round(max(pathway_pvals, na.rm = TRUE), 6), "\n")
    } else {
      "- Pathway p-values: Not available\n"
    },
    if (!is.null(pathway_fdr)) {
      paste0("- Pathway FDR range: ", round(min(pathway_fdr, na.rm = TRUE), 6), " to ", round(max(pathway_fdr, na.rm = TRUE), 6), "\n")
    } else {
      "- Pathway FDR: Not calculated (using raw p-values)\n"
    },
    "- Pathway enrichment FDR threshold: 0.05 (fixed)\n",
    "- Minimum pathway size: 3\n",
    "- Background permutations: 100\n\n",
    "**Input Data:**\n",
    "- Number of features: ", nrow(ttest_results), "\n",
    "- Output directory: ", output_dir, "\n"
  )

  # Write markdown file to output directory
  md_file <- file.path(output_dir, "analysis_parameters.md")
  cat("=== WRITING MARKDOWN FILE ===\n")
  cat("Attempting to write to:", md_file, "\n")
  cat("Current working directory:", getwd(), "\n")
  cat("Original working directory:", original_wd, "\n")

  writeLines(md_content, md_file)

  # Verify the file was created
  if (file.exists(md_file)) {
    cat("✓ Markdown file successfully created at:", md_file, "\n")
    cat("File size:", file.size(md_file), "bytes\n")
  } else {
    cat("✗ ERROR: Markdown file was NOT created!\n")
  }

  # Print comprehensive p-value and analysis summary
  cat("=== ANALYSIS SUMMARY ===\n")
  cat("Peak filtering threshold (rounded):", peak_filter_threshold, "(", pval_method_desc, ")\n")
  cat("Peak filtering threshold (precise):", actual_cutoff, "\n")
  cat("Peaks analyzed:", peaks_analyzed, "out of", nrow(ttest_results), "input features\n")
  cat("Pathways analyzed:", total_pathways, "\n")
  cat("Significant pathways (p < 0.05):", sig_pathways, "\n")
  cat("Pathway enrichment FDR threshold: 0.05 (fixed by PerformPSEA)\n")
  cat("Background permutations: 100\n")
  cat("Parameter summary saved to:", md_file, "\n")

  return(mSet)
}
