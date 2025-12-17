#' Ensure renv Environment is Ready
#'
#' Checks if the renv package environment is properly set up and restores
#' packages if needed. This function should be called at the beginning of
#' the analysis pipeline to ensure all dependencies are available.
#'
#' @details
#' The function performs the following checks:
#' 1. Verifies renv is active
#' 2. Checks if critical packages are available
#' 3. Automatically runs renv::restore() if packages are missing
#' 4. Uses renv::status() to detect any inconsistencies
#'
#' @return Invisible NULL. Prints status messages to console.
#' @export
ensure_renv_ready <- function() {
  # Activate renv if not already active
  if (!("renv" %in% loadedNamespaces())) {
    cat("📦 Activating renv environment...\n")
    source("renv/activate.R")
  }
  
  # Check critical packages that are always needed
  critical_packages <- c("yaml", "here", "dplyr", "readr", "ggplot2")
  missing_critical <- critical_packages[!sapply(critical_packages, requireNamespace, quietly = TRUE)]
  
  # If critical packages are missing, definitely need to restore
  if (length(missing_critical) > 0) {
    cat("⚠️  Critical packages missing:", paste(missing_critical, collapse = ", "), "\n")
    cat("🔄 Installing packages from renv lockfile...\n")
    cat("   (This may take 10-20 minutes on first run)\n\n")
    
    renv::restore(prompt = FALSE)
    cat("\n✅ Package installation complete!\n\n")
    
    # Verify restoration succeeded
    still_missing <- critical_packages[!sapply(critical_packages, requireNamespace, quietly = TRUE)]
    if (length(still_missing) > 0) {
      stop("❌ Failed to restore packages: ", paste(still_missing, collapse = ", "),
           "\n   Run renv::status() for details or try renv::restore() manually.")
    }
    
    return(invisible(NULL))
  }
  
  # Even if critical packages exist, check if library is out of sync
  # This catches cases where some packages are installed but others aren't
  status <- tryCatch({
    capture.output(renv::status(), type = "message")
  }, error = function(e) NULL)
  
  # Look for signs that packages need restoration
  if (!is.null(status) && any(grepl("inconsistent|out-of-sync|not installed", status, ignore.case = TRUE))) {
    cat("⚠️  renv library is out of sync with lockfile\n")
    cat("🔄 Synchronizing package library...\n")
    cat("   (This may take several minutes)\n\n")
    
    renv::restore(prompt = FALSE)
    cat("\n✅ Library synchronized!\n\n")
  } else {
    cat("✅ renv environment ready - all packages available\n\n")
  }
  
  invisible(NULL)
}

# Auto-run when sourced (but can be disabled by setting option)
if (!isTRUE(getOption("renv.helper.noautorun"))) {
  ensure_renv_ready()
}
