#' Comment Crawler
#'
#' Extracts structured comments from R scripts in a specific order.
#' Looks for comments starting with #*, #+, or #- and compiles them
#' into a single overview file.
#'
#' @return Creates comments_check.R in the root directory
#' @export
comment_crawl <- function() {
  # Define the file order as specified
  file_order <- c(
    "R/Scripts/00a_environment_setup.r",
    "R/Scripts/00b_setup.r", 
    "R/Scripts/00c_clinical_metadata.R",
    "R/Scripts/01_FTs.r",
    "R/Scripts/02_PCA_PLSDA_heatmaps.r",
    "R/Scripts/03_limma.r",
    "R/Scripts/04_pathway_enrichment.r",
    "R/Scripts/05_targeted_volcano_diverge.r",
    "R/Scripts/06_targeted_subject_based.r",
    "R/Scripts/07_results_numbers.r",
    "R/Scripts/08_assign_figures.R",
    "R/Scripts/09_render_figures.R",
    "R/Scripts/10_tables.R",
    "R/Scripts/11_supplementary.R",
    "R/Scripts/12_session_info.R"
  )
  
  # Initialize output vector
  output_lines <- c()
  files_processed <- 0
  files_with_comments <- 0
  
  cat("🔍 Crawling through R scripts for structured comments...\n")
  
  # Process each file in order
  for (file_path in file_order) {
    if (file.exists(file_path)) {
      # Read the file
      lines <- readLines(file_path, warn = FALSE)
      
      # Find lines that start with #*, #+, or #-
      comment_pattern <- "^#[*+-]"
      comment_lines <- lines[grepl(comment_pattern, lines)]
      
      if (length(comment_lines) > 0) {
        # Add file header with #! prefix
        output_lines <- c(output_lines, paste0("#! ", file_path))
        
        # Add all matching comment lines
        output_lines <- c(output_lines, comment_lines)
        
        # Add blank line for separation between files
        output_lines <- c(output_lines, "")
        
        files_with_comments <- files_with_comments + 1
        cat("✅", file_path, "- Found", length(comment_lines), "structured comments\n")
      } else {
        cat("⚪", file_path, "- No structured comments found\n")
      }
      
      files_processed <- files_processed + 1
      
    } else {
      cat("⚠️  Warning: File not found:", file_path, "\n")
    }
  }
  
  # Write output to root directory
  output_file <- "comments_check.R"
  writeLines(output_lines, output_file)
  
  # Summary
  cat("\n📝 Comment crawl complete!\n")
  cat("📁 Output written to:", output_file, "\n")
  cat("📊 Files processed:", files_processed, "/", length(file_order), "\n")
  cat("📋 Files with structured comments:", files_with_comments, "\n")
  cat("💬 Total comment lines extracted:", length(output_lines[!grepl("^#!", output_lines) & output_lines != ""]), "\n")
  
  invisible(output_lines)
}
