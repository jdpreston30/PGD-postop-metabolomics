#' Create Supplementary Material 3 Excel File with Conditional Formatting
#'
#' Creates a formatted Excel file for Supplementary Material 3 with highlighting
#' for significant p-values and corresponding metabolite names.
#'
#' @param SM3_ttest Data frame for t-test results (first tab)
#' @param SM3_limma Data frame for LIMMA results (second tab) 
#' @param file_path Path where the Excel file should be saved
#'
#' @return Invisibly returns the workbook object
#' @export
create_SM3 <- function(SM3_ttest, SM3_limma, file_path) {
  
  # Create workbook
  wb <- createWorkbook()
  
  # Add sheets
  addWorksheet(wb, "T-test Results")
  addWorksheet(wb, "LIMMA Results")
  
  # Write data to sheets
  writeData(wb, sheet = 1, SM3_ttest, startRow = 1, startCol = 1)
  writeData(wb, sheet = 2, SM3_limma, startRow = 1, startCol = 1)
  
  # Create styles to match create_SM2
  header_style <- createStyle(
    textDecoration = "bold",
    fgFill = "lightgray",
    border = c("top", "bottom", "left", "right")
  )
  
  data_style <- createStyle(
    border = c("top", "bottom", "left", "right")
  )
  
  significant_highlight <- createStyle(
    fgFill = "#FFB6C1",  # Light pink for significant rows
    border = c("top", "bottom", "left", "right")  # Keep borders
  )
  pvalue_highlight <- createStyle(
    fgFill = "#FFB6C1",  # Same light pink for significant p-values
    border = c("top", "bottom", "left", "right")  # Keep borders
  )
  
  # Apply header and data formatting to both sheets
  for (sheet_num in 1:2) {
    data_to_use <- if (sheet_num == 1) SM3_ttest else SM3_limma
    
    # Header formatting (grey background, bold, borders)
    addStyle(wb, sheet = sheet_num, header_style, 
             rows = 1, cols = 1:ncol(data_to_use), gridExpand = TRUE)
    
    # Data borders (universal borders on all data cells)
    if(nrow(data_to_use) > 0) {
      addStyle(wb, sheet = sheet_num, data_style, 
               rows = 2:(nrow(data_to_use) + 1), cols = 1:ncol(data_to_use), gridExpand = TRUE)
    }
    
    # Auto-adjust column widths
    setColWidths(wb, sheet = sheet_num, cols = 1:ncol(data_to_use), widths = "auto")
  }
  
  # Apply number formatting to LIMMA sheet p-value columns
  if (nrow(SM3_limma) > 0) {
    p_interaction_col <- which(names(SM3_limma) == "P (Interaction)")
    p_pgd_col <- which(names(SM3_limma) == "P (PGD)")  
    p_time_col <- which(names(SM3_limma) == "P (Time)")
    
    # Create scientific notation style for very small numbers, otherwise 3 decimals
    # Include borders so they don't get lost when number formatting is applied
    pvalue_style <- createStyle(
      numFmt = "0.000",
      border = c("top", "bottom", "left", "right")
    )
    
    # Apply to p-value columns
    if (length(p_interaction_col) > 0) {
      addStyle(wb, sheet = 2, pvalue_style, rows = 2:(nrow(SM3_limma) + 1), 
               cols = p_interaction_col, gridExpand = TRUE)
    }
    if (length(p_pgd_col) > 0) {
      addStyle(wb, sheet = 2, pvalue_style, rows = 2:(nrow(SM3_limma) + 1), 
               cols = p_pgd_col, gridExpand = TRUE)
    }
    if (length(p_time_col) > 0) {
      addStyle(wb, sheet = 2, pvalue_style, rows = 2:(nrow(SM3_limma) + 1), 
               cols = p_time_col, gridExpand = TRUE)
    }
  }
  
  # Special conditional formatting for both sheets
  
  # T-test sheet (sheet 1) - highlight p-values < 0.05
  if (nrow(SM3_ttest) > 0) {
    # Find column indices for p-values in t-test sheet
    p_12h_col <- which(names(SM3_ttest) == "P (12h)")
    p_24h_col <- which(names(SM3_ttest) == "P (24h)")
    
    # Manually highlight p-value cells < 0.05
    for (row in 1:nrow(SM3_ttest)) {
      data_row <- row + 1  # +1 for header
      
      # Check and highlight P (12h)
      if (length(p_12h_col) > 0 && !is.na(SM3_ttest[row, p_12h_col]) && SM3_ttest[row, p_12h_col] < 0.05) {
        addStyle(wb, sheet = 1, pvalue_highlight, rows = data_row, cols = p_12h_col)
      }
      
      # Check and highlight P (24h) 
      if (length(p_24h_col) > 0 && !is.na(SM3_ttest[row, p_24h_col]) && SM3_ttest[row, p_24h_col] < 0.05) {
        addStyle(wb, sheet = 1, pvalue_highlight, rows = data_row, cols = p_24h_col)
      }
    }
  }
  
  # LIMMA sheet (sheet 2) - highlight p-values < 0.05
  if (nrow(SM3_limma) > 0) {
    
    # Find column indices for p-values
    p_interaction_col <- which(names(SM3_limma) == "P (Interaction)")
    p_pgd_col <- which(names(SM3_limma) == "P (PGD)")  
    p_time_col <- which(names(SM3_limma) == "P (Time)")
    
    # Manually highlight p-value cells < 0.05 (conditional formatting wasn't working)
    for (row in 1:nrow(SM3_limma)) {
      data_row <- row + 1  # +1 for header
      
      # Check and highlight P (Interaction)
      if (length(p_interaction_col) > 0 && !is.na(SM3_limma[row, p_interaction_col]) && SM3_limma[row, p_interaction_col] < 0.05) {
        addStyle(wb, sheet = 2, pvalue_highlight, rows = data_row, cols = p_interaction_col)
      }
      
      # Check and highlight P (PGD) 
      if (length(p_pgd_col) > 0 && !is.na(SM3_limma[row, p_pgd_col]) && SM3_limma[row, p_pgd_col] < 0.05) {
        addStyle(wb, sheet = 2, pvalue_highlight, rows = data_row, cols = p_pgd_col)
      }
      
      # Check and highlight P (Time)
      if (length(p_time_col) > 0 && !is.na(SM3_limma[row, p_time_col]) && SM3_limma[row, p_time_col] < 0.05) {
        addStyle(wb, sheet = 2, pvalue_highlight, rows = data_row, cols = p_time_col)
      }
    }
  }
  
  # Save workbook
  saveWorkbook(wb, file_path, overwrite = TRUE)
  
  message(sprintf("Supplementary Material 3 saved to: %s", file_path))
  
  invisible(wb)
}