#' Create Formatted Excel Export for Supplementary Materials with Multiple Tabs
#'
#' This function exports two data frames to Excel with professional formatting including
#' gray headers, bold text, cell borders, and auto-sized columns for all data.
#'
#' @param data1 Data frame for first tab (Level 1)
#' @param data2 Data frame for second tab (Level 3/Annotations)
#' @param filename Character; output file path
#' @param sheet1_name Character; first worksheet name (default: "Level 1")
#' @param sheet2_name Character; second worksheet name (default: "Level 3/Annotations")
#'
#' @return Invisibly returns the file path
#' @export
create_SM2 <- function(data1, data2, filename, sheet1_name = "Level 1", sheet2_name = "Level 3") {
  
  # Load required library
  library(openxlsx)
  
  # Create workbook with minimal options
  wb <- createWorkbook()
  
  # Add worksheets
  addWorksheet(wb, sheet1_name)
  addWorksheet(wb, sheet2_name)
  
  # Write data only - no formatting at all
  writeData(wb, sheet1_name, data1)
  writeData(wb, sheet2_name, data2)
  
  # Create simple styles (conservative approach)
  header_style <- createStyle(
    textDecoration = "bold",
    fgFill = "lightgray",
    border = c("top", "bottom", "left", "right")
  )
  
  data_style <- createStyle(
    border = c("top", "bottom", "left", "right")
  )
  
  # Apply header styling
  addStyle(wb, sheet1_name, header_style, rows = 1, cols = 1:ncol(data1), gridExpand = TRUE)
  addStyle(wb, sheet2_name, header_style, rows = 1, cols = 1:ncol(data2), gridExpand = TRUE)
  
  # Apply data borders
  if(nrow(data1) > 0) {
    addStyle(wb, sheet1_name, data_style, rows = 2:(nrow(data1) + 1), cols = 1:ncol(data1), gridExpand = TRUE)
  }
  if(nrow(data2) > 0) {
    addStyle(wb, sheet2_name, data_style, rows = 2:(nrow(data2) + 1), cols = 1:ncol(data2), gridExpand = TRUE)
  }
  
  # Auto-fit column widths
  setColWidths(wb, sheet1_name, cols = 1:ncol(data1), widths = "auto")
  setColWidths(wb, sheet2_name, cols = 1:ncol(data2), widths = "auto")
  
  # Save with minimal options
  saveWorkbook(wb, filename, overwrite = TRUE)
  
  invisible(filename)
}
