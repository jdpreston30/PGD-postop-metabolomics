#' Calculate ISHLT Primary Graft Dysfunction (PGD) Classification
#'
#' Calculates Primary Graft Dysfunction grade according to International Society for 
#' Heart and Lung Transplantation (ISHLT) criteria based on mechanical support 
#' requirements and hemodynamic parameters within 24 hours post-transplant.
#'
#' @param clinical_metadata_i Data frame containing clinical metadata with postoperative variables
#' @param impella_dep_col Column name for Impella 5.5 dependence (default: "postop_MCS_Impella5.5_DEPENDENT")
#' @param rvad_col Column name for right ventricular assist device use (default: "postop_MCS_RVAD")
#' @param iabp_col Column name for intra-aortic balloon pump use (default: "postop_MCS_IABP")
#' @param ecmo_col Column name for VA-ECMO use (default: "postop_VA_ECMO")
#' @param cvp_col Column name for central venous pressure (default: "postop_CVP")
#' @param cardiac_index_col Column name for cardiac index (default: "postop_cardiac_index")
#' @param lvef_col Column name for left ventricular ejection fraction (default: "postop_LVEF_median")
#' @param inotrope_col Column name for inotrope score (default: "postop_inotrope_score")
#'
#' @return Data frame with additional columns:
#'   - postop_PGD_grade_ISHLT: ISHLT PGD grade ("N", "Mild", "Moderate", "Severe", or "")
#'   - postop_PGD_binary_ISHLT: Binary PGD classification ("N", "Y", or "")
#'
#' @details
#' ISHLT PGD Classification Criteria:
#' - **Severe PGD**: Any mechanical circulatory support (VA-ECMO, RVAD, or Impella dependence)
#' - **Moderate PGD**: Hemodynamic compromise AND (high inotrope score >10 OR IABP use)
#' - **Mild PGD**: Hemodynamic compromise AND moderate inotrope support (1-10)
#' - **No PGD**: No criteria met
#' 
#' Hemodynamic compromise defined as any of:
#' - LVEF ≤ 40%
#' - CVP > 15 mmHg  
#' - Cardiac Index < 2.0 L/min/m²
#'
#' @examples
#' \dontrun{
#'   # Calculate ISHLT PGD classification
#'   clinical_data_with_pgd <- calc_ISHLT(clinical_metadata)
#'   
#'   # View PGD distribution
#'   table(clinical_data_with_pgd$postop_PGD_grade_ISHLT)
#' }
#'
#' @references
#' International Society for Heart and Lung Transplantation Primary Graft Dysfunction 
#' consensus guidelines
#'
#' @export
calc_ISHLT <- function(clinical_metadata_i,
                       impella_dep_col = "postop_MCS_Impella5.5_DEPENDENT",
                       rvad_col = "postop_MCS_RVAD",
                       iabp_col = "postop_MCS_IABP",
                       ecmo_col = "postop_VA_ECMO",
                       cvp_col = "postop_CVP",
                       cardiac_index_col = "postop_cardiac_index",
                       lvef_col = "postop_LVEF_median",
                       inotrope_col = "postop_inotrope_score") {
  
  clinical_metadata <- clinical_metadata_i
  original_PGD <- clinical_metadata_i %>% select(postop_PGD_textbook_calc, postop_PGD_ISHLT)

  # Check column names exist
  missing_cols <- c()
  for (col in c(impella_dep_col,
                rvad_col,
                iabp_col,
                ecmo_col,
                cvp_col,
                cardiac_index_col,
                lvef_col,
                inotrope_col)) {
    if (!col %in% colnames(clinical_metadata)) {
      missing_cols <- c(missing_cols, col)
    }
  }
  
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following required columns are missing from the input data:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  # Calculate ISHLT PGD Grade
  clinical_metadata <- clinical_metadata %>%
    mutate(
      postop_PGD_grade_ISHLT = case_when(
        # Check if all values are empty/NA (return empty string)
        is.na(.data[[impella_dep_col]]) & is.na(.data[[rvad_col]]) & 
        is.na(.data[[ecmo_col]]) & is.na(.data[[iabp_col]]) & 
        is.na(.data[[cvp_col]]) & is.na(.data[[cardiac_index_col]]) & 
        is.na(.data[[lvef_col]]) & 
        (is.na(.data[[inotrope_col]]) | .data[[inotrope_col]] == 0) ~ "",
        
        # Check for Severe PGD first (any mechanical support)
        .data[[ecmo_col]] == "Y" | 
        .data[[rvad_col]] == "Y" | 
        .data[[impella_dep_col]] == "Y" ~ "Severe",
        
        # Moderate PGD: hemodynamic compromise AND (high inotrope OR IABP)
        ((.data[[lvef_col]] <= 40 & !is.na(.data[[lvef_col]])) |
         .data[[cvp_col]] > 15 |
         (if_else(is.na(.data[[cardiac_index_col]]), 1, .data[[cardiac_index_col]]) < 2)) &
        (.data[[inotrope_col]] > 10 | .data[[iabp_col]] == "Y") ~ "Moderate",
        
        # Mild PGD: hemodynamic compromise AND inotrope 1-10
        ((.data[[lvef_col]] <= 40 & !is.na(.data[[lvef_col]])) |
         .data[[cvp_col]] > 15 |
         (if_else(is.na(.data[[cardiac_index_col]]), 1, .data[[cardiac_index_col]]) < 2)) &
        (.data[[inotrope_col]] > 0 & .data[[inotrope_col]] <= 10) ~ "Mild",
        
        # Default: No PGD
        TRUE ~ "N"
      ),
      
      postop_PGD_binary_ISHLT = case_when(
        postop_PGD_grade_ISHLT == "" ~ "",
        postop_PGD_grade_ISHLT == "N" ~ "N", 
        TRUE ~ "Y"
      )
    )

  return(clinical_metadata)
}
