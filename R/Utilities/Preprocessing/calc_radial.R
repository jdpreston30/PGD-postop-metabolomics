#' Calculate RADIAL Risk Score for Heart Transplant Outcomes
#'
#' Calculates the RADIAL (RAP, Age, Diabetes, Inotrope, Donor Age, and Ischemic Time) 
#' risk score for predicting primary graft dysfunction and early mortality after 
#' heart transplantation. Provides detailed missing data reporting and validation.
#'
#' @param clinical_metadata Data frame containing clinical metadata with required variables
#' @param rap_col Column name for right atrial pressure in mmHg
#' @param age_col Column name for recipient age in years
#' @param dm_col Column name for diabetes status (Y/Yes/YES for positive)
#' @param inotrope_col Column name for inotrope support status (Y/Yes/YES for positive)
#' @param donor_age_col Column name for donor age in years
#' @param ischemic_time_col Column name for ischemic time in minutes
#'
#' @return Data frame with additional column:
#'   - preop_RADIAL_calc: RADIAL risk score (0-6 scale, NA if any required data missing)
#'
#' @details
#' The RADIAL score assigns 1 point for each of the following risk factors:
#' - **R**AP ≥ 10 mmHg (Right Atrial Pressure)
#' - **A**ge ≥ 60 years (Recipient)
#' - **D**iabetes present
#' - **I**notrope support required
#' - Donor **A**ge ≥ 30 years
#' - Ischemic Time (**L**ength) ≥ 240 minutes
#' 
#' Total score ranges from 0-6, with higher scores indicating greater risk.
#' 
#' **Missing Data Handling:**
#' - Provides comprehensive missing data reporting at column and patient levels
#' - Returns NA for patients with any missing required variables
#' - Issues warnings for missing data with affected patient identification
#'
#' @examples
#' \dontrun{
#'   # Calculate RADIAL scores
#'   clinical_data_with_radial <- calc_radial(
#'     clinical_metadata = clinical_data,
#'     rap_col = "preop_RAP",
#'     age_col = "recipient_age", 
#'     dm_col = "diabetes_status",
#'     inotrope_col = "preop_inotrope_support",
#'     donor_age_col = "donor_age",
#'     ischemic_time_col = "ischemic_time_minutes"
#'   )
#'   
#'   # View score distribution
#'   table(clinical_data_with_radial$preop_RADIAL_calc, useNA = "ifany")
#' }
#'
#' @references
#' RADIAL score for risk stratification in heart transplantation outcomes
#'
#' @export
calc_radial <- function(clinical_metadata, 
                        rap_col, 
                        age_col, 
                        dm_col, 
                        inotrope_col, 
                        donor_age_col, 
                        ischemic_time_col) {
  # Check for missing values in required columns
  required_cols <- c(rap_col, age_col, dm_col, inotrope_col, donor_age_col, ischemic_time_col)
  col_names <- c("RAP", "Age", "Diabetes", "Inotrope", "Donor Age", "Ischemic Time")
  
  # Track missing data details
  missing_data_found <- FALSE
  
  for (i in seq_along(required_cols)) {
    col <- required_cols[i]
    na_indices <- which(is.na(clinical_metadata[[col]]))
    na_count <- length(na_indices)
    
    if (na_count > 0) {
      missing_data_found <- TRUE
      affected_patients <- clinical_metadata$Patient[na_indices]
      
      cat("\n=== MISSING DATA ALERT ===\n")
      cat("Column:", col_names[i], "(", col, ")\n")
      cat("Missing values:", na_count, "\n")
      cat("Affected patients:", paste(affected_patients, collapse = ", "), "\n")
      cat("========================\n")
      
      warning(paste0("Missing values detected: ", na_count, " NA(s) in ", col_names[i], 
                     " column (", col, "). Radial score calculation may be incomplete for affected patients."))
    }
  }
  
  # Summary of all missing data
  if (missing_data_found) {
    # Create a matrix to check for missing values across all required columns
    missing_matrix <- is.na(clinical_metadata[required_cols])
    
    # Calculate missing count and identify missing columns for each patient
    missing_summary <- data.frame(
      Patient = clinical_metadata$Patient,
      missing_count = rowSums(missing_matrix),
      stringsAsFactors = FALSE
    )
    
    # Add missing column names for each patient
    missing_summary$missing_columns <- apply(missing_matrix, 1, function(row) {
      if (sum(row) > 0) {
        paste(col_names[which(row)], collapse = ", ")
      } else {
        ""
      }
    })
    
    # Filter to only patients with missing data
    missing_summary <- missing_summary[missing_summary$missing_count > 0, ]
    
    cat("\n=== PATIENT-LEVEL MISSING DATA SUMMARY ===\n")
    print(missing_summary)
    cat("==========================================\n\n")
  }
  
  # Calculate RADIAL score
  clinical_metadata_with_radial <- clinical_metadata |>
    mutate(
      preop_RADIAL_calc = ifelse(
        # Check if ANY required values are missing for this patient
        is.na(.data[[rap_col]]) | 
        is.na(.data[[age_col]]) | 
        is.na(.data[[dm_col]]) | 
        is.na(.data[[inotrope_col]]) | 
        is.na(.data[[donor_age_col]]) | 
        is.na(.data[[ischemic_time_col]]),
        NA,  # Return NA if any values are missing
        # Otherwise calculate the RADIAL score
        rowSums(
          cbind(
            .data[[rap_col]] >= 10,      # RAP ≥ 10 mmHg
            .data[[age_col]] >= 60,      # Age ≥ 60 years
            .data[[dm_col]] %in% c("Y", "Yes", "YES"), # Diabetes = Yes
            .data[[inotrope_col]] %in% c("Y", "Yes", "YES"), # Inotrope = Yes
            .data[[donor_age_col]] >= 30, # Donor age ≥ 30 years
            .data[[ischemic_time_col]] >= 240 # Ischemic time ≥ 240 min
          )
        )
      )
    )
  
  cat("RADIAL scores calculated successfully!\n")
  cat("RADIAL score range:", min(clinical_metadata_with_radial$preop_RADIAL_calc, na.rm = TRUE), 
      "to", max(clinical_metadata_with_radial$preop_RADIAL_calc, na.rm = TRUE), "\n\n")
  
  return(clinical_metadata_with_radial)
}
