#' Calculate Predicted Heart Mass (PHM) for Donors
#'
#' Calculates the predicted heart mass for organ donors using sex-specific regression
#' equations based on donor age, height, and weight. Used for donor-recipient matching
#' and size compatibility assessment in heart transplantation.
#'
#' @param clinical_metadata Data frame containing clinical metadata with donor characteristics
#' @param donor_age_col Column name for donor age in years (default: "donor_age")
#' @param donor_sex_col Column name for donor sex (default: "donor_sex")
#' @param donor_weight_col Column name for donor weight in kg (default: "donor_weight")
#' @param donor_height_col Column name for donor height in cm (default: "donor_height")
#'
#' @return Data frame with additional column:
#'   - donor_PHM_calc: Calculated predicted heart mass in grams
#'
#' @details
#' The predicted heart mass is calculated using sex-specific equations:
#' 
#' **For Female donors:**
#' PHM = (6.82 × (height/100)^0.54 × weight^0.61) + (10.59 × age^-0.32 × (height/100)^1.135 × weight^0.315)
#' 
#' **For Male donors (and all others):**
#' PHM = (8.25 × (height/100)^0.54 × weight^0.61) + (11.25 × age^-0.32 × (height/100)^1.135 × weight^0.315)
#' 
#' Where:
#' - height is in cm (converted to meters for calculation)
#' - weight is in kg
#' - age is in years
#'
#' @examples
#' \dontrun{
#'   # Calculate predicted heart mass for donors
#'   clinical_data_with_phm <- calc_PHM(clinical_metadata)
#'   
#'   # View distribution of predicted heart masses
#'   summary(clinical_data_with_phm$donor_PHM_calc)
#' }
#'
#' @references
#' Regression equations for predicted heart mass based on donor anthropometric data
#'
#' @export
calc_PHM <- function(clinical_metadata, 
                     donor_age_col = "donor_age",
                     donor_sex_col = "donor_sex",
                     donor_weight_col = "donor_weight",
                     donor_height_col = "donor_height") {

  clinical_metadata <- clinical_metadata |>
    mutate(
      donor_PHM_calc = case_when(
        .data[[donor_sex_col]] == "Female" ~ (6.82 * (.data[[donor_height_col]] / 100)^0.54 * .data[[donor_weight_col]]^0.61) + (10.59 * .data[[donor_age_col]]^-0.32 * (.data[[donor_height_col]] / 100)^1.135 * .data[[donor_weight_col]]^0.315),
        TRUE ~ (8.25 * (.data[[donor_height_col]] / 100)^0.54 * .data[[donor_weight_col]]^0.61) + (11.25 * .data[[donor_age_col]]^-0.32 * (.data[[donor_height_col]] / 100)^1.135 * .data[[donor_weight_col]]^0.315)
      )
    )

  return(clinical_metadata)
}