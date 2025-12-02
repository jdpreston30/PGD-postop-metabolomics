#' Calculate BiVAD Dependence Status
#'
#' Determines if a patient has biventricular assist device (BiVAD) dependence
#' based on the presence of both RVAD and IABP support.
#'
#' @param clinical_metadata Data frame containing clinical metadata with MCS variables
#'
#' @return Data frame with additional column:
#'   - postop_BiVAD_dependence: "Y" if both postop_MCS_RVAD and postop_MCS_IABP are "Y", otherwise "N"
#'
#' @details
#' BiVAD dependence is defined as having both:
#' - Right Ventricular Assist Device (RVAD) support
#' - Intra-Aortic Balloon Pump (IABP) support
#'
#' Both conditions must be marked as "Y" for BiVAD dependence to be coded as "Y".
#'
#' @examples
#' \dontrun{
#' # Calculate BiVAD dependence
#' clinical_data_with_bivad <- calc_BiVAD(clinical_metadata = clinical_data)
#'
#' # View BiVAD status distribution
#' table(clinical_data_with_bivad$postop_BiVAD_dependence, useNA = "ifany")
#' }
#'
#' @export
calc_BiVAD <- function(clinical_metadata) {
  clinical_metadata_with_bivad <- clinical_metadata %>%
    mutate(
      postop_BiVAD_dependence = if_else(
        postop_MCS_RVAD == "Y" & postop_MCS_IABP == "Y",
        "Y",
        "N"
      )
    )
  return(clinical_metadata_with_bivad)
}
