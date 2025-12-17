#' Find features matching target m/z values within a specified ppm tolerance
#'
#' This function searches for features that match target m/z values using 
#' ppm-based mass tolerance, following the same logic as create_identified_FT.
#' Retention time is not considered in the matching process.
#'
#' @param feature_df A tibble with columns: feature, column, mz, rt
#' @param target_mz_vec A numeric vector of target m/z values to search for
#' @param ppm_tolerance Mass tolerance in parts per million (default: 10)
#' @return A tibble of matching features with ppm error calculations
#'
#' @examples
#' \dontrun{
#'   target_mz <- c(71.0126, 53.0021, 68.997)
#'   matched <- find_features_by_mz(feature_names, target_mz, ppm_tolerance = 10)
#' }
#'
#' @export
find_features_by_mz <- function(feature_df, target_mz_vec, ppm_tolerance = 10) {
  
  # Load required libraries
  library(dplyr)
  library(purrr)
  
  # Initialize empty list to store matches
  all_matches <- list()
  
  # Iterate through each target m/z
  for (i in seq_along(target_mz_vec)) {
    target <- target_mz_vec[i]
    
    # Calculate m/z threshold in Daltons (ppm-based)
    # This matches the logic in create_identified_FT
    mz_threshold_da <- (ppm_tolerance * target) / 1000000
    
    # Find features within threshold
    matching_indices <- which(abs(feature_df$mz - target) <= mz_threshold_da)
    
    if (length(matching_indices) > 0) {
      # Create match results for this target
      for (idx in matching_indices) {
        match_result <- feature_df[idx, ] |>
          mutate(
            target_mz = target,
            mz_error_ppm = ((mz - target) / target) * 1000000,
            mz_error_abs_ppm = abs(mz_error_ppm)
          )
        
        all_matches <- append(all_matches, list(match_result))
      }
    }
  }
  
  # Combine all matches
  if (length(all_matches) > 0) {
    matched_features <- bind_rows(all_matches) |>
      arrange(target_mz, mz_error_abs_ppm)
    
    return(matched_features)
  } else {
    return(tibble())
  }
}
