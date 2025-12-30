#' Calculate Multi-Mode Detection (MMD)
#'
#' Determines if compounds are detected in both C18 and HILIC modes based on
#' hierarchical matching criteria: exact name match, KEGG ID, HMDB ID, then Formula.
#' Once a "Y" is assigned, it's never overridden to "N".
#'
#' @param matched_features Tibble with columns: compound_name, feature_mode, kegg_id, hmdb_id, formula
#' @return Tibble with added/updated MMD column ("Y" or "N")
#' @export
#'
#' @examples
#' \dontrun{
#'   # After create_identified_FT()
#'   TFT_confirmed_key <- identified$matched_features |>
#'     filter(!is.na(compound_name)) |>
#'     calculate_MMD() |>
#'     select(identified_name = compound_name, everything())
#' }
calculate_MMD <- function(matched_features) {
  
  # Initialize MMD column as "N"
  result <- matched_features |>
    mutate(MMD = "N")
  
  # Helper function to check if compound detected in both modes
  check_both_modes <- function(df, match_col) {
    df |>
      filter(!is.na(!!sym(match_col)), !!sym(match_col) != "") |>
      group_by(!!sym(match_col)) |>
      summarise(
        modes = list(unique(feature_mode)),
        n_modes = n_distinct(feature_mode),
        .groups = "drop"
      ) |>
      filter(n_modes > 1) |>  # Detected in multiple modes
      pull(!!sym(match_col))
  }
  
  # 1. Exact compound name matches (permanent Y)
  compound_both_modes <- check_both_modes(result, "compound_name")
  result <- result |>
    mutate(MMD = if_else(compound_name %in% compound_both_modes, "Y", MMD))
  
  # 2. KEGG ID matches (only update N -> Y, never Y -> N)
  kegg_both_modes <- check_both_modes(result, "kegg_id")
  result <- result |>
    mutate(MMD = if_else(MMD == "N" & kegg_id %in% kegg_both_modes, "Y", MMD))
  
  # 3. HMDB ID matches (only update N -> Y, never Y -> N)
  hmdb_both_modes <- check_both_modes(result, "hmdb_id")
  result <- result |>
    mutate(MMD = if_else(MMD == "N" & hmdb_id %in% hmdb_both_modes, "Y", MMD))
  
  # 4. Formula matches (only update N -> Y, never Y -> N)
  formula_both_modes <- check_both_modes(result, "formula")
  result <- result |>
    mutate(MMD = if_else(MMD == "N" & formula %in% formula_both_modes, "Y", MMD))
  
  return(result)
}