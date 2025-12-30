#' Find Other Adducts Matched
#'
#' Determines if compounds are detected with different adducts within the same 
#' chromatographic mode (C18 or HILIC) based on hierarchical matching criteria:
#' exact name match, KEGG ID, then HMDB ID. Formula matching is excluded.
#' Once a "Y" is assigned, it's never overridden to "N".
#'
#' @param matched_features Tibble with columns: compound_name, feature_mode, adduct, kegg_id, hmdb_id
#' @return Tibble with added/updated other_adducts_matched column ("Y" or "N")
#' @export
#'
#' @examples
#' \dontrun{
#'   # After create_identified_FT() and calculate_MMD()
#'   TFT_confirmed_key <- identified$matched_features |>
#'     filter(!is.na(compound_name)) |>
#'     calculate_MMD() |>
#'     adduct_finder() |>
#'     select(identified_name = compound_name, everything())
#' }
adduct_finder <- function(matched_features) {
  
  # Initialize other_adducts_matched column as "N"
  result <- matched_features |>
    mutate(other_adducts_matched = "N")
  
  # Helper function to check if compound detected with different adducts within same mode
  check_different_adducts <- function(df, match_col) {
    df |>
      filter(!is.na(!!sym(match_col)), !!sym(match_col) != "") |>
      group_by(!!sym(match_col), feature_mode) |>
      summarise(
        adducts = list(unique(adduct)),
        n_adducts = n_distinct(adduct),
        .groups = "drop"
      ) |>
      filter(n_adducts > 1) |>  # Different adducts within same mode
      select(!!sym(match_col), feature_mode) |>
      # Return compound-mode combinations that have multiple adducts
      unite("compound_mode", !!sym(match_col), feature_mode, sep = "|||") |>
      pull(compound_mode)
  }
  
  # 1. Exact compound name matches (permanent Y)
  compound_diff_adducts <- check_different_adducts(result, "compound_name")
  result <- result |>
    unite("temp_compound_mode", compound_name, feature_mode, sep = "|||", remove = FALSE) |>
    mutate(other_adducts_matched = if_else(temp_compound_mode %in% compound_diff_adducts, "Y", other_adducts_matched)) |>
    select(-temp_compound_mode)
  
  # 2. KEGG ID matches (only update N -> Y, never Y -> N)
  kegg_diff_adducts <- check_different_adducts(result, "kegg_id")
  result <- result |>
    unite("temp_kegg_mode", kegg_id, feature_mode, sep = "|||", remove = FALSE) |>
    mutate(other_adducts_matched = if_else(other_adducts_matched == "N" & temp_kegg_mode %in% kegg_diff_adducts, "Y", other_adducts_matched)) |>
    select(-temp_kegg_mode)
  
  # 3. HMDB ID matches (only update N -> Y, never Y -> N)
  hmdb_diff_adducts <- check_different_adducts(result, "hmdb_id")
  result <- result |>
    unite("temp_hmdb_mode", hmdb_id, feature_mode, sep = "|||", remove = FALSE) |>
    mutate(other_adducts_matched = if_else(other_adducts_matched == "N" & temp_hmdb_mode %in% hmdb_diff_adducts, "Y", other_adducts_matched)) |>
    select(-temp_hmdb_mode)
  
  # Note: Formula matching is excluded as requested
  
  return(result)
}