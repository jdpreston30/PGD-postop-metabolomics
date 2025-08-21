#' Run individual variable PERMANOVA analysis
#'
#' This function performs PERMANOVA (Permutational Multivariate Analysis of Variance)
#' for a single variable against a distance matrix. It's designed to be used with
#' lapply to test multiple variables individually.
#'
#' @param varname Character string specifying the variable name to test
#' @param feature_data Data frame containing the feature/metabolite data (samples x features)
#' @param meta_data Data frame containing the metadata with all variables
#' @param ctrl Permutation control object from permute::how()
#' @param distance_method Character string specifying distance method for vegan::vegdist. Default "euclidean"
#' @param seed Integer for random seed to ensure reproducibility. Default 2025
#'
#' @return A tibble with columns:
#'   \item{Variable}{Character, the variable name tested}
#'   \item{R2}{Numeric, the R-squared value from PERMANOVA}
#'   \item{p_value}{Numeric, the p-value from permutation test}
#'
#' @details
#' The function:
#' \itemize{
#'   \item Sets a seed for reproducibility
#'   \item Removes samples with missing values for the test variable
#'   \item Skips variables with insufficient factor levels (< 2)
#'   \item Calculates distance matrix using specified method
#'   \item Runs adonis2 PERMANOVA test
#'   \item Returns results in tidy format
#' }
#'
#' @examples
#' \dontrun{
#' # Set up data
#' features <- UFT_data[, feature_cols]
#' metadata <- UFT_data[, c("Sex", "Age", "Variant")]
#' ctrl <- permute::how(nperm = 999)
#' 
#' # Test single variable
#' result <- get_permanova("Sex", features, metadata, ctrl)
#' 
#' # Test multiple variables
#' variables <- c("Sex", "Age", "Variant")
#' results <- bind_rows(lapply(variables, get_permanova, 
#'                            feature_data = features, 
#'                            meta_data = metadata, 
#'                            ctrl = ctrl))
#' }
#'
#' @export
get_permanova <- function(varname, 
                         feature_data, 
                         meta_data, 
                         ctrl, 
                         distance_method = "euclidean", 
                         seed = 2025) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Check that variable exists in metadata
  if (!varname %in% names(meta_data)) {
    stop(paste("Variable", varname, "not found in meta_data"))
  }
  
  # Remove samples with missing values for the test variable
  keep <- complete.cases(meta_data[[varname]])
  sub_meta <- droplevels(meta_data[keep, , drop = FALSE])
  sub_features <- feature_data[keep, , drop = FALSE]
  
  # Skip variables with insufficient levels for factors
  if (is.factor(sub_meta[[varname]]) &&
      nlevels(droplevels(sub_meta[[varname]])) < 2) {
    return(tibble::tibble(
      Variable = varname, 
      R2 = NA_real_, 
      p_value = NA_real_
    ))
  }
  
  # Calculate distance matrix
  d <- vegan::vegdist(sub_features, method = distance_method)
  
  # Run PERMANOVA
  formula_str <- paste("d ~ `", varname, "`", sep = "")
  res <- vegan::adonis2(
    as.formula(formula_str), 
    data = sub_meta, 
    permutations = ctrl
  )
  
  # Return results in tidy format
  tibble::tibble(
    Variable = varname,
    R2 = as.numeric(res$R2[1]),
    p_value = as.numeric(res$`Pr(>F)`[1])
  )
}
