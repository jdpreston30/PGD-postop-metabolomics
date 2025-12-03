#' Preprocess Feature Table for PGD Analysis
#'
#' Performs standardized preprocessing of metabolomic feature tables including
#' sample ID correction, patient ID extraction, filtering, and clinical metadata
#' integration for primary graft dysfunction (PGD) analysis.
#'
#' @param df Data frame containing raw feature table with Sample_ID column
#' @param apply_unique_filter Logical (default FALSE). If TRUE, applies unique value 
#'   filtering to remove features with low variability.
#' @param unique_threshold Numeric (default 0.7). Only used if apply_unique_filter = TRUE.
#'   Features with unique percentage >= this threshold are retained.
#'
#' @return Processed data frame with:
#'   - Corrected sample IDs
#'   - Extracted Patient and Sample identifiers  
#'   - Filtered out baseline samples (S0)
#'   - Excluded problematic patients (H49)
#'   - Integrated PGD clinical metadata
#'   - Proper factor conversions for statistical analysis
#'   - Optional unique value filtering (if apply_unique_filter = TRUE)
#'
#' @details
#' Processing steps include:
#' 1. **Sample ID correction**: Fixes known ID error (H46SS0 -> H46S0)
#' 2. **Patient extraction**: Extracts patient numbers from Sample_ID (H##)
#' 3. **Sample filtering**: Keeps only baseline samples (S0)
#' 4. **Patient filtering**: Removes patient H49 (problematic data)
#' 5. **Metadata integration**: Joins with PGD_specifics clinical data
#' 6. **Factor conversion**: Converts categorical variables to factors
#' 7. **Column organization**: Reorders columns for analysis workflow
#' 8. **Unique filtering** (optional): Removes features with low variability
#'
#' **Required external data**: PGD_specifics table must be available in environment
#' 
#' **Expected Sample_ID format**: H##S# (e.g., H01S0, H23S1)
#'
#' @examples
#' \dontrun{
#'   # Standard preprocessing without unique filtering
#'   processed_FT <- preprocess_FT(raw_feature_table)
#'   
#'   # Preprocessing with unique filtering (70% threshold)
#'   processed_FT <- preprocess_FT(raw_feature_table, apply_unique_filter = TRUE)
#'   
#'   # Preprocessing with custom unique filtering threshold
#'   processed_FT <- preprocess_FT(raw_feature_table, 
#'                                apply_unique_filter = TRUE, 
#'                                unique_threshold = 0.8)
#' }
#'
#' @export
preprocess_FT <- function(df, apply_unique_filter = FALSE, unique_threshold = 0.8) {
  df |>
    mutate(Sample_ID = if_else(Sample_ID == "H46SS0", "H46S0", Sample_ID)) |>
    mutate(
      Patient_temp = str_extract(Sample_ID, "H\\d+"),
      Patient_num = as.numeric(str_extract(Patient_temp, "\\d+")),
      Sample = str_extract(Sample_ID, "S\\d+")
    ) |>
    arrange(Patient_num, Sample) |>
    mutate(
      Patient = paste0("H", Patient_num)
    ) |>
    select(Sample_ID, Patient, Sample, everything(), -Patient_temp, -Patient_num) |>
    filter(Sample != "S0") |>
    filter(Patient != "H49") |>
    left_join(PGD_specifics, by = "Patient") |>
    select(Patient, Sample, severe_PGD, PGD_grade_tier, any_PGD, everything(), -c(Sample_ID, postop_PGD_grade_ISHLT, postop_PGD_binary_ISHLT)) |>
    mutate(PGD_grade_tier = as.factor(PGD_grade_tier)) |>
    mutate(severe_PGD = as.factor(severe_PGD)) |>
    mutate(Patient = factor(Patient)) |>
    (\(x) if (apply_unique_filter) filter_unique_features(x, unique_threshold = unique_threshold) else x)()
}
