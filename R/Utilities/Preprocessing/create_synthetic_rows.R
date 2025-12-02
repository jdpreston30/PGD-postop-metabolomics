#' Create Synthetic Rows for Missing Timepoints
#'
#' @param feature_table A data frame containing the feature table with Patient, Sample, severe_PGD, PGD_grade_tier, and any_PGD columns
#' @param missing_samples A data frame with columns: Patient, Sample, severe_PGD, PGD_grade_tier, any_PGD
#' @return Feature table with synthetic rows added
#' @export
create_synthetic_rows <- function(feature_table, missing_samples) {
  
  # Create synthetic rows for each missing sample
  synthetic_rows <- lapply(seq_len(nrow(missing_samples)), function(i) {
    missing <- missing_samples[i, ]
    
    # Calculate mean of numeric features from same severe_PGD group and Sample timepoint
    feature_table %>%
      filter(severe_PGD == missing$severe_PGD, Sample == missing$Sample) %>%
      select(-Patient, -Sample) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
      mutate(
        Patient = as.character(missing$Patient),
        Sample = as.character(missing$Sample),
        severe_PGD = as.character(missing$severe_PGD),
        PGD_grade_tier = as.character(missing$PGD_grade_tier),
        any_PGD = as.character(missing$any_PGD)
      )
  }) %>%
    bind_rows()
  
  # Bind synthetic rows and sort
  # bind_rows will coerce factors to match existing levels automatically
  feature_table %>%
    bind_rows(synthetic_rows) %>%
    arrange(Patient, Sample)
}
