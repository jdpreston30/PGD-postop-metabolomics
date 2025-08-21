#' Perform t-tests between groups and format results for Mummichog analysis
#'
#' @param data Data frame with Patient_ID, Variant, and feature columns
#' @param group_assignments Data frame with Patient_ID and group assignments (e.g., Clade)
#' @param group_column Name of the group column in group_assignments (default: "Clade")
#' @param output_filename Filename for the exported CSV results
#' @param group1_value Value representing the first group (default: 1)
#' @param group2_value Value representing the second group (default: 2)
#' @return List containing the results tibble and summary statistics
#' @export
mummichog_ttests <- function(data,
                            group_assignments,
                            group_column = "Clade",
                            output_filename,
                            group1_value = 1,
                            group2_value = 2) {
  
  #_Prepare data for t-tests with group assignments
  ttest_data <- data %>%
    dplyr::left_join(
      group_assignments %>% dplyr::mutate(Group_Test = !!rlang::sym(group_column)),
      by = "Patient_ID"
    ) %>%
    dplyr::select(-Patient_ID) %>%
    dplyr::filter(!is.na(Group_Test))
  
  #_Remove the original grouping column if it exists (and is different from Group_Test)
  if (group_column %in% names(ttest_data) && group_column != "Group_Test") {
    ttest_data <- ttest_data %>% dplyr::select(-!!rlang::sym(group_column))
  }
  
  #_Get feature names (exclude Group_Test and any remaining non-numeric columns)
  feature_names <- names(ttest_data)[names(ttest_data) != "Group_Test"]
  
  #_Filter to only include columns that look like metabolite features (HILIC or C18 prefix)
  metabolite_feature_names <- feature_names[stringr::str_starts(feature_names, "HILIC|C18")]
  
  #_Remove any remaining categorical columns that shouldn't be tested
  #_Only keep columns that can be converted to numeric for t-tests
  numeric_feature_names <- c()
  for (col in metabolite_feature_names) {
    test_values <- ttest_data[[col]]
    if (is.numeric(test_values) || (!is.factor(test_values) && !all(is.na(suppressWarnings(as.numeric(as.character(test_values))))))) {
      numeric_feature_names <- c(numeric_feature_names, col)
    }
  }
  feature_names <- numeric_feature_names
  
  cat("Total columns in data:", ncol(ttest_data), "\n")
  cat("Metabolite feature columns (HILIC/C18):", length(metabolite_feature_names), "\n")
  cat("Numeric feature columns for t-tests:", length(feature_names), "\n")
  
  #_Initialize results list
  ttest_results <- list()
  
  #_Perform t-test for each feature between groups
  cat("Performing t-tests for", length(feature_names), "features...\n")
  
  for (i in seq_along(feature_names)) {
    feature <- feature_names[i]
    
    #_Progress indicator
    if (i %% 5000 == 0) {
      cat("Processed", i, "of", length(feature_names), "features...\n")
    }
    
    group1_values <- ttest_data[ttest_data$Group_Test == group1_value, feature]
    group2_values <- ttest_data[ttest_data$Group_Test == group2_value, feature]
    
    #_Convert to numeric if factor and remove NA values
    if (is.factor(group1_values)) group1_values <- as.numeric(as.character(group1_values))
    if (is.factor(group2_values)) group2_values <- as.numeric(as.character(group2_values))
    group1_values <- group1_values[!is.na(group1_values)]
    group2_values <- group2_values[!is.na(group2_values)]
    
    #_Check if there's sufficient data for t-test
    if (length(group1_values) < 2 || length(group2_values) < 2) {
      ttest_results[[feature]] <- 1.0  # Assign p-value of 1 for insufficient data
      next
    }
    
    #_Check for constant data (no variance) - assign p-value of 1
    #_Use a safer method to check for constant data
    group1_constant <- length(unique(group1_values)) == 1
    group2_constant <- length(unique(group2_values)) == 1
    if (group1_constant && group2_constant) {
      ttest_results[[feature]] <- 1.0  # Assign p-value of 1 for constant data
      next
    }
    
    #_Perform t-test with error handling
    tryCatch({
      test_result <- t.test(group1_values, group2_values)
      ttest_results[[feature]] <- test_result$p.value
    }, error = function(e) {
      cat("Warning: Assigning p-value of 1 to feature", feature, "due to error:", e$message, "\n")
      ttest_results[[feature]] <- 1.0  # Assign p-value of 1 for errors
    })
  }
  
  #_Parse feature names and create results tibble
  results_tibble <- tibble::tibble(
    Feature = names(ttest_results),
    p.value = unlist(ttest_results)
  ) %>%
    dplyr::mutate(
      #_Extract mode and convert: HILIC -> pos, C18 -> neg
      mode = dplyr::case_when(
        stringr::str_starts(Feature, "HILIC") ~ "positive",
        stringr::str_starts(Feature, "C18") ~ "negative",
        TRUE ~ NA_character_
      ),
      #_Extract m.z (first number after underscore)
      m.z = stringr::str_extract(Feature, "(?<=_)[0-9.]+"),
      #_Extract rt (second number - after second underscore)
      r.t = stringr::str_extract(Feature, "_[0-9.]+_([0-9.]+)") %>%
        stringr::str_extract("[0-9.]+$")
    ) %>%
    dplyr::select(m.z, p.value, mode, r.t) %>%
    dplyr::mutate(
      m.z = as.numeric(m.z),
      r.t = as.numeric(r.t)
    ) %>%
    #_Remove rows where feature parsing failed (invalid feature names)
    dplyr::filter(!is.na(m.z) & !is.na(mode))
  
  #_Export results
  readr::write_csv(results_tibble, paste0("Outputs/", output_filename))
  
  #_Display summary
  cat("T-test results exported to:", paste0("Outputs/", output_filename), "\n")
  cat("Total features processed:", length(feature_names), "\n")
  cat("Features with actual t-test p-values:", sum(unlist(ttest_results) < 1.0), "\n")
  cat("Features with assigned p-value of 1 (constant/insufficient data):", sum(unlist(ttest_results) == 1.0), "\n")
  cat("Total results:", nrow(results_tibble), "\n")
  cat("Group", group1_value, "vs Group", group2_value, "comparison\n")
  cat("First 10 results:\n")
  print(head(results_tibble, 10))
  
  #_Return results
  return(list(
    results = results_tibble,
    n_features = nrow(results_tibble),
    output_file = paste0("Outputs/", output_filename),
    group1_value = group1_value,
    group2_value = group2_value
  ))
}
