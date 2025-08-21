#' Perform statistical comparisons on targeted metabolomics data
#'
#' @param data Data frame with metabolite data, must include Patient_ID and grouping variable
#' @param grouping_var Character string specifying the column name for grouping (e.g., "Clade", "Variant", "T")
#' @param test_type Character string specifying test type: "t_test" or "anova"
#' @param exclude_cols Character vector of column names to exclude from analysis (default: c("Patient_ID"))
#' @param exclude_isotopes Logical, whether to exclude isotope metabolites (15N, 13C) from analysis (default: TRUE)
#' @return Tibble with statistical test results including p-values, FDR correction, and group means
#' @export
targeted_metabolite_comparison <- function(data, grouping_var, test_type = "t_test", exclude_cols = c("Patient_ID"), exclude_isotopes = TRUE) {
  
  # Load required libraries
  library(dplyr)
  library(purrr)
  library(broom)
  
  # Identify metabolite columns
  metabolite_cols <- names(data)[!names(data) %in% c(exclude_cols, grouping_var)]
  
  # Filter out isotopes if requested
  if (exclude_isotopes) {
    metabolite_cols <- metabolite_cols[!str_detect(metabolite_cols, "15N|13C")]
  }
  
  # Check if grouping variable exists
  if (!grouping_var %in% names(data)) {
    stop(paste("Grouping variable", grouping_var, "not found in data"))
  }
  
  # Perform statistical tests for each metabolite
  test_results <- map_dfr(metabolite_cols, ~{
    # Extract data for current metabolite
    metabolite_data <- data %>%
      filter(!is.na(.data[[grouping_var]]), !is.na(.data[[.x]])) %>%
      select(group = all_of(grouping_var), metabolite = all_of(.x))
    
    # Convert group to factor
    metabolite_data$group <- as.factor(metabolite_data$group)
    
    # Check if we have sufficient data
    if(nrow(metabolite_data) > 0 && 
       length(unique(metabolite_data$group)) >= 2 &&
       all(table(metabolite_data$group) >= 2)) {
      
      if(test_type == "t_test") {
        # Perform t-test (only works for 2 groups)
        if(length(unique(metabolite_data$group)) == 2) {
          test_result <- t.test(metabolite ~ group, data = metabolite_data)
          
          # Calculate group means
          group_means <- metabolite_data %>%
            group_by(group) %>%
            summarise(mean_val = mean(metabolite, na.rm = TRUE), .groups = "drop")
          
          # Create named list of means
          means_list <- setNames(group_means$mean_val, paste0("mean_", group_means$group))
          
          # Calculate fold change (difference for log2 data)
          if(nrow(group_means) == 2) {
            fold_change <- group_means$mean_val[2] - group_means$mean_val[1]
          } else {
            fold_change <- NA
          }
          
          result <- tibble(
            Metabolite = .x,
            test_statistic = test_result$statistic,
            p_value = test_result$p.value,
            fold_change = fold_change
          )
          
          # Add group means as separate columns
          for(i in seq_along(means_list)) {
            result[[names(means_list)[i]]] <- means_list[[i]]
          }
          
        } else {
          # More than 2 groups - cannot perform t-test
          result <- tibble(
            Metabolite = .x,
            test_statistic = NA,
            p_value = NA,
            fold_change = NA,
            error = "t-test requires exactly 2 groups"
          )
        }
        
      } else if(test_type == "anova") {
        # Perform ANOVA
        test_result <- aov(metabolite ~ group, data = metabolite_data)
        anova_summary <- summary(test_result)
        
        # Calculate group means
        group_means <- metabolite_data %>%
          group_by(group) %>%
          summarise(mean_val = mean(metabolite, na.rm = TRUE), .groups = "drop")
        
        # Create named list of means
        means_list <- setNames(group_means$mean_val, paste0("mean_", group_means$group))
        
        # Calculate overall range (max - min means)
        fold_change <- max(group_means$mean_val) - min(group_means$mean_val)
        
        result <- tibble(
          Metabolite = .x,
          test_statistic = anova_summary[[1]]$`F value`[1],
          p_value = anova_summary[[1]]$`Pr(>F)`[1],
          fold_change = fold_change
        )
        
        # Add group means as separate columns
        for(i in seq_along(means_list)) {
          result[[names(means_list)[i]]] <- means_list[[i]]
        }
        
      } else {
        stop("test_type must be either 't_test' or 'anova'")
      }
      
    } else {
      # Insufficient data
      result <- tibble(
        Metabolite = .x,
        test_statistic = NA,
        p_value = NA,
        fold_change = NA,
        error = "Insufficient data"
      )
    }
    
    return(result)
  })
  
  # Add FDR correction and significance flags
  final_results <- test_results %>%
    mutate(
      p_value_fdr = p.adjust(p_value, method = "fdr"),
      significant = p_value < 0.05 & !is.na(p_value),
      significant_fdr = p_value_fdr < 0.05 & !is.na(p_value_fdr)
    ) %>%
    arrange(p_value_fdr)
  
  # Report results
  n_significant <- sum(final_results$significant, na.rm = TRUE)
  n_significant_fdr <- sum(final_results$significant_fdr, na.rm = TRUE)
  
  cat("Statistical comparison results:\n")
  cat("Test type:", test_type, "\n")
  cat("Grouping variable:", grouping_var, "\n")
  cat("Isotopes excluded:", exclude_isotopes, "\n")
  cat("Number of metabolites tested:", nrow(final_results), "\n")
  cat("Significantly different (uncorrected p < 0.05):", n_significant, "\n")
  cat("Significantly different (FDR corrected p < 0.05):", n_significant_fdr, "\n\n")
  
  return(final_results)
}
