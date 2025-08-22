#* Set up targeted and untargeted data for pathway analysis - impute missing untargeted values and join metadata and key with targeted data
  #! Might need to increase vroom connection size to account for large size of untargeted data file
  #+ Read untargeted C18 and HILIC data and merge sets
    amshu_path <- "C://Users//amshu//OneDrive - Emory//Preston, Joshua's files - Amshu Josh PGD//Old Analysis"
    Sys.setenv("VROOM_CONNECTION_SIZE"=2^22)
    C18_unt <- read_csv(paste0(amshu_path, "/C18_untargeted.csv"))
    HILIC_unt <- read_csv(paste0(amshu_path, "/HILIC_untargeted.csv"))
    metadata <- read_csv(paste0(amshu_path, "/metadata.csv"))

  #+ Impute missing C18 values while accounting for PGD status and time associated with the sample ID
    C18_merged_data <- metadata %>%
      left_join(C18_unt, by = c("Sample_ID"= "Sample_Name"))  # Ensure Sample_ID is the matching key

    HILIC_merged_data <- metadata %>%
      left_join(HILIC_unt, by = c("Sample_ID" = "Sample_Name"))  # Ensure Sample_ID is the matching key
    common_cols <- setdiff(intersect(names(HILIC_merged_data), names(C18_merged_data)), "Sample_ID")
    combined_unt <- HILIC_merged_data %>%
      left_join(C18_merged_data %>% select(-all_of(common_cols)), by = "Sample_ID")
    HILIC_numeric_cols <- names(HILIC_merged_data)[7:ncol(HILIC_merged_data)]
    # Compute means grouped by Time and PGD for each numeric column
    HILIC_imp_means <- HILIC_merged_data %>%
      group_by(Time, Clinical_PGD) %>%
      summarise(across(all_of(HILIC_numeric_cols), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"), .groups = "drop")
    # Merge the means back into the dataset for missing rows
    HILIC_unt_imp <- HILIC_merged_data %>%
      left_join(HILIC_imp_means, by = c("Time", "Clinical_PGD")) %>%
      mutate(across(all_of(HILIC_numeric_cols), ~ ifelse(is.na(.x), get(paste0("mean_", cur_column())), .x))) %>%
      select(-starts_with("mean_"))
      # Identify common columns (except the join key)
    common_cols <- setdiff(intersect(names(HILIC_unt_imp), names(C18_unt_imp)), "Sample_ID")
    combined_unt <- HILIC_unt_imp %>%
      left_join(C18_unt_imp %>% select(-all_of(common_cols)), by = "Sample_ID")
    # Impute columns with all NA values
    imp_cols <- combined_unt %>%
      group_by(Patient) %>%
      filter(n_distinct(Time) == 2) %>%
      ungroup()
    pt_imp_id <- unique(imp_cols$Patient)
    for (i in pt_imp_id){
      subset <- imp_cols[imp_cols$Patient == i,]
      imp_time <- c(-12,12,24)[!c(-12,12,24) %in% subset$Time]
      rows_for_imp <- combined_unt[(combined_unt$Clinical_PGD == subset[1,]$Clinical_PGD) & (combined_unt$Time == imp_time),]
      mean_imp <- colMeans(rows_for_imp[,7:ncol(rows_for_imp)], na.rm = TRUE)
      imp_dat <- c(subset[1,1:6],mean_imp)
      imp_dat$Time <- imp_time
      imp_dat$Sample_ID <- paste0(i,"D0")
      combined_unt <- rbind(combined_unt,imp_dat)
    }
  #+ Read targeted data and merge with metadata and metabolite key
    metadata <- read_csv("metadata.csv")
    # Read C18 Raw data, and join metadata
    C18_targeted_raw <- as.tibble(read_csv("C://Users/amshu/OneDrive - Emory/Preston, Joshua's files - Amshu Josh PGD/C18_targeted.csv")) %>%
      filter(str_starts(Sample_ID, "H")) %>%
      mutate(across(where(is.numeric), ~.)) %>%
      select(Sample_ID, where(~ is.numeric(.) && mean(. == 0) <= 0.30)) %>%
      mutate(across(where(is.numeric), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .))) %>%
      # Log2 transform the data
      mutate(across(where(is.numeric), log2)) %>%
      left_join(metadata %>% select(Sample_ID,Clinical_PGD,Time,Patient), by = "Sample_ID") %>%
      select(Patient, Sample_ID, Time, Clinical_PGD,everything()) %>%
      arrange(Clinical_PGD)
    # Read HILIC Raw data, and join metadata
    HILIC_targeted_raw <- as.tibble(read_csv("C://Users/amshu/OneDrive - Emory/Preston, Joshua's files - Amshu Josh PGD/HILIC_targeted.csv")) %>%
      filter(str_starts(Sample_ID, "H")) %>%
      mutate(across(where(is.numeric), ~.)) %>%
      select(Sample_ID, where(~ is.numeric(.) && mean(. == 0) <= 0.30)) %>%
      mutate(across(where(is.numeric), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .))) %>%
      # Log2 transform the data
      mutate(across(where(is.numeric), log2)) %>%
      left_join(metadata %>% select(Sample_ID,Clinical_PGD,Time,Patient), by = "Sample_ID") %>%
      select(Patient, Sample_ID, Time, Clinical_PGD,everything()) %>%
      arrange(Clinical_PGD)
    # Read metabolite key
    feature_key <- read_csv("targeted_feature_key.csv")
    HILIC_key <- feature_key[grepl('HILIC',feature_key$Name),c(1,3)]
    C18_key <- feature_key[grepl('C18',feature_key$Name),c(1,3)]
    # Pivot the data to long format for easier merging
    C18_targeted <- C18_targeted_raw %>%
      pivot_longer(cols = -c(Sample_ID, Patient, Time, Clinical_PGD), names_to = "Metabolite", values_to = "Value")
    HILIC_targeted <- HILIC_targeted_raw %>%
      pivot_longer(cols = -c(Sample_ID, Patient, Time, Clinical_PGD), names_to = "Metabolite", values_to = "Value")
    # Merge the metabolite key with the C18 and HILIC data
    C18_targeted <- C18_targeted %>%
      left_join(C18_key, by = c("Metabolite"="Name")) %>%
      select(Patient, Sample_ID, Time, Clinical_PGD, everything())
    HILIC_targeted <- HILIC_targeted %>%  
      left_join(HILIC_key, by = c("Metabolite"="Name"))
    # Combine C18 and HILIC targeted data
    raw_combined_targeted <- C18_targeted_raw %>%
      left_join(HILIC_targeted_raw, by = c("Sample_ID", "Patient", "Time", "Clinical_PGD")) %>%
      select(Patient, Sample_ID, Time, Clinical_PGD, everything())
    combined_targeted <- rbind(C18_targeted, HILIC_targeted)

#* Figure 1: Pathway analysis and metabolite trends for variations in time only
  #+ Perform linear mixed effects model for each metabolite - Fig 1A and Fig 2A (put results through MetaboAnalyst Functional Analysis to obtain final figures)
    # Log2 transform the data
    combined_unt[,7:ncol(combined_unt)] <- log2(combined_unt[,7:ncol(combined_unt)]+1)
    sig_lme_res <- tibble(mode = character(), m.z = numeric(), r.t = numeric(), p.value = numeric(), comp = character(), met = character())
    # Apply linear mixed effects model for each metabolite
    for (i in colnames(combined_unt)[7:ncol(combined_unt)]) {
      # Create a formula for the linear mixed effects model
      dat <- combined_unt %>% select(all_of(i), Time, Patient, Clinical_PGD)
      tryCatch({model <- lmerTest::lmer(as.formula(paste(i, "~ Time * Clinical_PGD + (1|Patient)")), data = dat)
              res <- summary(model)}, error = function(e) {
        message(paste("Skipping LME for", i, "due to:", e$message))
      })
      # Extract coefficients and p-values
      coeffs <- res$coefficients
      comps <- rownames(res$coefficients)[2:length(res$coefficients[,1])]
      # Append all p-values along with metabolite info to results table
      for (j in 1:length(comps)) {
        mtd <- unlist(strsplit(i, split = "_"))
        sig_lme_res <- rbind(sig_lme_res, tibble(mode = ifelse(mtd[1]=="HILIC","positive","negative"), m.z = as.numeric(mtd[2]), r.t = as.numeric(mtd[3]), p.value = coeffs[j+1,5], comp = comps[j], met = i))
      }
    }
    # Separate results by comparison and sort by p-value
    time_lme <- sig_lme_res[sig_lme_res$comp == "Time",] %>%
      arrange(p.value)
    pgd_lme <- sig_lme_res[sig_lme_res$comp == "Clinical_PGDY",] %>%
      arrange(p.value)
    int_lme <- sig_lme_res[sig_lme_res$comp == "Time:Clinical_PGDY",] %>%
      arrange(p.value)
    write.csv(time_lme[,1:4], file = "time_lme.csv", row.names = FALSE)
    write.csv(pgd_lme[,1:4], file = "pgd_lme.csv", row.names = FALSE)
    write.csv(int_lme[,1:4], file = "int_lme.csv", row.names = FALSE)
    #! Perform MetaboAnalyst functional analysis using each csv file - Functional Analysis not included in this code
  #+ Linear mixed effects on time and time/PGD interaction on targeted data: plot time course of top 5 metabolites based on lowest pvals with fitted splines - Fig 1B, Fig 2B
    # Convert targeted data to z scores
    targ_zscores <- raw_combined_targeted[,5:ncol(raw_combined_targeted)] %>%
      mutate(across(everything(), ~ (.-mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)))
    targ_zscores <- raw_combined_targeted[,1:4] %>%
      cbind(targ_zscores)
    targ_zscores$Clinical_PGD[targ_zscores$Clinical_PGD == "Y'"] <- "Y"
    # Set up linear mixed effects model for each metabolite
    targ_sig_results <- tibble(Metabolite = character(), p.value = numeric(), comp = character())
    for (i in colnames(targ_zscores)[5:ncol(targ_zscores)]) {
      # Create a formula for the linear mixed effects model
      dat <- targ_zscores %>% select(all_of(i), Time, Patient, Clinical_PGD)
      tryCatch({model <- lmerTest::lmer(as.formula(paste(i, "~ Time * Clinical_PGD + (1|Patient)")), data = dat)
              res <- summary(model)}, error = function(e) {
        message(paste("Skipping LME for", i, "due to:", e$message))
      })
      # Extract coefficients and p-values
      coeffs <- res$coefficients
      comps <- rownames(res$coefficients)[2:length(res$coefficients[,1])]
      # Append all p-values along with metabolite info to results table
      for (j in 1:length(comps)) {
        targ_sig_results <- rbind(targ_sig_results, tibble(Metabolite = i, p.value = coeffs[j+1,5], comp = comps[j]))
      }
    }
    # Separate results by comparison and sort by p-value
    time_targ_ <- targ_sig_results[targ_sig_results$comp == "Time",] %>%
      filter(p.value < 0.05) %>%
      arrange(p.value) %>%
    int_targ <- targ_sig_results[targ_sig_results$comp == "Time:Clinical_PGDY",] %>%
      filter(p.value < 0.05) %>%
      arrange(p.value) %>%
    top6_int_targ <- slice_head(int_targ, n = 6)
    top5_time_targ <- slice_head(time_targ, n = 5)

    # Convert targeted z score table to long form and summarize by mean grouped by time and clinical PGD
    long_combined_targ <- targ_zscores %>%
      pivot_longer(cols = -c(Sample_ID, Patient, Time, Clinical_PGD), names_to = "Metabolite", values_to = "Value") %>%
      group_by(Metabolite, Time, Clinical_PGD) %>%
      summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")
    # Obtain top 5 metabolites for time and interaction effects based on p values
    top5_time <- long_combined_targ %>%
      filter(Metabolite %in% top5_time_targ$Metabolite) %>%
      mutate(Time = factor(Time, levels = c("-12", "12", "24"))) %>%
      group_by(Metabolite, Time) %>%
      summarise(mean_value = mean(mean_value, na.rm = TRUE), .groups = "drop")
    top6_int <- long_combined_targ %>%
      filter(Metabolite %in% top6_int_targ$Metabolite) %>%
      mutate(Time = factor(Time, levels = c("-12", "12", "24"))) %>%
      arrange(Metabolite, Time)
    ## Plotting top 5 metabolites for time
    # Convert time to numeric, maintaining evenly spaced positions
    
    time_summary_df <- top5_time %>%
      mutate(Time_numeric = as.numeric(factor(Time, levels = c("-12", "12", "24")))) %>%
      left_join(feature_key, by = c("Metabolite" = "Name"))

    # Nest by metabolite, fit smooth.spline on each group
    time_spline_df <- time_summary_df %>%
      group_by(Identified_Name) %>%
      nest() %>%
      mutate(
        interpolated = map(data, ~ {
          df <- .
          df <- df[order(df$Time_numeric), ]
          time_seq <- seq(min(df$Time_numeric), max(df$Time_numeric), length.out = 101)
          sp <- spline(x = df$Time_numeric, y = df$mean_value, xout = time_seq)
          tibble(Time_numeric = sp$x, Fitted = sp$y)
        })
      ) %>%
      select(Identified_Name, interpolated) %>%
      unnest(interpolated)

    # Plot with color-coded spline curves
    ggplot() +
      geom_point(data = time_summary_df, aes(x = Time_numeric, y = mean_value, color = Identified_Name)) +
      geom_line(data = time_spline_df, aes(x = Time_numeric, y = Fitted, color = Identified_Name), linewidth = 1.2) +
      theme_minimal() +
      labs(
        x = "Time",
        y = "Mean Metabolite Z-Score",
        title = "Spline-Smoothed Time Course for Top Metabolites",
        color = "Metabolite"
      ) +
      scale_color_brewer(palette = "Set1")
  #+ Plot volcano plots for time comparisons
    # Set up data for calculating fold changes
    targ12_PGD <- raw_combined_targeted %>%
      filter(Time == -12, Clinical_PGD == 'Y') %>%
      select(-Sample_ID, -Patient, -Time, -Clinical_PGD)
    targ12_NPGD <- raw_combined_targeted %>%
      filter(Time == 12, Clinical_PGD == 'N') %>%
      select(-Sample_ID, -Patient, -Time, -Clinical_PGD)
    targ24_PGD <- raw_combined_targeted %>%
      filter(Time == 24, Clinical_PGD == 'Y') %>%
      select(-Sample_ID, -Patient, -Time, -Clinical_PGD)
    targ24_NPGD <- raw_combined_targeted %>%
      filter(Time == 24, Clinical_PGD == 'N') %>%
      select(-Sample_ID, -Patient, -Time, -Clinical_PGD)
    # Undo log2 transformation for the targeted data
    undo_log2_transform <- function(df) {
      df <- 2^df # Apply the inverse transform to metabolite columns (4:end)
      return(df)
    }
    targ_list <- list(targ12_PGD, targ12_NPGD, targ24_PGD, targ24_NPGD)
    targ_list <- lapply(targ_list, undo_log2_transform)
    targ_list <- lapply(targ_list, function(df) colMeans(df, na.rm = TRUE))
    # Calculate fold changes and conduct t tests for p-values
    fc_targ12 <- data.frame(Value = log2(targ_list[[1]]/targ_list[[2]]))
    fc_targ12$Time <- "12h"
    fc_targ24 <- data.frame(Value = log2(targ_list[[3]]/targ_list[[4]]))
    fc_targ24$Time <- "24h"
    # Combine fold changes into a single data frame
    fc_combined <- rbind(fc_targ12, fc_targ24)
    # 12 vs -12 t test
    result <- lapply(seq_along(targ12_NPGD), function(i) {
      test <- t.test(targ12_NPGD[[i]], targ12_PGD[[i]])
      as.tibble(data.frame(
        Metabolite = colnames(targ12_PGD)[i],
        p_value = test$p.value,
        mean_difference = diff(test$estimate),
        mean_YPGD = mean(targ12_PGD[[i]]),
        mean_NPGD = mean(targ12_NPGD[[i]])
      ))
    })
    # Combine results into a single table
    targ_12PGD_ttest <- do.call(rbind, result)
    targ_12PGD_ttest$Time <- "12h"
    # 24 vs -12 t test
    result <- lapply(seq_along(targ24_NPGD), function(i) {
      test <- t.test(targ24_NPGD[[i]], targ24_PGD[[i]])
      as.tibble(data.frame(
        Metabolite = colnames(targ24_PGD)[i],
        p_value = test$p.value,
        mean_difference = diff(test$estimate),
        mean_YPGD = mean(targ24_PGD[[i]]),
        mean_NPGD = mean(targ24_NPGD[[i]])
      ))
    })
    # Combine results into a single table
    targ_24PGD_ttest <- do.call(rbind, result)
    targ_24PGD_ttest$Time <- "24h"
    # Combine t test results into a single data frame
    targ_combined_ttest <- rbind(targ_12PGD_ttest, targ_24PGD_ttest)
    # Create a volcano plot function for each comparison
    volc_plot <- function(FC_list, ttest_res, key){
      # Configure FC_table variable to merge fold-change of each metabolite with t-test results
      FC_table <- data.frame(
      Metabolite = unlist(rownames(FC_list)),
      Value = unlist(FC_list$Value),
      Time = unlist(FC_list$Time),
      row.names = NULL
      )
      FC_table$Time <- as.character(FC_table$Time)
      ttest_res$Time <- as.character(ttest_res$Time)
      View(FC_table)
      View(ttest_res)
      # Merge t-test results with fold-change data into the volcano_data variable for plotting
      volcano_data <- merge(FC_table, ttest_res, by = c("Metabolite", "Time"), all = TRUE)
      # View(volcano_data)
      volcano_data <- volcano_data %>% 
        left_join(key, by = c("Metabolite" = "Name"))
      # View(volcano_data)
      # Add a column for color coding
      volcano_data <- volcano_data %>%
        mutate(Legend = case_when(
          p_value < 0.05 & Value > 1 & Time == '24h' ~ "Up in PGD - 24h",
          p_value < 0.05 & Value < -1 & Time == '24h' ~ "Down in PGD - 24h",
          p_value < 0.05 & Value > 1 & Time == '12h' ~ "Up in PGD - 12h",
          p_value < 0.05 & Value < -1 & Time == '12h' ~ "Down in PGD - 12h",
          TRUE ~ "Not Significant"
        ))
      # Add a column for significance, used for labeling significant metabolites
      # volcano_data$Label <- ifelse(volcano_data$Legend != 'Not Significant', 'true', 'false')
      top10sig <- volcano_data %>%
        filter(p_value < 0.05) %>%
        arrange(Value) %>%
        slice_head(n = 10)
      bottom10sig <- volcano_data %>%
        filter(p_value < 0.05) %>%
        arrange(desc(Value)) %>%
        slice_head(n = 10)
      # Combine top and bottom metabolites for labeling
      top_bottom <- rbind(top10sig, bottom10sig)
      # Add a column for labeling metabolites
      volcano_data$Label <- ifelse(volcano_data$Metabolite %in% top_bottom$Metabolite, 'true', 'false')
      # Create the volcano plot
      ggplot(volcano_data, aes(x = Value, y = -log10(p_value), color = Legend)) +
        geom_point(size = 2, alpha = 0.7) +
        scale_color_manual(values = c('Not Significant' = "gray", 'Up in PGD - 12h' = '#800017', 'Down in PGD - 12h' = '#113d6a', 'Up in PGD - 24h' = '#ee516e', 'Down in PGD - 24h' = '#67adf5')) +
        # geom_text_repel(
        #   data = subset(volcano_data, as.logical(Label)),
        #   max.overlaps = 20,
        #   size = 3,
        #   show.legend = FALSE,
        # ) +
        theme_light(base_family = "Arial") +
        labs(x = "Log2 Fold Change", 
            y = "-Log10(p-value)") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
        theme(
          legend.position = "right",
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 13, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          ) +
          guides(color = guide_legend(override.aes = list(shape = 16, size=4)))
    }
    volc_plot(fc_combined, targ_combined_ttest, feature_key) +
      ggtitle("Volcano Plot for PGD vs No PGD at 12h and 24h")
    ggsave(
      filename = "./EACTS Figs/Volcano_combined.svg",
      plot = last_plot(),
      device = 'svg',
      dpi = 600
    )
  #+ Slopes of splines for top metabolites - Fig 1C
    # Create a new data frame for mean values of each metabolite
    time_means_df <- long_combined_targ %>%
      group_by(Metabolite, Time) %>%
      summarise(mean_value = mean(mean_value, na.rm = TRUE), .groups = "drop") %>%
      left_join(feature_key, by = c("Metabolite" = "Name")) %>%
      mutate(Time_numeric = as.numeric(factor(Time, levels = c("-12", "12", "24"))))
    # Calculate slopes at t = 12 for each metabolite
    # Get unique metabolite names
    time_metabolites <- unique(time_means_df$Identified_Name)

    # Store slopes in a data frame
    time_slope_df <- data.frame(Metabolite = character(), Slope = numeric(), stringsAsFactors = FALSE)

    for (met in time_metabolites) {
      sub_data <- subset(time_means_df, Identified_Name == met)
      # Check if we have 3 distinct time points
      if (length(unique(sub_data$Time_numeric)) >= 3) {
        spline_fit <- splinefun(x = sub_data$Time_numeric, y = sub_data$mean_value, method = "natural")
        slope_at_12 <- spline_fit(2, deriv = 1)
        
        time_slope_df <- rbind(time_slope_df, data.frame(Identified_Name = met, Slope = slope_at_12))
      }
    }
    # Sort slopes and select top 10 positive and negative slopes
    time_top_positive_slopes <- time_slope_df %>%
      arrange(desc(Slope)) %>%
      slice_head(n = 10) %>%
      mutate(Slope = round(Slope, 2))
    time_top_negative_slopes <- time_slope_df %>%
      arrange(Slope) %>%
      slice_head(n = 10) %>%
      mutate(Slope = round(Slope, 2))
    # Combine and plot
    time_combined_slopes <- bind_rows(time_top_positive_slopes, time_top_negative_slopes)

    ggplot(time_combined_slopes, aes(x = reorder(Identified_Name, Slope), y = Slope, fill = Slope > 0)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(
        x = "Metabolite",
        y = "Slope at t = 12",
        title = "Top Metabolites by Slope of Z-Score at t = 12",
        fill = "Slope > 0"
      ) +
      scale_fill_manual(values = c("red", "blue")) +
      theme_minimal()
#* Figure 2: Pathway analysis and metabolite trends for variations in time and PGD status
  #+ Plot time course of top 4 metabolites based on interaction effect with fitted splines - Fig 2B
    # Convert time to numeric, maintaining evenly spaced positions
    int_df <- top6_int %>%
      mutate(Time_numeric = case_when(
        Time == -12 ~ 1,
        Time == 12 ~ 2,
        Time == 24 ~ 3
      )) %>%
      left_join(feature_key, by = c("Metabolite" = "Name"))
    # Nest by metabolite and clinical PGD status, fit smooth.spline on each group
    int_spline_df <- int_df %>%
    group_by(Identified_Name, Clinical_PGD) %>%
      nest() %>%
      mutate(
        interpolated = map(data, ~ {
          df <- .
          df <- df[order(df$Time_numeric), ]
          time_seq <- seq(min(df$Time_numeric), max(df$Time_numeric), length.out = 101)
          sp <- spline(x = df$Time_numeric, y = df$mean_value, xout = time_seq)
          tibble(Time_numeric = sp$x, Fitted = sp$y)
        })
      ) %>%
      select(Identified_Name, interpolated) %>%
      unnest(interpolated)
    # Add Time column to spline df for plotting
    int_spline_df <- int_spline_df %>%
      mutate(Time = case_when(
        Time_numeric == 1 ~ "-12",
        Time_numeric == 2 ~ "12",
        Time_numeric == 3 ~ "24"
      )) %>%
      mutate(Time = factor(Time, levels = c("-12", "12", "24")))
    # Use abbreviations for metabolites with exceptionally long names
    int_spline_df$Identified_Name <- ifelse(int_spline_df$Identified_Name == '8-Amino-3,8-dideoxy-D-manno-octulosonate', 'Kdo8N', int_spline_df$Identified_Name)
    int_df$Identified_Name <- ifelse(int_df$Identified_Name == '8-Amino-3,8-dideoxy-D-manno-octulosonate', 'Kdo8N', int_df$Identified_Name)
    # Plot with color-coded spline curves
    ggplot(data = int_spline_df, aes(x = Time_numeric, y = Fitted, color = Clinical_PGD)) +
      geom_point(data = int_df, aes(x = Time_numeric, y = mean_value)) +
      geom_line(linewidth = 1.2, aes(group = Clinical_PGD)) +
      facet_wrap(~Identified_Name, scales = "free_y", ncol = 3) +
      theme_bw(base_family = "Arial") +
      labs(
        x = "Time",
        y = "Metabolite Z-Scores",
        title = "Time Course for Top Metabolites based on Interaction Effect sorted by PGD Status"
      ) +
      scale_color_manual(values = c("Y"= "#800017", "N" = "#113d6a")) +
      theme(
        legend.position = "right",
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      ) + 
      guides(color = guide_legend(title = "Clinical PGD"))
    ggsave(
      filename = "./EACTS Figs/Splines2B_edit.svg",
      plot = last_plot(),
      device = 'svg',
      dpi = 600
    )

#* Untargeted Figures: HCA with heatmaps, PLS-DA, and volcano plots
  # Divide data into groups based on clinical PGD status and time
  combined_unt$Clinical_PGD[combined_unt$Clinical_PGD == "Y'"] <- "Y"
  unt_24 <- combined_unt %>%
    filter(Time == 24) %>%
    select(-Patient, -Time, -Sex, -Age)
  unt_12 <- combined_unt %>%
    filter(Time == 12) %>%
    select(-Patient, -Time, -Sex, -Age)
  unt_YPGD <- combined_unt %>%
    filter(Clinical_PGD == "Y", Time != -12) %>%
    select(-Patient, -Clinical_PGD, -Sex, -Age)
  unt_NPGD <- combined_unt %>%
    filter(Clinical_PGD == "N", Time != -12) %>%
    select(-Patient, -Clinical_PGD, -Sex, -Age)
  # Metabolite values averaged between 12h and 24h
  metcols <- names(unt_24)[1:2]
  geomean <- function(x, na.rm = TRUE) {
    if(na.rm) x <- x[!is.na(x)]
    exp(mean(log(x+1))) - 1
  }
  unt_avg <- bind_rows(unt_12, unt_24) %>%
    group_by(Clinical_PGD) %>%
    summarise(across(-Sample_ID, geomean, na.rm = TRUE), .groups = "drop")
  #+ Heatmap of metabolite clusters
    # Function to create heatmap with hierarchical clustering
    create_heatmap <- function(data, title) {
      # Extract metabolite data and set row names
      df_sub <- data

      # Extract metabolite matrix only
      met_cols <- c("Sample_ID", "Clinical_PGD")
      data_cols <- setdiff(colnames(df_sub), met_cols)
      met_data <- as.matrix(df_sub[, data_cols])
      is.numeric(met_data)
      met_scaled <- as.matrix(scale(met_data))
      row_order <- order(df_sub$Clinical_PGD)
      rownames(met_scaled) <- df_sub$Sample_ID
      met_scaled <- met_scaled[, colSums(is.na(met_scaled)) == 0]
      met_scaled <- met_scaled[row_order, ]
      annotation <- data.frame(Clinical_PGD = df_sub$Clinical_PGD)
      rownames(annotation) <- df_sub$Sample_ID
      annotation <- annotation[row_order, , drop = FALSE]

      # Define colors
      ann_colors <- list(
        Clinical_PGD = c("Y" = "#800017", "N" = "#113d6a")
      )
      pheatmap(
        met_scaled,
        cluster_rows = FALSE,      # hierarchical clustering of samples
        cluster_cols = TRUE,      # hierarchical clustering of metabolites
        annotation_row = annotation,
        annotation_colors = ann_colors,
        scale = "row",            # optional: scale metabolites to mean=0, sd=1
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize_row = 8,
        fontsize_col = 8,
        main = title,
        filename = paste0("./EACTS Figs/Heatmap_", title, ".png"),
      )
    }
    # Create heatmaps for each group
    heatmap_24h <- create_heatmap(unt_24, "PGD vs Non-PGD at 24h")
  #+ PLS-DA analysis
    # Function to perform PLS-DA and plot results
    plot_pls <- function(data, title) {
      # Ensure the data is numeric and remove non-numeric columns
      pls_data <- data %>%
        select(-Sample_ID, -Clinical_PGD) %>%
        mutate(across(everything(), as.numeric)) %>%
        na.omit()  # Remove rows with NA values

      # Perform PLS-DA
      pls_model <- mixOmics::plsda(pls_data, data$Clinical_PGD, ncomp = 2)
      scores <- data.frame(pls_model$variates$X)
      scores$Clinical_PGD <- data$Clinical_PGD
      colnames(scores)[1:2] <- c("Comp1", "Comp2")
      ggplot(scores, 
        aes(x = Comp1, y = Comp2, fill = Clinical_PGD, color = Clinical_PGD)) +
        geom_point(shape = 21, size = 3, color = "black") + # black border points
        stat_ellipse(geom = "polygon", alpha = 0.2, color = NA) + # filled ellipse, no border
        scale_fill_manual(values = c("Y" = "#800017", "N" = "#113d6a")) +
        theme_minimal() +
        labs(title = "PLS-DA of 12h Metabolites by PGD") +
        theme(legend.position = "right"
      )
      ggsave(
        filename = paste0("./EACTS Figs/PLSDA_", title, ".png"),
        plot = last_plot(),
        device = 'png',
        dpi = 600
      )
    }
    pls_12 <- plot_pls(unt_12, "PLS-DA at 12h")
  #+ Volcano Plots
    # Set up data for calculating fold changes
    unt_12_PGD <- unt_12[unt_12$Clinical_PGD == 'Y', 3:ncol(unt_12)]
    unt_12_NPGD <- unt_12[unt_12$Clinical_PGD == 'N', 3:ncol(unt_12)]
    unt_24_PGD <- unt_24[unt_24$Clinical_PGD == 'Y', 3:ncol(unt_24)]
    unt_24_NPGD <- unt_24[unt_24$Clinical_PGD == 'N', 3:ncol(unt_24)]
    # Undo log2 transformation for the untargeted data
    undo_log2_transform <- function(df) {
      df <- 2^df - 1 # Apply the inverse transform to metabolite columns (3:end)
      return(df)
    }
    unt_list <- list(unt_12_PGD, unt_12_NPGD, unt_24_PGD, unt_24_NPGD)
    unt_list <- lapply(unt_list, undo_log2_transform)
    unt_list <- lapply(unt_list, function(df) colMeans(df, na.rm = TRUE))
    # Calculate fold changes and conduct t tests for p-values
    fc_unt12 <- data.frame(Value = log2(unt_list[[1]]/unt_list[[2]]))
    fc_unt12$Time <- "12h"
    fc_unt24 <- data.frame(Value = log2(unt_list[[3]]/unt_list[[4]]))
    fc_unt24$Time <- "24h"
    # Conduct t-tests for each metabolite
    result_12 <- lapply(seq_along(unt_12_NPGD), function(i) {
      test <- t.test(unt_12_NPGD[[i]], unt_12_PGD[[i]])
      as_tibble(data.frame(
        Metabolite = colnames(unt_12_PGD)[i],
        p_value = test$p.value,
        mean_difference = diff(test$estimate),
        mean_YPGD = mean(unt_12_PGD[[i]]),
        mean_NPGD = mean(unt_12_NPGD[[i]])
      ))
    })
    # Combine results into a single table
    unt_12PGD_ttest <- do.call(rbind, result_12)
    unt_12PGD_ttest$Time <- "12h"
    result_24 <- lapply(seq_along(unt_24_NPGD), function(i) {
      test <- t.test(unt_24_NPGD[[i]], unt_24_PGD[[i]])
      as_tibble(data.frame(
        Metabolite = colnames(unt_24_PGD)[i],
        p_value = test$p.value,
        mean_difference = diff(test$estimate),
        mean_YPGD = mean(unt_24_PGD[[i]]),
        mean_NPGD = mean(unt_24_NPGD[[i]])
      ))
    })
    # Combine results into a single table
    unt_24PGD_ttest <- do.call(rbind, result_24)
    unt_24PGD_ttest$Time <- "24h"

    volc_PGD <- function(FC_table, ttest_res, timechar){
      # Configure FC_table variable to merge fold-change of each metabolite with t-test results
      FC_table <- data.frame(
        Metabolite = unlist(rownames(FC_table)),
        Value = unlist(FC_table$Value),
        Time = unlist(FC_table$Time),
        row.names = NULL
      )
      FC_table$Time <- as.character(FC_table$Time)
      ttest_res$Time <- as.character(ttest_res$Time)
      # FC_list <- list()
      # for (metabolite in colnames(baseline)[3:ncol(baseline)]) {
      #   fc <- log2(mean(comparison[, metabolite]) / mean(baseline[, metabolite]))
      #   FC_list[[metabolite]] <- data.frame(Value = fc)
      # }
      # View(FC_list)
      # FC_table <- data.frame(
      # Metabolite = unlist(rownames(FC_list)),
      # Value = unlist(FC_list$Value),
      # row.names = NULL
      # )
      # Merge t-test results with fold-change data into the volcano_data variable for plotting
      volcano_data <- merge(FC_table, ttest_res, by = "Metabolite", all = TRUE)
        View(volcano_data)
        # View(volcano_data)
        # Add a column for color coding
        volcano_data <- volcano_data %>%
          mutate(Legend = case_when(
            p_value < 0.05 & Value > 1 ~ "Up in PGD",
            p_value < 0.05 & Value < -1 ~ "Down in PGD",
            TRUE ~ "Not Significant"
          ))
        # Add a column for significance, used for labeling significant metabolites
        # volcano_data$Label <- ifelse(volcano_data$Legend != 'Not Significant', 'true', 'false')
      top10sig <- volcano_data %>%
        filter(p_value < 0.05) %>%
        arrange(Value) %>%
        slice_head(n = 10)
      bottom10sig <- volcano_data %>%
        filter(p_value < 0.05) %>%
        arrange(desc(Value)) %>%
        slice_head(n = 10)
      # Combine top and bottom metabolites for labeling
      top_bottom <- rbind(top10sig, bottom10sig)
      # Add a column for labeling metabolites
      volcano_data$Label <- ifelse(volcano_data$Metabolite %in% top_bottom$Metabolite, 'true', 'false')
      # Create the volcano plot
      ggplot(volcano_data, aes(x = Value, y = -log10(p_value), color = Legend)) +
        geom_point(size = 2, alpha = 0.7) +
        scale_color_manual(values = c('Not Significant' = "gray", 'Up in PGD' = '#800017', 'Down in PGD' = '#113d6a')) +
        # geom_text_repel(
        #   data = subset(volcano_data, as.logical(Label)),
        #   max.overlaps = 20,
        #   size = 3,
        #   show.legend = FALSE,
        # ) +
        theme_light(base_family = "Arial") +
        labs(x = "Log2 Fold Change", 
            y = "-Log10(p-value)") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
        theme(
          legend.position = "right",
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 13, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          ) +
          guides(color = guide_legend(override.aes = list(shape = 16, size=4)))
    }
    volc_PGD(fc_unt12, unt_12PGD_ttest, "12h") +
      ggtitle("Volcano Plot for PGD vs Non-PGD at 12h")
    ggsave(
      filename = "./EACTS Figs/Volcano_PGD_vs_NonPGD_12h.png",
      plot = last_plot(),
      device = 'png',
      dpi = 600
    )
