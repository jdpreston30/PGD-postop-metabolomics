#8 7.1: Targeted volcano analysis
#+ 7.1.1: Set up data for calculating fold changes
    #- 7.1.1.1: Prepare data subsets
    targ12_PGD <- combined_TFT %>%
      filter(Time == 12, Clinical_PGD == 'Y') %>%
      select(-Sample_ID, -Patient, -Time, -Clinical_PGD, -Sex, -Age)
    targ12_NPGD <- combined_TFT %>%
      filter(Time == 12, Clinical_PGD == 'N') %>%
      select(-Sample_ID, -Patient, -Time, -Clinical_PGD, -Sex, -Age)
    targ24_PGD <- combined_TFT %>%
      filter(Time == 24, Clinical_PGD == 'Y') %>%
      select(-Sample_ID, -Patient, -Time, -Clinical_PGD, -Sex, -Age)
    targ24_NPGD <- combined_TFT %>%
      filter(Time == 24, Clinical_PGD == 'N') %>%
      select(-Sample_ID, -Patient, -Time, -Clinical_PGD, -Sex, -Age)
    targcomb_PGD <- combined_TFT %>%
      filter(Clinical_PGD == 'Y') %>%
      select(-Sample_ID, -Patient, -Time, -Clinical_PGD, -Sex, -Age)
    targcomb_NPGD <- combined_TFT %>%
      filter(Clinical_PGD == 'N') %>%
      select(-Sample_ID, -Patient, -Time, -Clinical_PGD, -Sex, -Age)
#+ 7.1.2: Undo log2 transformation for the targeted data
    #- 7.1.2.1: Apply undo_log2_transform and calculate means
    targ_list <- list(targ12_PGD, targ12_NPGD, targ24_PGD, targ24_NPGD, targcomb_PGD, targcomb_NPGD)
    targ_list <- lapply(targ_list, undo_log2_transform)
    targ_list <- lapply(targ_list, function(df) colMeans(df, na.rm = TRUE))
#+ 7.1.3: Calculate fold changes and conduct t tests for p-values
    #- 7.1.3.1: Calculate fold changes
    fc_targ12 <- data.frame(Value = log2(targ_list[[1]]/targ_list[[2]]))
    fc_targ12$Time <- "12h"
    fc_targ24 <- data.frame(Value = log2(targ_list[[3]]/targ_list[[4]]))
    fc_targ24$Time <- "24h"
    fc_targcomb <- data.frame(Value = log2(targ_list[[5]]/targ_list[[6]]))
    fc_targcomb$Time <- "12h and 24h"
#+ 7.1.4: 12h t test
    #- 7.1.4.1: Run t tests for 12h
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
#+ 7.1.5: Combine results into a single table (12h)
    #- 7.1.5.1: Bind 12h t test results
    targ_12PGD_ttest <- do.call(rbind, result)
    targ_12PGD_ttest$Time <- "12h"
#+ 7.1.6: 24h t test
    #- 7.1.6.1: Run t tests for 24h
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
#+ 7.1.7: Combine results into a single table (24h)
    #- 7.1.7.1: Bind 24h t test results
    targ_24PGD_ttest <- do.call(rbind, result)
    targ_24PGD_ttest$Time <- "24h"
#+ 7.1.8: Combined time t-test
    #- 7.1.8.1: Run t tests for combined time
    result <- lapply(seq_along(targcomb_NPGD), function(i) {
      test <- t.test(targcomb_NPGD[[i]], targcomb_PGD[[i]])
      as.tibble(data.frame(
        Metabolite = colnames(targcomb_PGD)[i],
        p_value = test$p.value,
        mean_difference = diff(test$estimate),
        mean_YPGD = mean(targcomb_PGD[[i]]),
        mean_NPGD = mean(targcomb_NPGD[[i]])
      ))
    })
#+ 7.1.9: Combine results into a single table (combined)
    #- 7.1.9.1: Bind combined t test results
    targ_comb_ttest <- do.call(rbind, result)
    targ_comb_ttest$Time <- "12h and 24h"
#+ 7.1.10: Create a volcano plot function for each comparison
    #- 7.1.10.1: Volcano plot for 12h
    volc_12 <- volc_targ(fc_targ12, targ_12PGD_ttest, feature_key)
    #- 7.1.10.2: Balloon plot for 12h
    balloon_12 <- balloon_plot(volc_12)
    View(balloon_12$ball_plot)
    ggsave(
      filename = "./Bubble_12h.png",
      plot = last_plot(),
      device = 'png',
      dpi = 600
    )
    # Combined time t-test
    result <- lapply(seq_along(targcomb_NPGD), function(i) {
      test <- t.test(targcomb_NPGD[[i]], targcomb_PGD[[i]])
      as.tibble(data.frame(
        Metabolite = colnames(targcomb_PGD)[i],
        p_value = test$p.value,
        mean_difference = diff(test$estimate),
        mean_YPGD = mean(targcomb_PGD[[i]]),
        mean_NPGD = mean(targcomb_NPGD[[i]])
      ))
    })
    # Combine results into a single table
    targ_comb_ttest <- do.call(rbind, result)
    targ_comb_ttest$Time <- "12h and 24h"
                        
    # Create a volcano plot function for each comparison
    volc_12 <- volc_targ(fc_targ12, targ_12PGD_ttest, feature_key)
    balloon_12 <- balloon_plot(volc_12)
    View(balloon_12$ball_plot)
    ggsave(
      filename = "./Bubble_12h.png",
      plot = last_plot(),
      device = 'png',
      dpi = 600
    )

#+ 7.2: Untargeted volcano analysis
#+ 7.2.1: Set up data for calculating fold changes
    #- 7.2.1.1: Prepare data subsets
    unt_12_PGD <- unt_12[unt_12$Clinical_PGD == 'Y', 3:ncol(unt_12)]
    unt_12_NPGD <- unt_12[unt_12$Clinical_PGD == 'N', 3:ncol(unt_12)]
    unt_24_PGD <- unt_24[unt_24$Clinical_PGD == 'Y', 3:ncol(unt_24)]
    unt_24_NPGD <- unt_24[unt_24$Clinical_PGD == 'N', 3:ncol(unt_24)]
    unt_YPGD <- UFT_YPGD[, 3:ncol(UFT_YPGD)]
    unt_NPGD <- UFT_NPGD[, 3:ncol(UFT_NPGD)]
#+ 7.2.2: Undo log2 transformation for the untargeted data
    #- 7.2.2.1: Apply undo_log2_transform and calculate means
    unt_list <- list(unt_12_PGD, unt_12_NPGD, unt_24_PGD, unt_24_NPGD, unt_YPGD, unt_NPGD)
    unt_list <- lapply(unt_list, undo_log2_transform)
    unt_list <- lapply(unt_list, function(df) colMeans(df, na.rm = TRUE))
#+ 7.2.3: Calculate fold changes
    #- 7.2.3.1: Calculate fold changes
    fc_unt12 <- data.frame(Value = log2(unt_list[[1]]/unt_list[[2]]))
    fc_unt12$Time <- "12h"
    fc_unt24 <- data.frame(Value = log2(unt_list[[3]]/unt_list[[4]]))
    fc_unt24$Time <- "24h"
    fc_comb <- data.frame(Value = log2(unt_list[[5]]/unt_list[[6]]))
    fc_comb$Time <- "12h and 24h"
#+ 7.2.4: Conduct t-tests for each metabolite
    #- 7.2.4.1: 12h t test
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
#+ 7.2.5: Combine results into a single table (12h)
    #- 7.2.5.1: Bind 12h t test results
    unt_12PGD_ttest <- do.call(rbind, result_12)
    unt_12PGD_ttest$Time <- "12h"
    #- 7.2.4.2: 24h t test
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
#+ 7.2.6: Combine results into a single table (24h)
    #- 7.2.6.1: Bind 24h t test results
    unt_24PGD_ttest <- do.call(rbind, result_24)
    unt_24PGD_ttest$Time <- "24h"
    #- 7.2.4.3: Combined time t test
    result_comb <- lapply(seq_along(unt_NPGD), function(i) {
      test <- t.test(unt_NPGD[[i]], unt_YPGD[[i]])
      as_tibble(data.frame(
        Metabolite = colnames(unt_YPGD)[i],
        p_value = test$p.value,
        mean_difference = diff(test$estimate),
        mean_YPGD = mean(unt_YPGD[[i]]),
        mean_NPGD = mean(unt_NPGD[[i]])
      ))
    })
    unt_comb_ttest <- do.call(rbind, result_comb)
    unt_comb_ttest$Time <- "12h and 24h"
#+ 7.2.7: Create volcano plot for untargeted data
    #- 7.2.7.1: Volcano plot for 24h
    volc_24 <- volc_unt(fc_unt24, unt_24PGD_ttest, "24h") +
      ggtitle("Volcano Plot for PGD vs Non-PGD at 24h")
#+ 7.3: Balloon plot for untargeted
  #- 7.3.1: Generate balloon plot - stored in balloon$ball_plot
  balloon <- balloon_plot(volc_24)
