#* 5: Targeted Balloon and Volcano Plots
#+ 5.1: Set up data for calculating fold changes
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
#+ 5.2: Undo log2 transformation for the targeted data
  targ_list <- list(targ12_PGD, targ12_NPGD, targ24_PGD, targ24_NPGD, targcomb_PGD, targcomb_NPGD)
  targ_list <- lapply(targ_list, undo_log2_transform)
  targ_list <- lapply(targ_list, function(df) colMeans(df, na.rm = TRUE))
#+ 5.3: Calculate fold changes and conduct t tests for p-values
  fc_targ12 <- data.frame(Value = log2(targ_list[[1]]/targ_list[[2]]))
  fc_targ12$Time <- "12h"
  fc_targ24 <- data.frame(Value = log2(targ_list[[3]]/targ_list[[4]]))
  fc_targ24$Time <- "24h"
  fc_targcomb <- data.frame(Value = log2(targ_list[[5]]/targ_list[[6]]))
  fc_targcomb$Time <- "12h and 24h"
#+ 5.4: 12h t test
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
#+ 5.5: Combine results into a single table (12h)
  targ_12PGD_ttest <- do.call(rbind, result)
  targ_12PGD_ttest$Time <- "12h"
#+ 5.6: 24h t test
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
#+ 5.7: Combine results into a single table (24h)
  targ_24PGD_ttest <- do.call(rbind, result)
  targ_24PGD_ttest$Time <- "24h"
#+ 5.8: Combined time t-test
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
#+ 5.9: Combine results into a single table (combined)
  targ_comb_ttest <- do.call(rbind, result)
  targ_comb_ttest$Time <- "12h and 24h"                      
#+ 5.10: Create a volcano plot function for each comparison
  volc_12 <- volc_targ(fc_targ12, targ_12PGD_ttest, feature_key)
#+ 5.11: Create a balloon plot based on volcano plot results
  balloon_12 <- balloon_plot(volc_12$volcano_data)