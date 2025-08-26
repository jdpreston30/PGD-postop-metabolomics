#* Untargeted volcano plots
    #+ Set up data for calculating fold changes
    unt_12_PGD <- unt_12[unt_12$Clinical_PGD == 'Y', 3:ncol(unt_12)]
    unt_12_NPGD <- unt_12[unt_12$Clinical_PGD == 'N', 3:ncol(unt_12)]
    unt_24_PGD <- unt_24[unt_24$Clinical_PGD == 'Y', 3:ncol(unt_24)]
    unt_24_NPGD <- unt_24[unt_24$Clinical_PGD == 'N', 3:ncol(unt_24)]
    unt_YPGD <- UFT_YPGD[, 3:ncol(UFT_YPGD)]
    unt_NPGD <- UFT_NPGD[, 3:ncol(UFT_NPGD)]
    #+ Undo log2 transformation for the untargeted data
    unt_list <- list(unt_12_PGD, unt_12_NPGD, unt_24_PGD, unt_24_NPGD, unt_YPGD, unt_NPGD)
    unt_list <- lapply(unt_list, undo_log2_transform)
    unt_list <- lapply(unt_list, function(df) colMeans(df, na.rm = TRUE))
    #+ Calculate fold changes
    fc_unt12 <- data.frame(Value = log2(unt_list[[1]]/unt_list[[2]]))
    fc_unt12$Time <- "12h"
    fc_unt24 <- data.frame(Value = log2(unt_list[[3]]/unt_list[[4]]))
    fc_unt24$Time <- "24h"
    fc_comb <- data.frame(Value = log2(unt_list[[5]]/unt_list[[6]]))
    fc_comb$Time <- "12h and 24h"
    #+ Conduct t-tests for each metabolite
    #+ 12h t test
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
    #+ Combine results into a single table
    unt_12PGD_ttest <- do.call(rbind, result_12)
    unt_12PGD_ttest$Time <- "12h"
    #+ 24h t test
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
    #+ Combine results into a single table
    unt_24PGD_ttest <- do.call(rbind, result_24)
    unt_24PGD_ttest$Time <- "24h"
    #+ Combined time t test
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
    #+ Create volcano plot for untargeted data
    volc_24 <- volc_unt(fc_unt24, unt_24PGD_ttest, "24h") +
      ggtitle("Volcano Plot for PGD vs Non-PGD at 24h")