#* 5: Subject-based analysis and limma workflow
  #+ 1.1: Fix imputed subjects to set them back to NAs
    #- 1.1.1: Create NA rows for missing samples
      new_rows <- tibble::tribble(
        ~Patient, ~Time, ~Clinical_PGD,
        "H2", 12, "Y",
        "H5", 12, "Y",
        "H19", 24, "Y"
      ) %>%
        tibble::add_column(!!!setNames(rep(list(NA_real_), length(met_cols)), met_cols)) %>%
        mutate(across(c(Patient, Time, Clinical_PGD), as.factor))
    #- 1.1.2: Drop any existing rows for those combos
      combined_TFT_LMM_i <- combined_TFT %>%
        select(-c(Sex, Age, Sample_ID)) %>%
        filter(
          !(Patient == "H2" & Time == 12 & Clinical_PGD == "Y"),
          !(Patient == "H5" & Time == 12 & Clinical_PGD == "Y"),
          !(Patient == "H19" & Time == 24 & Clinical_PGD == "Y")
        )
    #- 1.1.3: Bind them back
      combined_TFT_LMM <- bind_rows(combined_TFT_LMM_i, new_rows) %>%
        arrange(Patient, Time)
  #+ 1.2: limma orientation and full workflow
    #- 1.2.1: Ensure factors have the right baseline levels
      combined_TFT_LMM$Time <- factor(combined_TFT_LMM$Time, levels = c("12", "24"))
      combined_TFT_LMM$Clinical_PGD <- factor(combined_TFT_LMM$Clinical_PGD, levels = c("N", "Y"))
    #- 1.2.2: Identify metabolite columns (4+)
      met_cols <- names(combined_TFT_LMM)[4:ncol(combined_TFT_LMM)]
    #- 1.2.3: Build expression matrix as features x samples (rows = metabolites, cols = samples)
      expr <- t(as.matrix(combined_TFT_LMM[, met_cols, drop = FALSE]))
      colnames(expr) <- paste0(combined_TFT_LMM$Patient, "_", combined_TFT_LMM$Time)
    #- 1.2.4: Design matrix (rows must match number of columns in expr == nrow(combined_TFT_LMM))
      design <- model.matrix(~ Time * Clinical_PGD, data = combined_TFT_LMM)
    #- 1.2.5: Blocking factor for repeated measures (length must equal ncol(expr))
      block <- combined_TFT_LMM$Patient
    #- 1.2.6: Estimate within-subject correlation and fit
      corfit <- limma::duplicateCorrelation(expr, design, block = block)
      fit <- limma::lmFit(expr, design, block = block, correlation = corfit$consensus)
      fit <- limma::eBayes(fit)
    #- 1.2.7: Extract Group x Time interaction results (flat vs non-flat between groups)
      tt_inter <- limma::topTable(fit, coef = "Time24:Clinical_PGDY", number = Inf, sort.by = "P")
    #- 1.2.8: Tidy tibble with metabolite, estimate, p, FDR
      interaction_results <- tibble::tibble(
        Metabolite = rownames(tt_inter),
        logFC = tt_inter[, "logFC"],
        t = tt_inter[, "t"],
        P = tt_inter[, "P.Value"],
        FDR = tt_inter[, "adj.P.Val"]
      ) %>%
        dplyr::arrange(FDR) %>%
        left_join(combined_key, by = "Metabolite") %>%
        filter(P < 0.05)
  write.csv(interaction_results, "interaction_results.csv", row.names = FALSE)
