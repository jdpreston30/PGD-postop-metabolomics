#* 7: Numbers for results section
#+ 7.1: Untargeted paragraph
#- 7.1.1: Count features in complete UFT
hilic_count <- UFT |>
  select(starts_with("HILIC")) |>
  ncol()  
c18_count <- UFT |>
  select(starts_with("C18")) |>
  ncol()
combined_count <- hilic_count + c18_count
#- 7.1.2: Remaining features following elimination of 20% detection
UFT_filtered_count <- UFT_filtered |>
  select(starts_with("C18"), starts_with("HILIC")) |>
  ncol()
#- 7.1.3: LIMMA significant
sig_limma_untarg <- tibble(
  contrast = c("interaction", "time_main", "group_main"),
  n_sig = c(
    sum(limma_untarg$interaction$p.value < 0.05, na.rm = TRUE),
    sum(limma_untarg$time_main$p.value < 0.05, na.rm = TRUE),
    sum(limma_untarg$group_main$p.value < 0.05, na.rm = TRUE)
  )
)
#- 7.1.4: Permutation test results
perm_12h_p <- plsda_12h$permutation$p_value
perm_24h_p <- plsda_24h$permutation$p_value
perm_combined_p <- plsda_combined$permutation$p_value
#+ 7.2: Targeted paragraph
#- 7.2.1: Count total features in TFT_combined by ion mode
TFT_hilic_count <- TFT_combined |>
  select(starts_with("HILIC")) |>
  ncol()
TFT_c18_count <- TFT_combined |>
  select(starts_with("C18")) |>
  ncol()
TFT_combined_count <- TFT_hilic_count + TFT_c18_count
#- 7.2.2: Count Level 1 (library confirmed) and Level 3 (annotated) features in final dataset
# Get feature names from TFT_combined
TFT_feature_names <- names(TFT_combined)[!names(TFT_combined) %in% c("Patient", "Sample", "Sample_ID", "Time", "severe_PGD", "PGD_grade_tier", "any_PGD")]
# Count how many are Level 1 vs Level 3 by matching with TFT_merged_features
TFT_level1_count <- TFT_merged_features |>
  filter(feature %in% TFT_feature_names, lib_conf == "Y") |>
  nrow()
TFT_level3_count <- TFT_merged_features |>
  filter(feature %in% TFT_feature_names, lib_conf == "N") |>
  nrow()
#+ 7.3: Results Summary and Exports
cat("\n", crayon::bgBlue(crayon::white("\n === RESULTS SUMMARY ===\n")), "\n")
cat(  
  crayon::blue(sprintf("• A total of %d untargeted features were detected in the analyzed samples (%d HILIC+ and %d C18-).\n", combined_count, hilic_count, c18_count)),
  crayon::yellow(sprintf("• In contrast, partial least squares discriminant analysis (PLS-DA) revealed separation between severe PGD and no severe PGD at 12 (Figure 1B) and 24 hours post-transplant (Figure 1C), though permutation testing only demonstrated statistically robust model performance at 24 hours (p=%.3f) compared to 12 hours (p=%.3f). A PLS-DA of combined timepoints did result in a statistically robust model (p=%.3f).\n", perm_24h_p, perm_12h_p, perm_combined_p)),
  crayon::green(sprintf("• To transition from pathway-level findings to metabolite-level patterns, we analyzed a subset of identified or annotated metabolites. This encompassed %d features (%d HILIC+ and %d C18-), with %d Level 1 identifications and %d Level 3 annotations.\n", TFT_combined_count, TFT_hilic_count, TFT_c18_count, TFT_level1_count, TFT_level3_count)))

