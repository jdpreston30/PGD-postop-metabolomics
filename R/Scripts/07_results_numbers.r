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
UFT_filtered_count <- UFT |>
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
# #+ 7.2: Targeted paragraph
# #- 7.2.1: Count features in complete TFT
# HILIC_count_targ <- HILIC_targeted_raw |>
#   select(-Sample_ID) |>
#   ncol()
# C18_count_targ <- C18_targeted_raw |>
#   select(-Sample_ID) |>
#   ncol()
# combined_count_targ <- HILIC_count_targ + C18_count_targ
# #- 7.2.2: Remaining features following elimination of 20% detection
# targ_post_elim <- TFT |>
#   select(-meta_cols) |>
#   ncol()
# #- 7.2.3: Detection within both modes based on KEGG
# both_modes <- feature_key |>
#   group_by(KEGGID) |>
#   summarise(n_modes = n_distinct(ion_mode), .groups = "drop") |>
#   filter(n_modes == 2)
# n_both_modes <- nrow(both_modes)
# #- 7.2.4: Number of isomers
# collapsed_all_isomers_counts <- collapsed_all_isomers |>
#   filter(n_candidates > 1) |>
#   select(Mode, mz, time) |>
#   unique() |>
#   nrow()
# #- 7.2.5: Volcano counts (Up)
# balloon_data_12_up <- balloon_data_12 |>
#   filter(Legend == "Up in PGD") |>
#   nrow()
# balloon_data_24_up <- balloon_data_24 |>
#   filter(Legend == "Up in PGD") |>
#   nrow()
# #- 7.2.6: Volcano counts (Down)
# balloon_data_12_down <- balloon_data_12 |>
#   filter(Legend == "Down in PGD") |>
#   nrow()
# balloon_data_24_down <- balloon_data_24 |>
#   filter(Legend == "Down in PGD") |>
#   nrow()
# #- 7.2.7: LIMMA significant
# sig_limma_targ <- tibble(
#   contrast = c("interaction", "time_main", "group_main"),
#   n_sig = c(
#     sum(limma_targ$interaction$p.value < 0.05, na.rm = TRUE),
#     sum(limma_targ$time_main$p.value < 0.05, na.rm = TRUE),
#     sum(limma_targ$group_main$p.value < 0.05, na.rm = TRUE)
#   )
# )
#+ 7.3: Results Summary and Exports
#- 7.3.1: Print colored results summary as bulleted sentences
cat("\n", crayon::bgBlue(crayon::white("\n === RESULTS SUMMARY ===\n")), "\n")
cat(crayon::blue(sprintf("• A total of %d untargeted features were detected in the analyzed samples (%d HILIC+ and %d C18−). Following elimination of features detected in <20%% of samples, %d features remained for analysis.\n\n", combined_count, hilic_count, c18_count, UFT_filtered_count)))
cat(crayon::cyan(sprintf("• A total of %d targeted features were quantified (%d HILIC and %d C18−). Following elimination of features with >20%% missing values, %d features remained for analysis.\n\n", combined_count_targ, HILIC_count_targ, C18_count_targ, targ_post_elim)))
cat(crayon::green(sprintf("• Of the %d targeted metabolites identified via accurate mass, %d were detected in both LC modes (HILIC and C18−).\n\n", targ_post_elim, n_both_modes)))
cat(crayon::yellow(sprintf("• Structural elucidation revealed %d putative isomeric features based on identical m/z within the same analytical mode.\n\n", collapsed_all_isomers_counts)))
cat(crayon::magenta(sprintf("• At 12 hours post-transplant, %d metabolites were significantly increased and %d were significantly decreased in patients with severe PGD (nominal p < 0.05).\n\n", balloon_data_12_up, balloon_data_12_down)))
cat(crayon::red(sprintf("• At 24 hours post-transplant, %d metabolites were significantly increased and %d were significantly decreased in patients with severe PGD (nominal p < 0.05).\n\n", balloon_data_24_up, balloon_data_24_down)))
cat(crayon::bold("\nStatistical validation (LIMMA, nominal p < 0.05):\n"))
cat(crayon::bgGreen(crayon::black(sprintf("   Untargeted: interaction effect %d, time-main %d, group-main %d\n", sig_limma_untarg$n_sig[1], sig_limma_untarg$n_sig[2], sig_limma_untarg$n_sig[3]))))
cat(crayon::bgMagenta(crayon::black(sprintf("   Targeted: interaction effect %d, time-main %d, group-main %d\n\n", sig_limma_targ$n_sig[1], sig_limma_targ$n_sig[2], sig_limma_targ$n_sig[3]))))


