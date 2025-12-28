# #+ 7.2: Targeted paragraph

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



  crayon::cyan(sprintf("• A total of %d targeted features were quantified (%d HILIC and %d C18−). Following elimination of features with >20%% missing values, %d features remained for analysis.\n", combined_count_targ, HILIC_count_targ, C18_count_targ, targ_post_elim)),
  crayon::green(sprintf("• Of the %d targeted metabolites identified via accurate mass, %d were detected in both LC modes (HILIC and C18−).\n", targ_post_elim, n_both_modes)),
  crayon::yellow(sprintf("• Structural elucidation revealed %d putative isomeric features based on identical m/z within the same analytical mode.\n", collapsed_all_isomers_counts)),
  crayon::magenta(sprintf("• At 12 hours post-transplant, %d metabolites were significantly increased and %d were significantly decreased in patients with severe PGD (nominal p < 0.05).\n", balloon_data_12_up, balloon_data_12_down)),
  crayon::red(sprintf("• At 24 hours post-transplant, %d metabolites were significantly increased and %d were significantly decreased in patients with severe PGD (nominal p < 0.05).\n", balloon_data_24_up, balloon_data_24_down)),
  crayon::bold("\n\nStatistical validation (LIMMA, nominal p < 0.05):\n"),
  crayon::bgGreen(crayon::black(sprintf("   Untargeted: interaction effect %d, time-main %d, group-main %d\n", sig_limma_untarg$n_sig[1], sig_limma_untarg$n_sig[2], sig_limma_untarg$n_sig[3]))),
  crayon::bgMagenta(crayon::black(sprintf("   Targeted: interaction effect %d, time-main %d, group-main %d\n\n", sig_limma_targ$n_sig[1], sig_limma_targ$n_sig[2], sig_limma_targ$n_sig[3])))
)
