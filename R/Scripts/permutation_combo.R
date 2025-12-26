 
 
   plsda_combined <- make_PCA(UFT_12and24, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PLSDA", run_permutation = TRUE)
   
   
    saveRDS(plsda_combined$scores, config$paths$permutation$plsda_combined_scores)
  saveRDS(plsda_combined$scores_df, config$paths$permutation$plsda_combined_scores_df)
  saveRDS(plsda_combined$explained, config$paths$permutation$plsda_combined_explained)
  saveRDS(plsda_combined$Y, config$paths$permutation$plsda_combined_Y)
  saveRDS(plsda_combined$cv, config$paths$permutation$plsda_combined_cv)
  saveRDS(plsda_combined$permutation, config$paths$permutation$plsda_combined_permutation)