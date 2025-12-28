 
 
   plsda_combined <- make_PCA(UFT_12and24, comp_x = 1, comp_y = 2, group_var = "severe_PGD", method = "PLSDA", run_permutation = TRUE)

# Remove plot object and save as list
plsda_combined$plot <- NULL
saveRDS(plsda_combined, config$paths$permutation$plsda_combined)