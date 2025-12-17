#+ Filter down
#- Remove metadata
UFT_features <- UFT |>
  select(-c(Patient, Sample, Sample_ID, Time, severe_PGD, PGD_grade_tier, any_PGD))
#- Make a tibble of features
feature_names <- tibble(
  feature = colnames(UFT_features)
) |>
  separate(feature, into = c("column", "mz", "rt"), sep = "_", remove = FALSE) |>
  mutate(
    mz = as.numeric(mz),
    rt = as.numeric(rt)
  ) |>
  arrange(mz)


  feature_names |>
    filter(column == "C18")
#- Set vectors
malondialdehyde <- c(157.1223)

#+ Find matching features
matched_features <- find_features_by_mz(feature_names, malondialdehyde, ppm_tolerance = 10)
matched_features

