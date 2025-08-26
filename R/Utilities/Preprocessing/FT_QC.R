#' QC-aware half-min + log2 processing for metabolomics feature tables
#'
#' @title FT_QC: Prepare wide metabolomics tables for analysis (MetaboAnalyst-ready and filtered)
#'
#' @description
#' Processes a wide feature table (rows = samples, columns = features) by:
#' (1) computing per-feature detection rate as the proportion of **zeros among non-NA** values;
#' (2) keeping features with detection rate ≥ (1 - `zero_thresh`);
#' (3) imputing **zeros only** via half of the smallest positive value (computed across all rows, including QC);
#' (4) log-transforming features (default `log2`);
#' (5) optionally removing QC rows (e.g., NIST, q1–q9) from outputs.
#'
#' **Semantics:** `0` = measured but not detected (LOD) → counted as "missing" for detection
#' and imputed with half-min; `NA` = not present / not measured (e.g., missing block) → **excluded**
#' from the detection denominator and **not** imputed.
#'
#' @param df A data frame / tibble with one row per sample. Must contain an ID column and any metadata
#'   listed in `common_cols`. All other columns are treated as numeric feature columns.
#' @param id_col Character scalar. Name of the sample identifier column (default: `"Sample_ID"`).
#' @param common_cols Character vector of columns to preserve as metadata (default includes `id_col`).
#'   All columns not in `common_cols` are assumed numeric features.
#' @param qc_remove Logical. If `TRUE`, QC rows (matching `qc_pattern`) are dropped **after**
#'   half-min + log transform and after feature filtering. QC rows are still used to compute
#'   half-min and detection rates. Default `TRUE`.
#' @param qc_pattern Character regex used to identify QC samples by `id_col`
#'   (default: `"^(NIST1|NIST2|q[1-9]+)$"`).
#' @param zero_thresh Numeric in [0,1]. A feature is kept if
#'   `(zeros among non-NA) / (non-NA count) <= zero_thresh`. Default `0.20` (keep features with ≥80% detection).
#' @param min_nonNA Integer. Minimum number of non-NA observations required per feature to be considered.
#'   Features with `non-NA < min_nonNA` are treated as failing detection (i.e., dropped). Default `3`.
#' @param log_fun Function to apply for the log transform (default `log2`). Must accept numeric input and return numeric.
#'
#' @details
#' - **Zeros vs NAs**:
#'   - Zeros are interpreted as *undetected* and contribute to the detection-rate numerator; they are
#'     imputed using half of the smallest positive value per feature (half-min).
#'   - NAs indicate *not measured / not present* (e.g., missing mode) and are excluded from the detection denominator;
#'     they remain NA (no half-min imputation).
#' - **Half-min calculation**: The smallest positive value is computed per feature across **all rows, including QC**.
#'   If a feature has no positive values at all, its half-min is set to `NA_real_` and zeros remain NA after imputation.
#' - **QC rows**: When `qc_remove = TRUE`, QC rows are excluded from the returned tibbles but still influence
#'   detection-rate and half-min estimation.
#'
#' @return A list with two tibbles:
#' \describe{
#'   \item{`full`}{All metadata + **all** feature columns after half-min and log transform (no feature dropping).
#'                 QC rows removed if `qc_remove = TRUE`.}
#'   \item{`filtered`}{Same as `full`, but restricted to features passing the detection filter
#'                     (≤ `zero_thresh` zeros among non-NA and `>= min_nonNA` non-NA observations).}
#' }
#'
#' @examples
#' \dontrun{
#' meta_cols <- c("Sample_ID", "Patient", "Time", "Sex", "Age", "Clinical_PGD")
#' res <- FT_QC(
#'   df          = C18_HILIC_targeted,
#'   id_col      = "Sample_ID",
#'   common_cols = meta_cols,
#'   qc_remove   = TRUE
#' )
#' metabo_full <- res$full # MetaboAnalyst-ready, keeps all features, no QC rows
#' metabo_filtered <- res$filtered # Same but drops features failing the 20% rule
#' }
#'
#' @seealso
#'   \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{across}}, \code{\link[tidyr]{pivot_longer}}
#'
#' @importFrom dplyr mutate across all_of summarise select filter if_else
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_detect
#'
FT_QC <- function(
    df,
    id_col = "Sample_ID",
    common_cols = c("Sample_ID"),
    qc_remove = TRUE,
    qc_pattern = "^(NIST1|NIST2|q[1-9]+)$",
    zero_thresh = 0.20, # keep if (zeros / non-NA) <= 20%
    min_nonNA = 3, # require at least this many non-NA values per feature
    log_fun = log2) {
  stopifnot(id_col %in% names(df))
  common_cols <- unique(c(common_cols, id_col))
  feature_cols <- setdiff(names(df), common_cols)

  # Ensure feature cols are numeric
  df <- df %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(feature_cols), ~ suppressWarnings(as.numeric(.x))))

  # 1) % zeros among non-NA (dynamic denominator); enforce min_nonNA
  zero_pct <- df %>%
    dplyr::summarise(dplyr::across(
      dplyr::all_of(feature_cols),
      ~ {
        n_nonNA <- sum(!is.na(.x))
        if (n_nonNA < min_nonNA) {
          1 # treat as 100% "missing" so it fails the filter
        } else {
          sum(.x == 0, na.rm = TRUE) / n_nonNA
        }
      },
      .names = "{.col}"
    )) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "feature", values_to = "zero_pct")

  keep_features <- zero_pct %>%
    dplyr::filter(zero_pct <= zero_thresh) %>%
    dplyr::pull(feature)

  # 2) Half-min replacement ONLY for zeros (NAs remain NA)
  halfmin_vec <- vapply(feature_cols, function(col) {
    v <- df[[col]]
    mp <- suppressWarnings(min(v[v > 0], na.rm = TRUE))
    if (is.infinite(mp)) NA_real_ else 0.5 * mp
  }, numeric(1))
  names(halfmin_vec) <- feature_cols

  df_halfmin <- df %>%
    dplyr::mutate(dplyr::across(
      dplyr::all_of(feature_cols),
      ~ dplyr::if_else(.x == 0, halfmin_vec[cur_column()], .x)
    ))

  # 3) Log transform (NAs stay NA)
  df_log <- df_halfmin %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(feature_cols), ~ log_fun(.x)))

  # 4) Outputs
  full_out <- df_log
  filtered_out <- df_log %>%
    dplyr::select(
      dplyr::all_of(id_col), # always first
      dplyr::all_of(setdiff(common_cols, id_col)), # other metadata
      dplyr::all_of(keep_features) # features
    )

  # 5) Optionally remove QC rows at the end
  if (isTRUE(qc_remove)) {
    full_out <- full_out %>% filter(!stringr::str_detect(.data[[id_col]], qc_pattern))
    filtered_out <- filtered_out %>% filter(!stringr::str_detect(.data[[id_col]], qc_pattern))
  }

  list(full = full_out, filtered = filtered_out)
}
