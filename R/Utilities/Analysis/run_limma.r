#' Run limma (duplicateCorrelation) for repeated-measures metabolomics
#' and emit contrast tables, with optional MetaboAnalyst/Mummichog p-lists
#' and optional targeted-name annotation.
#'
#' @param data            data.frame/tibble, rows = samples, numeric feature cols + metadata.
#' @param time_var        character. Time column (e.g., "Time").
#' @param group_var       character. Group column (e.g., "Clinical_PGD").
#' @param block_var       character. Subject/Patient column (e.g., "Patient").
#' @param time_levels     character(2). Baseline first, e.g. c("12","24").
#' @param group_levels    character(2). Baseline first, e.g. c("N","Y").
#' @param meta_cols       character vector of metadata columns to EXCLUDE from features.
#'                        If NULL, defaults to c(time_var, group_var, block_var).
#' @param export_prefix   optional character. If non-NULL, CSVs are written using this prefix.
#' @param export_dir      character. Folder for outputs when exporting. Default "Outputs".
#' @param metaboanalyst   logical. If TRUE, produce MetaboAnalyst/Mummichog-ready tables (m.z, p.value, mode, r.t)
#'                        and export them when export_prefix is set. Default TRUE.
#' @param targeted_rename logical. If TRUE, left-join a rename_key onto all stats outputs. Default FALSE.
#' @param rename_key      tibble with columns c("Metabolite", <any other cols to append>), e.g.
#'                        ("Metabolite","Identified_Name"). Used only if targeted_rename = TRUE.
#'
#' @return list with elements:
#'   $interaction, $time_main, $group_main (stats tables: Metabolite, t.score, p.value, FDR, logFC, Rank [+ any joined cols])
#'   If metaboanalyst=TRUE: $mummichog (list of three tibbles), else NULL.
#'   $coef_names (named character vector)
run_limma_three_contrasts <- function(
    data,
    time_var = "Time",
    group_var = "Clinical_PGD",
    block_var = "Patient",
    time_levels = c("12", "24"),
    group_levels = c("N", "Y"),
    meta_cols = NULL,
    export_prefix = NULL,
    export_dir = "Outputs",
    metaboanalyst = TRUE,
    targeted_rename = FALSE,
    rename_key = NULL) {
  stopifnot(all(c(time_var, group_var, block_var) %in% names(data)))

  df <- data
  df[[time_var]] <- factor(df[[time_var]], levels = time_levels)
  df[[group_var]] <- factor(df[[group_var]], levels = group_levels)
  df[[block_var]] <- factor(df[[block_var]])

  if (is.null(meta_cols)) meta_cols <- c(time_var, group_var, block_var)
  feat_cols <- setdiff(names(df), meta_cols)

  # coerce feature cols to numeric
  df[feat_cols] <- lapply(df[feat_cols], function(x) suppressWarnings(as.numeric(x)))

  # expression (features x samples) and weights for NAs
  E <- t(as.matrix(df[, feat_cols, drop = FALSE]))
  W <- ifelse(is.na(E), 0, 1) # 1 = observed, 0 = missing
  E0 <- E
  E0[is.na(E0)] <- 0 # values at weight=0 are ignored

  # design: ~ time * group
  fml <- stats::as.formula(paste("~", time_var, "*", group_var))
  design <- model.matrix(fml, data = df)
  block <- df[[block_var]]

  # expected coefficient names
  t2 <- levels(df[[time_var]])[2]
  g2 <- levels(df[[group_var]])[2]
  coef_time <- paste0(time_var, t2) # e.g. "Time24"
  coef_group <- paste0(group_var, g2) # e.g. "Clinical_PGDY"
  coef_interact <- paste0(coef_time, ":", coef_group) # e.g. "Time24:Clinical_PGDY"

  miss <- setdiff(c(coef_time, coef_group, coef_interact), colnames(design))
  if (length(miss)) stop("Expected coefficients not found in design: ", paste(miss, collapse = ", "))

  # limma with duplicateCorrelation + weights
  corfit <- limma::duplicateCorrelation(E0, design, block = block, weights = W)
  fit <- limma::lmFit(E0, design, block = block, correlation = corfit$consensus, weights = W)
  fit <- limma::eBayes(fit)

  # tidy topTable
  tidy_tt <- function(fit, coef_name) {
    tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
    out <- tibble::tibble(
      Metabolite = rownames(tt),
      t.score    = unname(tt[, "t"]),
      p.value    = unname(tt[, "P.Value"]), # raw p: what MetaboAnalyst/Mummichog expects
      FDR        = unname(tt[, "adj.P.Val"]),
      logFC      = unname(tt[, "logFC"])
    ) |>
      dplyr::mutate(Rank = dplyr::row_number()) |>
      dplyr::select(Metabolite, t.score, p.value, FDR, logFC, Rank)
    out
  }

  # optional targeted rename joiner
  maybe_annotate <- function(tbl) {
    if (isTRUE(targeted_rename)) {
      if (is.null(rename_key) || !"Metabolite" %in% names(rename_key)) {
        stop("targeted_rename=TRUE but rename_key is NULL or missing 'Metabolite' column.")
      }
      # Join ALL extra columns from rename_key (beyond 'Metabolite')
      extra_cols <- setdiff(names(rename_key), "Metabolite")
      if (length(extra_cols) == 0) {
        return(tbl)
      }
      tbl <- dplyr::left_join(tbl, rename_key, by = "Metabolite")
      # move the extra columns right after Metabolite
      tbl <- dplyr::relocate(tbl, dplyr::any_of(extra_cols), .after = Metabolite)
    }
    tbl
  }

  # MetaboAnalyst/Mummichog converter
  to_mummichog <- function(stats_tbl) {
    stats_tbl |>
      dplyr::transmute(
        Feature = Metabolite,
        p.value = p.value,
        mode = dplyr::case_when(
          stringr::str_starts(Feature, "HILIC") ~ "positive",
          stringr::str_starts(Feature, "C18") ~ "negative",
          TRUE ~ NA_character_
        ),
        m.z = suppressWarnings(as.numeric(stringr::str_extract(Feature, "(?<=_)[0-9.]+"))),
        r.t = suppressWarnings(as.numeric(
          stringr::str_extract(Feature, "_[0-9.]+_([0-9.]+)") |>
            stringr::str_extract("[0-9.]+$")
        ))
      ) |>
      dplyr::filter(!is.na(mode) & !is.na(m.z) & !is.na(r.t)) |>
      dplyr::select(m.z, p.value, mode, r.t)
  }

  # build stats tables
  interaction_tbl <- tidy_tt(fit, coef_interact) |> maybe_annotate()
  time_tbl <- tidy_tt(fit, coef_time) |> maybe_annotate()
  group_tbl <- tidy_tt(fit, coef_group) |> maybe_annotate()

  # metaboanalyst outputs (optional)
  mummi <- NULL
  if (isTRUE(metaboanalyst)) {
    mummi <- list(
      interaction = to_mummichog(interaction_tbl),
      time_main   = to_mummichog(time_tbl),
      group_main  = to_mummichog(group_tbl)
    )
  }

  # optional export
  if (!is.null(export_prefix)) {
    dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

    # stats tables
    readr::write_csv(interaction_tbl, file.path(export_dir, paste0(export_prefix, "_interaction.csv")))
    readr::write_csv(time_tbl, file.path(export_dir, paste0(export_prefix, "_time.csv")))
    readr::write_csv(group_tbl, file.path(export_dir, paste0(export_prefix, "_group.csv")))

    # MetaboAnalyst/Mummichog lists
    if (isTRUE(metaboanalyst)) {
      readr::write_csv(mummi$interaction, file.path(export_dir, paste0(export_prefix, "_mummichog_interaction.csv")))
      readr::write_csv(mummi$time_main, file.path(export_dir, paste0(export_prefix, "_mummichog_time.csv")))
      readr::write_csv(mummi$group_main, file.path(export_dir, paste0(export_prefix, "_mummichog_group.csv")))
    }
  }

  list(
    interaction = interaction_tbl,
    time_main   = time_tbl,
    group_main  = group_tbl,
    mummichog   = mummi,
    coef_names  = c(time = coef_time, group = coef_group, interaction = coef_interact)
  )
}
