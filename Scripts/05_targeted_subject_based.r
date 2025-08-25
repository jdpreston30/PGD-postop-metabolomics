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
  #+ 1.3: Compute and plot fold change of fold changes
    #- 1.3.1: Compute within-subject log2FC (24 − 12) per patient, per metabolite
      fc_df <- combined_TFT_LMM %>%
        select(Patient, Time, Clinical_PGD, all_of(met_cols)) %>%
        pivot_longer(cols = all_of(met_cols), names_to = "Metabolite", values_to = "Intensity") %>%
        mutate(
          Time = factor(Time, levels = c("12", "24")),
          Clinical_PGD = factor(Clinical_PGD, levels = c("N", "Y"))
        ) %>%
        pivot_wider(names_from = Time, values_from = Intensity) %>%
        filter(!is.na(`12`), !is.na(`24`)) %>%
        mutate(log2FC = `24` - `12`) %>%
        select(Patient, Clinical_PGD, Metabolite, log2FC)
    #- 1.3.2: Join names (for plotting labels)
      fc_df <- fc_df %>%
        left_join(combined_key %>% select(Metabolite, Identified_Name), by = "Metabolite") %>%
        mutate(Identified_Name = if_else(is.na(Identified_Name), Metabolite, Identified_Name))
    #- 1.3.3: For each metabolite, compute group means (ΔY, ΔN) and Δ = ΔY − ΔN, plus p-value
      summary_dfi <- fc_df %>%
        group_by(Metabolite, Identified_Name, Clinical_PGD) %>%
        summarise(
          mean_fc = mean(log2FC, na.rm = TRUE),
          .groups = "drop_last"
        ) %>%
        tidyr::pivot_wider(names_from = Clinical_PGD, values_from = mean_fc, names_prefix = "mean_fc_") %>%
        mutate(
          mean_fc_N = ifelse(is.na(mean_fc_N), NA_real_, mean_fc_N),
          mean_fc_Y = ifelse(is.na(mean_fc_Y), NA_real_, mean_fc_Y)
        ) %>%
        mutate(delta_log2FC = mean_fc_Y - mean_fc_N) %>%
        ungroup()
    #- 1.3.4: Compute a p-value per metabolite by testing log2FC between groups (Y vs N)
      pvals <- fc_df %>%
        group_by(Metabolite) %>%
        do({
          x <- .
          y_vals <- x$log2FC[x$Clinical_PGD == "Y"]
          n_vals <- x$log2FC[x$Clinical_PGD == "N"]
          if (length(y_vals) >= 2 && length(n_vals) >= 2) {
            pv <- tryCatch(t.test(y_vals, n_vals, var.equal = FALSE)$p.value, error = function(e) NA_real_)
          } else {
            pv <- NA_real_
          }
          tibble(P = pv)
        }) %>%
        ungroup()
    #- 1.3.5: Define color categories based on sign and p<0.05
      summary_df <- summary_dfi %>%
        left_join(pvals, by = "Metabolite") %>%
        mutate(
          Direction = case_when(
            delta_log2FC > 0 ~ "PGD up",
            delta_log2FC < 0 ~ "PGD down",
            TRUE ~ "flat"
          ),
          Sig = if_else(!is.na(P) & P < 0.05, "sig", "ns"),
          ColorCat = case_when(
            Sig == "sig" & Direction == "PGD up" ~ "PGD up (p<0.05)",
            Sig == "sig" & Direction == "PGD down" ~ "PGD down (p<0.05)",
            TRUE ~ "Other"
          )
        )      
    #- 1.3.6: Rank by |Δ log2FC| (largest magnitude first)
      plot_df <- summary_df %>%
        arrange(delta_log2FC, P, Identified_Name, Metabolite) %>%
        mutate(Rank = row_number())
    # --- build the two name lists with line breaks ---
      red_hits <- plot_df %>%
        filter(ColorCat == "PGD up (p<0.05)") %>%
        arrange(Rank) %>%
        distinct(Identified_Name) %>%
        pull(Identified_Name)
      blue_hits <- plot_df %>%
        filter(ColorCat == "PGD down (p<0.05)") %>%
        arrange(Rank) %>%
        distinct(Identified_Name) %>%
        pull(Identified_Name)
      red_str <- if (length(red_hits)) paste(red_hits, collapse = "\n") else "none"
      blue_str <- if (length(blue_hits)) paste(blue_hits, collapse = "\n") else "none"
    # --- base plot ---
      p_fcfc <- ggplot(plot_df, aes(x = Rank, y = delta_log2FC)) +
        geom_point(
          data = subset(plot_df, ColorCat == "Other"),
          color = "grey80", size = 2, alpha = 0.9
        ) +
        geom_point(
          data = subset(plot_df, ColorCat == "PGD up (p<0.05)"),
          color = "#94001E", size = 2.4, alpha = 1
        ) +
        geom_point(
          data = subset(plot_df, ColorCat == "PGD down (p<0.05)"),
          color = "#03507D", size = 2.4, alpha = 1
        ) +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
        scale_x_continuous(expand = c(0.01, 0.01)) +
        labs(
          x = "Rank (most No-PGD ←  Δlog2FC  → most PGD)",
          y = "Δ log2FC = (mean log2FC_Y − mean log2FC_N)",
          title = "Fold-change of fold-changes ranked from No-PGD to PGD"
        ) +
        theme_bw(base_size = 12) +
        theme(panel.grid.minor = element_blank()) +
        annotate("text",
          x = max(plot_df$Rank) * 0.8, # right-ish
          y = max(plot_df$delta_log2FC) * 0.9, # near the top
          label = red_str,
          hjust = 1, vjust = 1,
          color = "#94001E", size = 3.2, fontface = "bold"
        ) +
          annotate("text",
            x = max(plot_df$Rank) * 0.2, # left-ish
            y = min(plot_df$delta_log2FC) * 0.9, # near the bottom
            label = blue_str,
            hjust = 0, vjust = 0,
            color = "#03507D", size = 3.2, fontface = "bold"
          )
# GOAL
# Make TWO bar plots centered at 0 (bars go up if +, down if −), with exposed-side error bars.
# Plot A = metabolites where PGD Y is + (paired with N = + or −)  → "++ / +-"
# Plot B = metabolites where PGD Y is − (paired with N = − or +)  → "-- / -+"
# Colors: bars Y = "#94001E" (deep red), N = "#03507D" (deep blue).
# Also add in-panel name lists (red list top-right; blue list bottom-left), using the same colors.

library(dplyr)
library(tidyr)
library(ggplot2)

# --- 0) Build within-subject log2FC for each patient/metabolite ---
# fc_df must exist already; if not, (re)build it:
# fc_df: Patient, Clinical_PGD, Metabolite, Identified_Name, log2FC (24−12)

# --- 1) Mean ± SE per group & metabolite ---
sum_df <- fc_df %>%
  group_by(Clinical_PGD, Metabolite, Identified_Name) %>%
  summarise(
    mean_fc = mean(log2FC, na.rm = TRUE),
    se      = sd(log2FC, na.rm = TRUE) / sqrt(sum(!is.na(log2FC))),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = ifelse(mean_fc < 0, mean_fc - se, mean_fc),
    ymax = ifelse(mean_fc < 0, mean_fc, mean_fc + se)
  )

# --- 2) Classify each metabolite by sign of mean FC in Y and N ---
sum_wide <- sum_df %>%
  pivot_wider(
    id_cols = c(Metabolite, Identified_Name),
    names_from = Clinical_PGD,
    values_from = c(mean_fc, se),
    names_sep = "_"
  ) %>%
  mutate(
    dir_Y = case_when(mean_fc_Y > 0 ~ "up", mean_fc_Y < 0 ~ "down", TRUE ~ "flat"),
    dir_N = case_when(mean_fc_N > 0 ~ "up", mean_fc_N < 0 ~ "down", TRUE ~ "flat"),
    Category = case_when(
      dir_Y == "up" & dir_N == "up" ~ "++",
      dir_Y == "up" & dir_N == "down" ~ "+-",
      dir_Y == "down" & dir_N == "down" ~ "--",
      dir_Y == "down" & dir_N == "up" ~ "-+",
      TRUE ~ "other"
    )
  ) %>%
  select(Metabolite, Identified_Name, Category)

sum_long <- sum_df %>% left_join(sum_wide, by = c("Metabolite", "Identified_Name"))

# --- 3) Order metabolites within each plot by max |mean_fc| across groups ---
sum_long <- sum_long %>%
  group_by(Category, Metabolite, Identified_Name) %>%
  mutate(ord = max(abs(mean_fc), na.rm = TRUE)) %>%
  ungroup()

# --- 4) Helper plotting function (bars with exposed-side error bars, centered at 0) ---
plot_bars <- function(dat, title_txt, red_names = character(0), blue_names = character(0)) {
  # order x by effect size
  dat <- dat %>%
    arrange(desc(ord)) %>%
    mutate(Identified_Name = factor(Identified_Name, levels = unique(Identified_Name)))

  p <- ggplot(
    dat,
    aes(x = Identified_Name, y = mean_fc, fill = Clinical_PGD)
  ) +
    geom_col(
      position = position_dodge(width = 0.72),
      width = 0.64, linewidth = 0.25, color = "grey20"
    ) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      position = position_dodge(width = 0.72),
      width = 0.18, linewidth = 0.4
    ) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
    scale_fill_manual(
      values = c("N" = "#03507D", "Y" = "#94001E"),
      labels = c("N" = "No PGD", "Y" = "PGD"),
      name = NULL
    ) +
    labs(
      x = NULL,
      y = "Within-subject log2FC (24 − 12)",
      title = title_txt
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "top"
    )

  # Add in-plot lists (each name on its own line)
  if (length(red_names) == 0) red_names <- "none"
  if (length(blue_names) == 0) blue_names <- "none"
  red_str <- paste(red_names, collapse = "\n")
  blue_str <- paste(blue_names, collapse = "\n")

  # compute panel extents from data actually shown
  xr <- ggplot_build(p)$layout$panel_params[[1]]$x.range
  yr <- range(dat$mean_fc, dat$ymin, dat$ymax, na.rm = TRUE)

  p +
    annotate(
      "text",
      x = mean(xr) + 0.3 * diff(xr), # right-ish middle
      y = yr[2] - 0.1 * diff(yr), # near top
      label = red_str, hjust = 1, vjust = 1,
      color = "#94001E", size = 3.2, fontface = "bold"
    ) +
    annotate(
      "text",
      x = mean(xr) - 0.3 * diff(xr), # left-ish middle
      y = yr[1] + 0.1 * diff(yr), # near bottom
      label = blue_str, hjust = 0, vjust = 0,
      color = "#03507D", size = 3.2, fontface = "bold"
    )
}

# --- 5) Build the two datasets ---
dat_pos <- sum_long %>% filter(Category %in% c("++", "+-"))
dat_neg <- sum_long %>% filter(Category %in% c("--", "-+"))
