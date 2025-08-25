#+ Subject based analysis
#- Create NA rows for missing samples
  new_rows <- tibble::tribble(
    ~Patient, ~Time, ~Clinical_PGD,
    "H2", 12, "Y",
    "H5", 12, "Y",
    "H19", 24, "Y"
  ) %>%
    tibble::add_column(!!!setNames(rep(list(NA_real_), length(met_cols)), met_cols)) %>%
    mutate(across(c(Patient, Time, Clinical_PGD), as.factor))
#- drop any existing rows for those combos
  combined_TFT_LMM_i <- combined_TFT %>%
    select(-c(Sex, Age, Sample_ID)) %>%
    filter(
      !(Patient == "H2" & Time == 12 & Clinical_PGD == "Y"),
      !(Patient == "H5" & Time == 12 & Clinical_PGD == "Y"),
      !(Patient == "H19" & Time == 24 & Clinical_PGD == "Y")
    )
#- bind them back
  combined_TFT_LMM <- bind_rows(combined_TFT_LMM_i, new_rows) %>%
    arrange(Patient, Time)
    
# --- FIX FOR limma ORIENTATION + FULL WORKFLOW ---

# Ensure factors have the right baseline levels
combined_TFT_LMM$Time <- factor(combined_TFT_LMM$Time, levels = c("12", "24"))
combined_TFT_LMM$Clinical_PGD <- factor(combined_TFT_LMM$Clinical_PGD, levels = c("N", "Y"))

# Identify metabolite columns (4+)
met_cols <- names(combined_TFT_LMM)[4:ncol(combined_TFT_LMM)]

# Build expression matrix as features x samples (rows = metabolites, cols = samples)
# limma expects columns to be samples (arrays)
expr <- t(as.matrix(combined_TFT_LMM[, met_cols, drop = FALSE]))
colnames(expr) <- paste0(combined_TFT_LMM$Patient, "_", combined_TFT_LMM$Time)

# Design matrix (rows must match number of columns in expr == nrow(combined_TFT_LMM))
design <- model.matrix(~ Time * Clinical_PGD, data = combined_TFT_LMM)

# Blocking factor for repeated measures (length must equal ncol(expr))
block <- combined_TFT_LMM$Patient

# Estimate within-subject correlation and fit
corfit <- limma::duplicateCorrelation(expr, design, block = block)
fit <- limma::lmFit(expr, design, block = block, correlation = corfit$consensus)
fit <- limma::eBayes(fit)

# Inspect coefficient names if needed:
# colnames(design)

# Extract Group x Time interaction results (flat vs non-flat between groups)
tt_inter <- limma::topTable(fit, coef = "Time24:Clinical_PGDY", number = Inf, sort.by = "P")

# Tidy tibble with metabolite, estimate, p, FDR
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


# Compute "fold change of fold changes" = log2FC_Y − log2FC_N (i.e., log2((24/12)_Y / (24/12)_N))
# Plot all metabolites as points ranked by |Δ log2FC|.
# Color: p<0.05 & Δ>0 = red (“PGD up”); p<0.05 & Δ<0 = blue (“PGD down”); others = light gray.
# Use Identified_Name for labels/ordering.

library(dplyr)
library(tidyr)
library(ggplot2)

# --- 1) Identify metabolite columns (everything from col 4 on) ---
met_cols <- names(combined_TFT_LMM)[4:ncol(combined_TFT_LMM)]

# --- 2) Long -> within-subject log2FC (24 − 12) per patient, per metabolite ---
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

# --- 3) Join names (for plotting labels) ---
fc_df <- fc_df %>%
  left_join(combined_key %>% select(Metabolite, Identified_Name), by = "Metabolite") %>%
  mutate(Identified_Name = if_else(is.na(Identified_Name), Metabolite, Identified_Name))

# --- 4) For each metabolite, compute group means (ΔY, ΔN) and Δ = ΔY − ΔN, plus p-value ---
summary_df <- fc_df %>%
  group_by(Metabolite, Identified_Name, Clinical_PGD) %>%
  summarise(
    mean_fc = mean(log2FC, na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  tidyr::pivot_wider(names_from = Clinical_PGD, values_from = mean_fc, names_prefix = "mean_fc_") %>%
  # ensure both columns exist; if any missing, set to NA
  mutate(
    mean_fc_N = ifelse(is.na(mean_fc_N), NA_real_, mean_fc_N),
    mean_fc_Y = ifelse(is.na(mean_fc_Y), NA_real_, mean_fc_Y)
  ) %>%
  # Δ log2FC (log2 FC-of-FC)
  mutate(delta_log2FC = mean_fc_Y - mean_fc_N) %>%
  ungroup()

# compute a p-value per metabolite by testing log2FC between groups (Y vs N)
pvals <- fc_df %>%
  group_by(Metabolite) %>%
  do({
    x <- .
    # need both groups present with >=2 total obs ideally
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

summary_df <- summary_df %>%
  left_join(pvals, by = "Metabolite")

# --- 5) Define color categories based on sign and p<0.05 ---
summary_df <- summary_df %>%
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

# --- 6) Rank by |Δ log2FC| (largest magnitude first) ---
plot_df <- summary_df %>%
  # Primary sort: delta_log2FC ascending (most negative → most positive)
  # Tie-breakers: smaller p first, then Identified_Name (stable/deterministic)
  arrange(delta_log2FC, P, Identified_Name, Metabolite) %>%
  mutate(Rank = row_number())
# --- PLOT (same coloring as before) ---
ggplot(plot_df, aes(x = Rank, y = delta_log2FC)) +
  geom_point(
    data = subset(plot_df, ColorCat == "Other"),
    color = "grey80", size = 2, alpha = 0.9
  ) +
  geom_point(
    data = subset(plot_df, ColorCat == "PGD up (p<0.05)"),
    color = "#d62728", size = 2.4, alpha = 1
  ) +
  geom_point(
    data = subset(plot_df, ColorCat == "PGD down (p<0.05)"),
    color = "#1f77b4", size = 2.4, alpha = 1
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  labs(
    x = "Rank (most No‑PGD ←  Δlog2FC  → most PGD)",
    y = "Δ log2FC  = (mean log2FC_Y − mean log2FC_N)",
    title = "Fold‑change of fold‑changes ranked from No‑PGD to PGD"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())


# --- 6) Rank by |Δ log2FC| (largest magnitude first) ---
plot_df_clay <-  summary_df %>%
  # Primary sort: delta_log2FC ascending (most negative → most positive)
  # Tie-breakers: smaller p first, then Identified_Name (stable/deterministic)
filter(P < 0.05)
