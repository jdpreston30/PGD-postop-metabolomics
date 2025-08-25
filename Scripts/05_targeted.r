#* 5: Targeted  
#+ 5.1: Merge TFT with path data
  anova_results <- combined_TFT %>%
    # make sure factors are factors
    dplyr::mutate(
      Time = as.factor(Time),
      Clinical_PGD = as.factor(Clinical_PGD)
    ) %>%
    tidyr::pivot_longer(
      cols = -all_of(meta_cols),
      names_to = "Metabolite",
      values_to = "Value"
    ) %>%
    dplyr::group_by(Metabolite) %>%
    dplyr::summarise(
      {
        fit <- aov(Value ~ Time * Clinical_PGD, data = dplyr::cur_data())
        tbl <- broom::tidy(fit)
        tibble::tibble(
          p_time  = tbl$p.value[tbl$term == "Time"],
          p_pgd   = tbl$p.value[tbl$term == "Clinical_PGD"],
          p_inter = tbl$p.value[tbl$term == "Time:Clinical_PGD"]
        )
      },
      .groups = "drop"
    ) %>%
    filter(p_inter < 0.05) %>%
    # filter(p_time < 0.05 & p_pgd < 0.05) %>% # filter to those sig for time AND PGD
    left_join(combined_key, by = "Metabolite") %>%
    arrange(p_pgd)
#+ 5.2: Filter down to those chosen metabolites above; plot
  #- 
  combined_TFT_plot <- combined_TFT %>%
    select(all_of(meta_cols), any_of(anova_results$Metabolite))

# 1) Long-format + keep only the metabolites you want to plot
long_df <- combined_TFT_plot %>%
  pivot_longer(
    cols = -all_of(meta_cols),
    names_to = "Metabolite",
    values_to = "Value"
  ) %>%
  mutate(
    Time = factor(Time, levels = c("12", "24")), # or c(12,24) if numeric
    Clinical_PGD = factor(Clinical_PGD, levels = c("N", "Y"))
  )

# 2) Attach identified names from ANOVA results (assumes `anova_results` has Metabolite + Identified_Name)
long_df <- long_df %>%
  left_join(
    anova_results %>% select(Metabolite, Identified_Name),
    by = "Metabolite"
  ) %>%
  mutate(
    # facet label: "Identified_Name (Metabolite)" so you can tell them apart if names repeat
    Facet = if_else(!is.na(Identified_Name),
      paste0(Identified_Name, " (", Metabolite, ")"),
      Metabolite
    )
  )

# 3) Compute means by Time x PGD x Metabolite
means_df <- long_df %>%
  group_by(Metabolite, Facet, Time, Clinical_PGD) %>%
  summarise(mean_val = mean(Value, na.rm = TRUE), .groups = "drop")
library(ggforce)

# How many pages you’ll need (depends on your nrow/ncol layout)
n_pages <- n_pages(p_lines + facet_wrap_paginate(~Facet, ncol = 3, nrow = 4))

# Loop over pages
for (i in seq_len(n_pages)) {
  g <- ggplot(means_df, aes(x = Time, y = mean_val, group = Clinical_PGD, colour = Clinical_PGD)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("N" = "#1f77b4", "Y" = "#d62728"), name = "PGD") +
    labs(x = NULL, y = "Mean (original units)", title = NULL) +
    ggforce::facet_wrap_paginate(~Facet, scales = "free_y", ncol = 3, nrow = 4, page = i) +
    theme_minimal(base_family = "Arial") +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )

  ggsave(paste0("p_lines_page_", i, ".png"), g, width = 10, height = 8, dpi = 600, bg = "white")
}


p_lines <- ggplot(means_df, aes(x = Time, y = mean_val, group = Clinical_PGD, colour = Clinical_PGD)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("N" = "#1f77b4", "Y" = "#d62728"), name = "PGD") +
  labs(x = NULL, y = "Mean (original units)", title = NULL) +
  facet_wrap(~Facet, scales = "free_y") +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

p_lines
