#* 6: Subject-specific Analysis
#+ 6.1: Clean up interaction metabolites for graphing
  plot_dat <- limma_targ$interaction %>%
    arrange(p.value) %>%
    filter(p.value < 0.05) %>%
    select(Metabolite, Identified_Name, logFC) %>%
    # drop explicit duplicates you don't want
    # C16309_C18 is detected in HILIC and stronger magnitude
    filter(!Metabolite %in% c("C11384_HILIC", "C11383_HILIC", "C20893_C18", "C16309_C18")) %>%
    # start from the given name (fallback to ID if missing)
    mutate(label = coalesce(na_if(Identified_Name, ""), Metabolite)) %>%
    mutate(
      label = dplyr::case_when(
        Metabolite == "C00079_HILIC" ~ "Phenylalanine",
        Metabolite == "C20892_C18*" ~ "3-(5-oxoisoxazolin-2-yl)-alanine",
        Metabolite == "C09848_HILIC*" ~ "Citronellal",
        Metabolite == "C01767_HILIC*" ~ "Carvone",
        Metabolite == "C22300_C18" ~ "3-Hydroxyisoleucine", # keep one version
        Metabolite == "C10462_HILIC" ~ "Gingerol",
        Metabolite == "C21087_HILIC" ~ "Epoxy-4S-H4HPP",
        Metabolite == "C20818_C18" ~ "Carbapenem biosynthesis intermediate 2",
        Metabolite == "C20324_HILIC" ~ "4-OH-TMCP acetate", # short name you chose
        TRUE ~ label
      )
    ) %>%
    # append exactly one * if metabolite ID ends with * and label lacks one, then strip all * per your preference
    mutate(
      label = if_else(stringr::str_detect(Metabolite, "\\*$") & !stringr::str_detect(label, "\\*$"),
        paste0(label, "*"), label
      ),
      label = stringr::str_replace_all(label, "\\*", "")
    ) %>%
    # FINAL ordering and factor levels (do this LAST and don’t touch `label` afterwards)
    arrange((logFC)) %>% # top = most positive
    mutate(
      dir   = if_else(logFC >= 0, "Positive", "Negative"),
      label = factor(label, levels = label) # lock order
    )
#+.6.2: Graph
  #- 8.2.1: Create graph
    p_diverging <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = logFC, y = label, fill = dir)) +
      ggplot2::geom_col(width = 0.7, color = "black", linewidth = 0.6, show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = c(Positive = "#800017", Negative = "#113d6a")) +
      ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.76) +
      ggplot2::scale_x_continuous(
        limits = c(-3.1, 3.1),
        breaks = seq(-3, 3, by = 1),
        labels = seq(-3, 3, by = 1)
      ) +
      ggplot2::labs(x = expression(bold(log[2]("Fold Change"))), y = NULL) +
      ggplot2::theme_minimal(base_family = "Arial") +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.2),
        axis.text.x = ggplot2::element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = ggplot2::element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = ggplot2::element_text(size = 15, face = "bold", color = "black"),
        axis.ticks = ggplot2::element_blank()
      )

      plot_dat %>%
        select(Metabolite,label)
  #- 8.2.2: Save graph
    ggsave(
      "Outputs/LIMMA/p_diverging.png",
      p_diverging,
      width = 10,
      height = 8,
      units = "in",
      dpi = 600
    )