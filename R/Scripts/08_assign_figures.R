#* 8: Assign Figures
#+ 8.1: Figure 1
p1A <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/fig1A.png")))
# PLSDA plots - no layer removal needed (no p-value annotations in these plots)
p1B <- plsda_12h$plot
p1C <- plsda_24h$plot
#+ 8.2: Figure 2
p2A <- MFN_plot
#+ 8.3: Figure 3
p3A <- volc_12$volcano_plot
p3B <- volc_24$volcano_plot
p3C <- div_bars_12
p3D <- div_bars_24
#+ 8.3: Figure 4
p4A <- div_bars_intx
p4B <- succinate_summary_plot
p4C <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/fig4C.png")))
#+ 8.5: Supplemental Figure 1 (and shorten labels)
#- 8.5.1: Relabel PGD (no need to remove layers - these plots don't have p-value annotations)
s1A <- relabel_pgd(plsda_12h$plot)
s1B <- relabel_pgd(plsda_24h$plot)
s1C <- relabel_pgd(plsda_combined$plot)
s1D <- relabel_pgd(pca_12h$plot)
s1E <- relabel_pgd(pca_24h$plot)
s1F <- relabel_pgd(pca_combined$plot)
#+ 8.6: Pull combined PLSDA for graphical abstract
ggsave("Outputs/Figures/Raw/graphical_abstract_figure.png", s1C, width = 3, height = 3, dpi = 800, bg = "transparent")
