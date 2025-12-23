#* 8: Assign Figures
#+ 8.1: Figure 1
p1A <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/fig1A.png")))
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

#+ 8.5: Supplemental Figure 1
s1A <- plsda_12h$plot
s1B <- plsda_24h$plot
s1C <- plsda_combined$plot
s1D <- pca_12h$plot
s1E <- pca_24h$plot
s1F <- pca_combined$plot
