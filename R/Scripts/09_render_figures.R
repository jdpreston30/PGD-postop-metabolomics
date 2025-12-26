#* 9: Render Manuscript Figures
#+ 9.1: Figure 1
fig1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
#- 1A
draw_plot(ggdraw() + draw_grob(p1A), x = 2.25, y = 7.303333, width = 4, height = 2.78) +
#- 1B
draw_plot(p1B, x = 1.18, y = 4.133333, width = 3, height = 3) +
#- 1C
draw_plot(p1C, x = 4.35, y = 4.133333, width = 3, height = 3) +
#- Labels
figure_labels(list(
A = c(2.296667, 9.98),
B = c(1.276667, 7.058333),
C = c(4.446667, 7.058333),
"Figure 1" = c(0.49, 10.43)
))
#+ 9.2: Figure 2
fig2 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
draw_plot(ggdraw() + draw_grob(grid::grobTree(ggplotGrob(p2A))), 
          x = 0.1150003, y = 5.491666, width = 7, height = 4.8) +
#- Labels
figure_labels(list(
"Figure 2" = c(0.49, 10.43)
))
#+ 9.3: Figure 3
fig3 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
#- 3A
draw_plot(p3A, x = 1.308333, y = 7.276667, width = 2.89, height = 2.89) +
#- 3B
draw_plot(p3B, x = 4.806666, y = 7.276667, width = 2.89, height = 2.89) +
#- 3C
draw_plot(p3C, x = 0.655, y = 2.263334, width = 3.6, height = 5) +
#- 3D
draw_plot(p3D, x = 4.151666, y = 2.263334, width = 3.6, height = 5) +
#- Labels
figure_labels(list(
A = c(1.341667, 9.98),
B = c(4.84, 9.98),
C = c(1.341667, 6.918334),
D = c(4.84, 6.918334),
"Figure 3" = c(0.49, 10.43)
))
#+ 9.4: Figure 4
fig4 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
#- 4A
draw_plot(p4A, x = 0.485, y = 5.349667, width = 4.5, height = 5) +
#- 4B
draw_plot(p4B, x = 5.196666, y = 7.456667-0.61, width = 2.89, height = 3.5) +
# centered
# draw_plot(p4B, x = 5.196666, y = 6.132167, width = 2.89, height = 3.5) +
#- 4C
draw_plot(ggdraw() + draw_grob(p4C), x = 1.25, y = 0.84, width = 6, height = 4.210884) +
#- Labels
figure_labels(list(
A = c(0.726667, 9.98),
B = c(5.323333, 9.98),
# centered
# B = c(5.323333, 8.928),
C = c(0.726667, 4.981667),
"Figure 4" = c(0.49, 10.43)
))
#+ 9.5: Supplementary Figure 1 
#- 9.5.1: Create separate rows for independent positioning
# PCA row (top)
pca_row <- plot_grid(
  s1D, s1E, s1F,
  nrow = 1, ncol = 3,
  align = "hv",
  axis = "tblr"
)
# PLS-DA row (bottom)
plsda_row <- plot_grid(
  s1A, s1B, s1C,
  nrow = 1, ncol = 3,
  align = "hv",
  axis = "tblr"
)
#- 9.5.2: Add column labels on top
col_labels <- plot_grid(
  ggdraw() + draw_label("12h", size = 14, fontface = "bold", fontfamily = "Arial"),
  ggdraw() + draw_label("24h", size = 14, fontface = "bold", fontfamily = "Arial"),
  ggdraw() + draw_label("12h + 24h", size = 14, fontface = "bold", fontfamily = "Arial"),
  nrow = 1, ncol = 3
)
#- 9.5.3: Manual positioning with specified shifts
sf1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(col_labels, x = 1.502167, y = 9.860, width = 7.2, height = 0.3) +
  draw_plot(pca_row, x = 0.87, y = 7.416667, width = 7.2, height = 2.5) +
  draw_label("PCA", x = 0.62, y = 8.773333, size = 14, fontface = "bold", fontfamily = "Arial", angle = 90) +
  draw_plot(plsda_row, x = 0.87, y = 5.033334, width = 7.2, height = 2.5) +
  draw_label("PLS-DA", x = 0.65, y = 6.408334, size = 14, fontface = "bold", fontfamily = "Arial", angle = 90) +
  figure_labels(list(
    "Supplementary Figure 1" = c(0.49, 10.43)
  ))
#+ 9.6: Save Figures
print_to_png(fig1, "PNG/fig1.png")
print_to_png(fig2, "PNG/fig2.png")
print_to_png(fig3, "PNG/fig3.png")
print_to_png(fig4, "PNG/fig4.png")
print_to_png(sf1, "PNG/sf1.png")
#+ 9.6: Compile as single PDF
{
  pdf("Outputs/Figures/Final/PDF/Figs1-4.pdf", width = 8.5, height = 11)
  # Page 1: Fig1
  img1 <- readPNG("Outputs/Figures/Final/PNG/fig1.png")
  grid.newpage()
  grid.raster(img1, width = unit(8.5, "inches"), height = unit(11, "inches"))
  # Page 2: Fig2
  img2 <- readPNG("Outputs/Figures/Final/PNG/fig2.png")
  grid.newpage()
  grid.raster(img2, width = unit(8.5, "inches"), height = unit(11, "inches"))
  # Page 3: Fig3
  img3 <- readPNG("Outputs/Figures/Final/PNG/fig3.png")
  grid.newpage()
  grid.raster(img3, width = unit(8.5, "inches"), height = unit(11, "inches"))
  # Page 4: Fig4
  img4 <- readPNG("Outputs/Figures/Final/PNG/fig4.png")
  grid.newpage()
  grid.raster(img4, width = unit(8.5, "inches"), height = unit(11, "inches"))
  dev.off()
}
#+ 9.7: Create individual PDFs
{
  # Figure 1
  pdf("Outputs/Figures/Final/PDF/fig1.pdf", width = 8.5, height = 11)
  img1 <- readPNG("Outputs/Figures/Final/PNG/fig1.png")
  grid.newpage()
  grid.raster(img1, width = unit(8.5, "inches"), height = unit(11, "inches"))
  dev.off()
  
  # Figure 2
  pdf("Outputs/Figures/Final/PDF/fig2.pdf", width = 8.5, height = 11)
  img2 <- readPNG("Outputs/Figures/Final/PNG/fig2.png")
  grid.newpage()
  grid.raster(img2, width = unit(8.5, "inches"), height = unit(11, "inches"))
  dev.off()
  
  # Figure 3
  pdf("Outputs/Figures/Final/PDF/fig3.pdf", width = 8.5, height = 11)
  img3 <- readPNG("Outputs/Figures/Final/PNG/fig3.png")
  grid.newpage()
  grid.raster(img3, width = unit(8.5, "inches"), height = unit(11, "inches"))
  dev.off()
  
  # Figure 4
  pdf("Outputs/Figures/Final/PDF/fig4.pdf", width = 8.5, height = 11)
  img4 <- readPNG("Outputs/Figures/Final/PNG/fig4.png")
  grid.newpage()
  grid.raster(img4, width = unit(8.5, "inches"), height = unit(11, "inches"))
  dev.off()
  
  # Supplementary Figure 1
  pdf("Outputs/Figures/Final/PDF/sf1.pdf", width = 8.5, height = 11)
  img_sf1 <- readPNG("Outputs/Figures/Final/PNG/sf1.png")
  grid.newpage()
  grid.raster(img_sf1, width = unit(8.5, "inches"), height = unit(11, "inches"))
  dev.off()
}
