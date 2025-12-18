#* 9: Render Manuscript Figures
#+ 9.1: Figure 1
{
fig1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # Panels
  draw_plot(ggdraw() + draw_grob(p1A), x = 2.25, y = 7.303333, width = 4, height = 2.78) +
  draw_plot(p1B, x = 1.18, y = 4.133333, width = 2.8, height = 2.8) +
  draw_plot(p1C, x = 4.35, y = 4.133333, width = 2.8, height = 2.8) +
  # Labels
  figure_labels(list(
  A = c(2.296667, 9.98),
  B = c(1.276667, 7.058333),
  C = c(4.446667, 7.058333),
  "Figure 1" = c(0.49, 10.43)
))
print_to_png(fig1+grdgd(), "fig1.png", dpi = 300)
}
#+ 9.3: Figure 3
{
fig3 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # Panels
  draw_plot(p3A, x = 0.755, y = 6.646667, width = 3.5, height = 3.5) +
  draw_plot(p3B, x = 4.275, y = 6.646667, width = 3.5, height = 3.5) +
  draw_plot(p3C, x = 0.451667, y = 1.786667, width = 3.6, height = 5) +
  draw_plot(p3D, x = 4.066667-131/600+23.5/300, y = 1.786667, width = 3.6, height = 5) +
  # Labels
  figure_labels(list(
  A = c(0.815, 9.98),
  B = c(4.335, 9.98),
  C = c(0.815, 6.441667),
  D = c(4.335, 6.441667),
  "Figure 3" = c(0.49, 10.43)
))
print_to_png(fig3+grdgd(), "fig3.png", dpi = 300)
}
