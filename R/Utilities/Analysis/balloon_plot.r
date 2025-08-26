balloon_plot <- function(balloon_dat){
 # Create columns for magnitude (balloon size) and direction (balloon color)
 balloon_dat <- balloon_dat %>%
        filter(Label == TRUE) %>%
        select(Identified_Name, Value) %>%
        mutate(
          magn = abs(Value),
          direction = ifelse(Value > 0, "Up in PGD", "Down in PGD")
        ) %>%
        arrange(Value) %>%
        mutate(Identified_Name = factor(Identified_Name), levels = Identified_Name)
    # Create the balloon plot
    ball_plot <- ggplot(balloon_dat, aes(x = 1, y = reorder(Identified_Name, Value))) +
        geom_point(aes(size = magn, color = direction)) +
        scale_size(range = c(3, 10)) +  # adjust balloon size scaling
        scale_color_manual(values = c("Up in PGD" = "#800017", "Down in PGD" = "#113d6a")) +
        theme_minimal() +
        theme(
          axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank()
        ) +
          labs(
            title = "Top 10 Up and Down Fold Changes in PGD at 12h",
            x = 'Metabolite',
            y = 'Fold change',
            size = 'Magnitude',
            color = 'Direction'
          )
        list(
            balloon_dat = balloon_dat,
            ball_plot = ball_plot
        )
}
