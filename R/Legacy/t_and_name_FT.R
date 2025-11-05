t_and_name_FT <- function(data_tibble, sequence_tibble) {
  data_tibble %>%
    pivot_longer(cols = -Feature, names_to = "File_Name", values_to = "Value") %>%
    pivot_wider(names_from = Feature, values_from = Value) %>%
    left_join(sequence_tibble, by = "File_Name") %>% # Join based on File_Name
    select(Sample_ID, everything(), -File_Name) %>% # Make Sample_ID the first column and remove File_Name
    arrange(Sample_ID)
}