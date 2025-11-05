process_feature_table <- function(data_tibble, sequence_tibble, data_type) {
  data_tibble %>%
    mutate(Feature = paste(data_type, mz, time, sep = "_")) %>%
    select(Feature, everything(), -c(mz, time)) %>%
    rename_with(~ gsub("\\.mzXML$", "", .)) %>%
    pivot_longer(cols = -Feature, names_to = "File_Name", values_to = "Value") %>%
    pivot_wider(names_from = Feature, values_from = Value) %>%
    left_join(sequence_tibble, by = "File_Name") %>%
    select(Sample_ID, everything(), -File_Name) %>%
    mutate(Sample_ID = sub("_.*", "", Sample_ID)) %>%
    arrange(Sample_ID)
}
