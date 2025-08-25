undo_log2_transform <- function(df) {
      df <- 2^df # Apply the inverse transform to metabolite columns (4:end)
      return(df)
    }