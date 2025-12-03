#' Read and filter clinical metadata sheets
#'
#' @param sheet_name Name of the Excel sheet to read
#' @param patients Vector of patient IDs to filter by (optional)
#' @param suppress_warnings Whether to suppress warnings (default: FALSE)
#' @return Tibble with clinical data, optionally filtered by patients
#' @export
read_clinical_sheet <- function(sheet_name, patients = NULL, suppress_warnings = FALSE) {
  read_fn <- if (suppress_warnings) {
    function(...) suppressWarnings(readxl::read_xlsx(...))
  } else {
    readxl::read_xlsx
  }

  result <- read_fn(
    config$paths$clinical_metadata,
    sheet = sheet_name,
    na = c("", "NA", "-")
  )

  # Filter by patients if provided
  if (!is.null(patients)) {
    result <- result |> dplyr::filter(Patient %in% patients)
  }

  return(result)
}