#' Centralized plot color themes and styling
#'
#' Returns standardized color palettes and theme elements for consistent
#' visualization across all plots in the project.
#'
#' @param theme_name Character; name of the theme to return. Options:
#'   - "pgd_status": Colors for PGD status (Y/N, PGD/No PGD)
#'   - "time_points": Colors for time point comparisons
#'   - "all": Returns a list of all available themes
#'
#' @details
#' The function returns named vectors where names correspond to the data values:
#' - For PGD status: names are "Y" (red) and "N" (blue)
#' - For display labels: use the returned vector with custom names
#'
#' @return Named vector of colors or list of all themes
#' @export
#'
#' @examples
#' # Get PGD status colors (for data with Y/N values)
#' pgd_colors <- plot_themes("pgd_status")
#'
#' # Use with custom labels in ggplot
#' color_mapping <- c("Y" = plot_themes("pgd_status")["Y"],
#'                    "N" = plot_themes("pgd_status")["N"])
plot_themes <- function(theme_name = "pgd_status") {
  # Define all available themes
  themes <- list(
    pgd_status = c(
      "Y" = "#94001E",  # Red for Severe PGD
      "N" = "#03507D"   # Blue for No Severe PGD
    ),
    pgd_status_labels = c(
      "PGD" = "#94001E",      # Red for Severe PGD
      "No PGD" = "#03507D"    # Blue for No Severe PGD
    ),
    time_points = c(
      "12" = "#1b9e77",
      "24" = "#d95f02",
      "Combined" = "#7570b3"
    )
  )

  if (theme_name == "all") {
    return(themes)
  }

  if (!theme_name %in% names(themes)) {
    stop("Unknown theme: '", theme_name, "'. Available themes: ",
         paste(names(themes), collapse = ", "))
  }

  themes[[theme_name]]
}
