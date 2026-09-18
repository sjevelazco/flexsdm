#' Build a presence-absence database from multi-species occurrence data
#'
#' @description
#' Converts a database containing occurrence records for multiple species
#' into a presence-absence database, either for every species present in
#' the data or for a specified subset of target species. For each target
#' species, its own records are coded as presences and the records of all
#' other species are recoded as absences for that species.
#'
#' @param data data.frame or tibble. A table with occurrence records for
#' one or more species, including columns for species identity and x/y
#' coordinates.
#' @param x character. Name of the column in \code{data} with the
#' x (longitude) coordinate. Default \code{"x"}.
#' @param y character. Name of the column in \code{data} with the
#' y (latitude) coordinate. Default \code{"y"}.
#' @param species character. Name of the column in \code{data} with the
#' species name/identifier. Default \code{"species"}.
#' @param target_species character vector. One or more species names to
#' build presence-absence data for. If \code{NULL} (default), a
#' presence-absence database is built for every unique species found in
#' \code{data[[species]]}.
#' @param pr_ab_name character. Name to give the presence-absence column
#' in the output, coded \code{1} for presences and \code{0} for absences.
#' Default \code{"pr_ab"}.
#'
#' @return
#' A tibble with columns \code{species}, \code{x}, \code{y} (using the
#' names supplied in the \code{species}, \code{x}, and \code{y} arguments)
#' and the presence-absence column named according to \code{pr_ab_name}.
#' For each species in \code{target_species} (or each unique species in
#' \code{data}, if \code{target_species = NULL}), the output contains one
#' presence-absence block: records belonging to that species are coded
#' \code{1}, and records of all other species are recoded as belonging to
#' that species with a value of \code{0}. Blocks for all target species are
#' stacked row-wise, so the returned tibble has
#' \code{nrow(data) * length(sp)} rows, where \code{sp} is the set of
#' target species.
#'
#' @details
#' Only the \code{species}, \code{x}, and \code{y} columns are retained in
#' the output; any other columns present in \code{data} (e.g.,
#' environmental covariates or metadata) are dropped. If you need those
#' columns downstream, re-join them after calling this function.
#'
#' Because absences for a given target species are drawn from every other
#' species' occurrence records, duplicate x/y coordinates across species
#' in the input will produce duplicate absence coordinates in the output.
#'
#' @export
#'
#' @importFrom dplyr bind_rows as_tibble
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   species = c("sp1", "sp1", "sp2", "sp2", "sp3"),
#'   x = c(-74.1, -74.2, -73.9, -73.8, -74.0),
#'   y = c(4.6, 4.65, 4.7, 4.72, 4.68)
#' )
#'
#' # Presence-absence database for every species in the data
#' get_absences(data, x = "x", y = "y", species = "species", target_species = NULL, pr_ab_name = "pr_ab")
#'
#' # Presence-absence database for a single target species
#' get_absences(data, x = "x", y = "y", species = "species", target_species = "sp1", pr_ab_name = "pr_ab")
#'
#' # Presence-absence database for a subset of species, with a custom
#' # presence-absence column name
#' get_absences(data, x = "x", y = "y", species = "species", target_species = c("sp1", "sp2"), pr_ab_name = "occ")
#' }
get_absences <- function(data, x = "x", y = "y", species = "species",
                         target_species = NULL, pr_ab_name = "pr_ab") {
  if (!all(c(x, y, species) %in% names(data))) {
    stop("'data' must contain the columns specified in 'x', 'y', and 'species'.")
  }

  data <- as.data.frame(data)
  data <- data[c(species, x, y)]
  data[[pr_ab_name]] <- 1

  if (is.null(target_species)) {
    sp <- sort(unique(data[[species]]))
  } else {
    sp <- target_species
  }

  data_new <- vector("list", length(sp))
  for (i in seq_along(sp)) {
    pr <- data[data[[species]] == sp[i], ]
    ab <- data[data[[species]] != sp[i], ]
    ab[[pr_ab_name]] <- 0
    ab[[species]] <- sp[i]
    data_new[[i]] <- dplyr::bind_rows(pr, ab)
  }

  data_new <- dplyr::bind_rows(data_new)
  return(dplyr::as_tibble(data_new))
}