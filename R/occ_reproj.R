#' Reproject location coordinates to a new coordinate reference system
#'
#' @description
#' Reprojects the x/y (longitude/latitude) coordinates of a set of species
#' occurrence records from one coordinate reference system (CRS) to another,
#' and returns the results either as new columns appended to the input data
#' or as replacements for the original coordinate columns.
#'
#' @param data data.frame or tibble. A table containing occurrence records
#' with at least two columns for the x (longitude) and y (latitude)
#' coordinates.
#' @param x character. Name of the column in \code{data} containing the
#' x (longitude) coordinate. Default \code{"x"}.
#' @param y character. Name of the column in \code{data} containing the
#' y (latitude) coordinate. Default \code{"y"}.
#' @param from character. CRS of the input coordinates, provided as an EPSG
#' code (e.g., \code{"EPSG:4326"}) or any other CRS string accepted by
#' \code{terra::vect()}. Default \code{"EPSG:4326"} (WGS84 unprojected).
#' @param to character. CRS to project the coordinates to, provided as an
#' EPSG code (e.g., \code{"EPSG:32617"}) or any other CRS string accepted
#' by \code{terra::project()}. This argument is mandatory; the function
#' will stop with an error if \code{to = NULL}.
#' @param replace_cols logical. If \code{TRUE}, the original \code{x} and
#' \code{y} columns in \code{data} are overwritten with the reprojected
#' coordinates. If \code{FALSE} (default), the reprojected coordinates are
#' appended as new columns \code{x_new} and \code{y_new}, leaving the
#' original columns unchanged.
#'
#' @return
#' A tibble (or data.frame, if \code{replace_cols = TRUE}) with the
#' reprojected coordinates. If \code{replace_cols = FALSE}, two new columns,
#' \code{x_new} and \code{y_new}, are added to \code{data}. If
#' \code{replace_cols = TRUE}, the original \code{x} and \code{y} columns
#' are replaced in place and the object class of \code{data} is preserved.
#'
#' @details
#' This function is a convenience wrapper around \code{terra::vect()} and
#' \code{terra::project()} for reprojecting tabular occurrence coordinates
#' without needing to first convert \code{data} into a spatial object.
#' It is useful, for example, when occurrence records are stored in
#' geographic coordinates (e.g., \code{"EPSG:4326"}) but need to be
#' converted to a projected CRS (e.g., a UTM zone) for distance-based
#' calculations or to match the CRS of environmental raster layers used
#' elsewhere in a species distribution modeling workflow.
#'
#' @seealso \code{\link[terra]{project}}, \code{\link[terra]{vect}}
#'
#' @export
#'
#' @importFrom terra vect project geom
#' @importFrom dplyr select as_tibble
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   species = c("sp1", "sp2", "sp3"),
#'   x = c(-74.1, -73.9, -74.0),
#'   y = c(4.6, 4.7, 4.65)
#' )
#'
#' # Append reprojected coordinates as new columns
#' occ_reproj(data, x = "x", y = "y", from = "EPSG:4326", to = "EPSG:32618")
#'
#' # Replace original coordinate columns with reprojected values
#' occ_reproj(
#'   data,
#'   x = "x", y = "y",
#'   from = "EPSG:4326", to = "EPSG:32618",
#'   replace_cols = TRUE
#' )
#' }
occ_reproj <- function(data, x = "x", y = "y", from = "EPSG:4326", to = NULL, replace_cols = FALSE) {
  if (is.null(to)) {
    stop("A projection must be specified in 'to' argument")
  } else {
    v <- data[c(x, y)]
    names(v) <- c("lon", "lat")
    v <- terra::vect(v, crs = from)
    v <- terra::project(v, to)
    v <- terra::geom(v) %>%
      as.data.frame() %>%
      dplyr::select(x, y)
    names(v) <- c("x_new", "y_new")
    if (replace_cols) {
      data[c(x, y)] <- v
      return(data)
    } else {
      return(cbind(data, v) %>% dplyr::as_tibble())
    }
  }
}
