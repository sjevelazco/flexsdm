#' Delimit calibration area for constructing species distribution models
#'
#' @description This function offers different methods to define de calibration area. The output could be used with other flexsdm functions like sample_backgroud, sample_pseudoabs, enm_predict, and enm_ensemble, among others
#'
#' @param data data.frame or tibble. Database with presences
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param method character. Pseudo-absence allocation method. Could be necessary concatenate (c()) different objects for this argument. The next methods are implemented:
#' \itemize{
#' \item buffer: calibration area is defined by a buffer around presences. Usage method = c('buffer', width=40000).
#' \item mcp: calibration area is defined by a minimum convex polygon. Usage method = 'mcp'.
#' \item bmcp: calibration area is defined by buffed minimum convex polygon. Usage method = c('bmcp', width=40000).
#' \item mask: calibration area is defined by those polygons intersected by presences. Usage method = c("mask", clusters, "DN"). The second element concatenated must be a SpatialPolygonDataFrame,  the third element is a character with the column name from SpatialPolygonDataFrame used for filtering polygons.
#' }
#' @param groups character. Column name with that differentiate set of points. This could be used with mcp and bmcp method. Default NULL
#' @param crs character. Coordinate reference system used for transforming occurrences and outputs. In case it is set as NULL, crs of result will be NA for buffer, mcp, and bmcp methods. For mask method, the result will have the same crs as SpatialPolygonDataFrame used
#'
#' @return
#' A SpatialPolygon or SpatialPolygonDataFrame
#' @export
#'
#' @importFrom grDevices chull
#' @importFrom sp coordinates Polygon Polygons SpatialPolygons spTransform
#' @importFrom terra vect buffer crs extract
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#' data("spp")
#' clusters <- system.file("external/clusters.shp", package = "flexsdm")
#' clusters <- terra::vect(clusters)
#'
#' single_spp <-
#'   spp %>%
#'   dplyr::filter(species == "sp1") %>%
#'   dplyr::filter(pr_ab == 1) %>%
#'   dplyr::select(-pr_ab)
#'
#'
#' # buffer method
#' ca_1 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = c("buffer", width = 40000),
#' )
#' plot(ca_1)
#' points(single_spp[, 2:3], pch = 19, cex = 0.5)
#'
#' # mcp method
#' ca_2 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = "mcp"
#' )
#' plot(ca_2)
#' points(single_spp[, 2:3], pch = 19, cex = 0.5)
#'
#' # mcp method for different groups
#' single_spp <- single_spp %>% mutate(groups = ifelse(x > 150000, "a", "b"))
#'
#' plot(single_spp[, 2:3], pch = 19, col = "blue")
#' points(single_spp[single_spp$groups == "a", 2:3], col = "red", pch = 19)
#' points(single_spp[, 2:3])
#'
#' ca_2.1 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = c("mcp"),
#'   groups = "groups"
#' )
#' plot(ca_2.1)
#' points(single_spp[, 2:3], pch = 19, cex = 0.5)
#'
#' # bmcp method
#' ca_3 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = c("bmcp", width = 30000)
#' )
#' plot(ca_3)
#' points(single_spp[, 2:3], pch = 19, cex = 0.5)
#'
#' # bmcp method for different groups
#' ca_3.1 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = c("bmcp", width = 30000),
#'   groups = "groups"
#' )
#' plot(ca_3.1)
#' points(single_spp[, 2:3], pch = 19, cex = 0.5)
#'
#' # mask method
#' plot(clusters)
#' names(clusters)
#'
#' ca_3.1 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = c("mask", clusters, "clusters"),
#' )
#' plot(ca_3.1)
#' points(single_spp[, 2:3], pch = 19, cex = 0.5, col = "red")
#' }
calib_area <- function(data, x, y, method, groups = NULL, crs = NULL) {
  if (!method[1] %in% c("buffer", "mcp", "bmcp", "mask")) {
    stop("argument 'method' was misused, available methods buffer, mcp, bmpc, and maks")
  }

  if (method[1] %in% c("bmcp", "buffer")) {
    if (!"width" %in% names(method)) {
      stop("provide width value for ", method[1], " method", ", e.g. method = c('bmcp', width = 70000)")
    }
  }

  if (method[1] %in% c("mask")) {
    if (!class(method[[2]]) %in% c("SpatialPolygonsDataFrame", "SpatVector")) {
      stop("provide a SpatVector or SpatialPolygonDataFrame in method argument", ", e.g. method = c('mask', clusters)")
    }
    if (class(method[[2]]) != "SpatVector") {
      method[[2]] <- terra::vect(method[[2]])
    }
  }

  data <- data.frame(data[, c(x, y, groups)])
  names(data) <- c("x", "y", "groups")[!sapply(list(x, y, groups), is.null)]


  if (method[1] == "buffer") {
    data <- data[, c("x", "y")]
    data_sp <- data
    sp::coordinates(data_sp) <- ~ x + y
    result <- terra::buffer(data_sp, width = as.numeric(method["width"]))
    result <- terra::vect(result)
  }

  if (method[1] == "mcp") {
    if (is.null(groups)) {
      data$groups <- 1
    }
    data <- split(data, data$groups)
    result <- list()
    for (i in 1:length(data)) {
      data_pl <- data.frame(data[[i]][, c("x", "y")])
      data_pl <- data_pl[grDevices::chull(data_pl), ]
      data_pl <- sp::Polygon(data_pl)
      data_pl <- sp::Polygons(list(data_pl), ID = 1)
      data_pl <- sp::SpatialPolygons(list(data_pl))
      result[[i]] <- terra::vect(data_pl)
    }
    if (length(result) > 1) {
      result <- do.call(rbind, result)
    } else {
      result <- result[[1]]
    }
  }

  if (method[1] == "bmcp") {
    if (is.null(groups)) {
      data$groups <- 1
    }
    data <- split(data, data$groups)
    result <- list()
    for (i in 1:length(data)) {
      data_pl <- data.frame(data[[i]][, c("x", "y")])
      data_pl <- data_pl[grDevices::chull(data_pl), ]
      data_pl <- sp::Polygon(data_pl)
      data_pl <- sp::Polygons(list(data_pl), ID = 1)
      data_pl <- sp::SpatialPolygons(list(data_pl))
      data_pl <- terra::buffer(data_pl, width = as.numeric(method["width"]))
      result[[i]] <- data_pl
    }
    if (length(result) > 1) {
      result <- sapply(result, terra::vect)
      result <- do.call(terra:::rbind.SpatVector, result)
    } else {
      result <- result[[1]]
    }
  }

  if (method[1] == "mask") {
    polyc <- method[[2]]
    cname <- method[[3]]
    data <- data[, c("x", "y")]
    data_sp <- data
    sp::coordinates(data_sp) <- ~ x + y
    terra::crs(data_sp) <- terra::crs(polyc)
    data_sp <- sp::spTransform(data_sp, terra::crs(polyc))
    result <- terra::extract(polyc, vect(data_sp))[, cname] %>% unique()
    result <- polyc[polyc[[cname]][, 1] %in% result, ]
  }
  return(result)
}
