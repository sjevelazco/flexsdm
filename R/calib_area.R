##' Delimit calibration area for constructing species distribution models
#'
#' @description This function offers different methods to define the calibration area. The output could be used with other flexsdm functions like sample_backgroud, sample_pseudoabs, and sdm_predict, among others
#'
#' @param data data.frame or tibble. Database with presences
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param method character. Method used for delimiting a calibration area. Could be necessary to concatenate (c()) different objects for this argument. The following methods are implemented:
#' \itemize{
#' \item buffer: calibration area is defined by a buffer around presences. Usage method = c('buffer', width=40000). A value of the buffer width in m must be provided if CRS has a longitude/latitude, or in map units in other cases
#' \item mcp: calibration area is defined by a minimum convex polygon. Usage method = 'mcp'.
#' \item bmcp: calibration area is defined by buffered minimum convex polygon with buffer width. Usage method = c('bmcp', width=40000). A value of the buffer width in m must be provided if CRS has a longitude/latitude, or in map units in other cases
#' \item mask: calibration area is defined by selected polygons in a spatial vector object intersected by presences. Usage method = c("mask", clusters, "DN"). The second concatenated element must be a SpatVector, the third element is a character with the column name from SpatVector used for filtering polygons.
#' }
#' @param groups character. Column name indicating differentiated subsets of points. This could be used with mcp and bmcp method. Default NULL
#' @param crs character. Coordinate reference system used for transforming occurrences and outputs. If set as NULL, the result of mask method will have the same crs as the SpatVector used. Define a crs  is mandatory for buffer, mcp
#' and bmcp method.
#'
#' @return
#' A SpatVector
#' @export
#'
#' @importFrom grDevices chull
#' @importFrom methods as
#' @importFrom terra vect buffer aggregate union crs extract
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#' data("spp")
#' clusters <- system.file("external/clusters.shp", package = "flexsdm")
#' clusters <- terra::vect(clusters)
#'
#'
#' single_spp <-
#'   spp %>%
#'   dplyr::filter(species == "sp1") %>%
#'   dplyr::filter(pr_ab == 1) %>%
#'   dplyr::select(-pr_ab)
#'
#'
#' plot(clusters)
#' points(single_spp[-1], col="red")
#' crs(clusters, proj=TRUE) # coordinate reference system (CRS) used for this points database
#' # note that the unit of this CRS is in m, consequently the buffer width
#' # will be interpreted in m too
#'
#' # buffer method
#' ca_1 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = c("buffer", width = 40000),
#'   crs = crs(clusters)
#' )
#' plot(ca_1)
#' points(single_spp[, 2:3], pch = 19, cex = 0.5)
#'
#' # mcp method
#' ca_2 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = "mcp",
#'   crs = crs(clusters)
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
#'   crs = crs(clusters),
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
#'   method = c("bmcp", width = 30000),
#'   crs = crs(clusters)
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
#'   crs = crs(clusters),
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
  . <- NULL
  if (!method[1] %in% c("buffer", "mcp", "bmcp", "mask")) {
    stop("argument 'method' was misused, available methods buffer, mcp, bmpc, and mask")
  }

  if (method[1] %in% c("bmcp", "buffer")) {
    if (!"width" %in% names(method)) {
      stop("provide width value for ", method[1], " method", ", e.g. method = c('bmcp', width = 70000)")
    }

    if (method[1] %in% c("bmcp", "buffer", "mcp") & is.null(crs)) {
      stop("A coordinate reference system is needed in 'crs' agument for this method")
    }
  }


  if (method[1] %in% c("mask")) {
    if (is.na(method[2])) {
      stop("provide a SpatVector or SpatialPolygonDataFrame in method argument", ", e.g. method = c('mask', clusters)")
    }
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
    data_sp <- terra::vect(data_sp, geom = names(data_sp), crs = crs)
    result <- terra::buffer(data_sp, width = as.numeric(method["width"])) %>%
      terra::aggregate()
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
      data_pl <- data.frame(object = 1, part = 1, data_pl, hole = 0)
      data_pl <- terra::vect(as.matrix(data_pl), type = "polygons", crs = crs)
      result[[i]] <- data_pl
    }
    if (length(result) > 1) {
      result <- do.call(terra::union, result)
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
      data_pl <- data.frame(object = 1, part = 1, data_pl, hole = 0)
      data_pl <- terra::vect(as.matrix(data_pl), type = "polygons", crs = crs)
      data_pl <- terra::buffer(data_pl, width = as.numeric(method["width"]))
      result[[i]] <- data_pl
    }
    if (length(result) > 1) {
      result <- do.call(terra::union, result)
    } else {
      result <- result[[1]]
    }
  }

  if (method[1] == "mask") {
    polyc <- method[[2]]
    cname <- method[[3]]
    data <- data[, c("x", "y")]
    data_sp <- data
    data_sp <- terra::vect(data_sp, geom = names(data_sp), crs = terra::crs(polyc))
    result <- terra::extract(polyc, data_sp)[, cname] %>% unique()
    result <- polyc[polyc[[cname]][, 1] %in% result, ]
  }
  return(result)
}
