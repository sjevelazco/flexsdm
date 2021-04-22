#' Delimit calibration area for ecological niche models
#'
#' @description This function offers different methods to define de calibration area. The output could be used with other flexsdm functions like sample_backgroud, sample_pseudoabs, enm_predict, and enm_ensemble, among others
#'
#' @param data data.frame or tibble. Database with presences
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param method character. Pseudo-absence allocation method. Could be necessary concatenate (c()) different obje for this argument. The next methods are implemented:
#' \itemize{
#' \item buffer: calibration area is defined by a buffer around presences. Usage method = c('buffer', width=40000).
#' \item mcp: calibration area is defined by a minimum convex polygon. Usage method = 'mcp'.
#' \item bmcp: calibration area is defined by buffed minimum convex polygon. Usage method = c('bmcp', width=40000).
#' \item mask: calibration area is defined by those polygons intersected by presences. Usage method = c("mask", clusters, "DN"). The second element concatenated must be a SpatialPolygonDataFrame
#' }
#' @param groups character. Column name with that differentiate set of points. This could be used with mcp and bmcp method. Default NULL
#' @param crs character. Coordinate reference system used for transforming occurrences and outputs. In case it is set as NULL, crs of result will be NA for buffer, mcp, and bmcp methods. For mask method, the result will have the same crs as SpatialPolygonDataFrame used
#'
#' @return
#' A SpatialPolygon or SpatialPolygonDataFrame
#' @export
#'
#' @importFrom grDevices chull
#' @importFrom raster buffer bind projection extract
#' @importFrom sp coordinates
#'
#' @examples
#' \dontrun{
#' require(raster)
#' require(dplyr)
#' data("spp")
#' data("clusters")
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
#' plot(single_spp[, 2:3], pch = 19)
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
#'   method = c("mask", clusters, "DN"),
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
    if (!"SpatialPolygonsDataFrame" %in% class(clusters)) {
      stop("provide a SpatialPolygonDataFrame in method argument", ", e.g. method = c('mask', clusters)")
    }
  }

  data <- data.frame(data[, c(x, y, groups)])
  names(data) <- c("x", "y", "groups")[!sapply(list(x, y, groups), is.null)]


  if (method[1] == "buffer") {
    data <- data[, c("x", "y")]
    data_sp <- data
    sp::coordinates(data_sp) <- ~ x + y
    result <- raster::buffer(data_sp, width = as.numeric(method["width"]))
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
      result[[i]] <- data_pl
    }
    if (length(result) > 1) {
      result <- do.call(raster::bind, result)
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
      data_pl <- raster::buffer(data_pl, width = as.numeric(method["width"]))
      result[[i]] <- data_pl
    }
    if (length(result) > 1) {
      result <- do.call(raster::bind, result)
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
    raster::projection(data_sp) <- raster::crs(polyc)
    data_sp <- sp::spTransform(data_sp, raster::crs(polyc))
    result <- raster::extract(polyc, data_sp)[, cname]
    result <- polyc[polyc[[cname]] %in% result, ]
  }
  return(result)
}
