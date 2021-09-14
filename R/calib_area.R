##' Delimit calibration area for constructing species distribution models
#'
#' @description This function offers different methods to define the calibration area. The output could be used with other flexsdm functions like sample_backgroud, sample_pseudoabs, and sdm_predict, among others
#'
#' @param data data.frame or tibble. Database with presences
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param method character. Method used for delimiting a calibration area. Could be necessary to concatenate (c()) different objects for this argument. The following methods are implemented:
#' \itemize{
#' \item buffer: calibration area is defined by a buffer around presences. Usage method = c('buffer', width=40000).
#' \item mcp: calibration area is defined by a minimum convex polygon. Usage method = 'mcp'.
#' \item bmcp: calibration area is defined by buffered minimum convex polygon with buffer width. Usage method = c('bmcp', width=40000).
#' \item mask: calibration area is defined by selected polygons in a spatial vector object intersected by presences. Usage method = c("mask", clusters, "DN"). The second concatenated element must be a SpatVector, the third element is a character with the column name from SpatVector used for filtering polygons.
#' }
#' @param groups character. Column name indicating differentiated subsets of points. This could be used with mcp and bmcp method. Default NULL
#' @param crs character. Coordinate reference system used for transforming occurrences and outputs. If set as NULL, crs of result will be NA for buffer, mcp, and bmcp methods. For mask method, the result will have the same crs as the SpatVector used
#'
#' @return
#' A SpatVector
#' @export
#'
#' @importFrom grDevices chull
#' @importFrom methods as
#' @importFrom raster buffer
#' @importFrom terra vect union crs extract
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
  . <- NULL
  if (!method[1] %in% c("buffer", "mcp", "bmcp", "mask")) {
    stop("argument 'method' was misused, available methods buffer, mcp, bmpc, and mask")
  }

  if (method[1] %in% c("bmcp", "buffer")) {
    if (!"width" %in% names(method)) {
      stop("provide width value for ", method[1], " method", ", e.g. method = c('bmcp', width = 70000)")
    }
  }

  if (method[1] %in% c("mask")) {
    if(is.na(method[2])){
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
    if (is.null(crs)) {
      data_sp <- terra::vect(data_sp, geom = names(data_sp))
      result <- raster::buffer(methods::as(data_sp, "Spatial"), width = as.numeric(method["width"]))
      result <- terra::vect(result)
    } else {
      data_sp <- terra::vect(data_sp, geom = names(data_sp), crs = crs)
      result <- raster::buffer(methods::as(data_sp, "Spatial"), width = as.numeric(method["width"]))
      result <- terra::vect(result)
    }
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
      data_pl <- apply(data_pl, 1, function(x) {
        paste(x[1], x[2])
      }) %>%
        paste(., collapse = ", ") %>%
        paste("POLYGON", "((", ., "))")
      if (is.null(crs)) {
        data_pl <- terra::vect(data_pl)
      } else {
        data_pl <- terra::vect(data_pl, crs = crs)
      }
      result[[i]] <- data_pl
    }
    if (length(result) > 1) {
      result <- do.call(
        terra::union,
        sapply(result, function(x) {
          methods::as(x, "Spatial")
        })
      ) %>%
        terra::vect()
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
      data_pl <- apply(data_pl, 1, function(x) {
        paste(x[1], x[2])
      }) %>%
        paste(., collapse = ", ") %>%
        paste("POLYGON", "((", ., "))")
      if (is.null(crs)) {
        data_pl <- terra::vect(data_pl)
      } else {
        data_pl <- terra::vect(data_pl, crs = crs)
      }
      data_pl <- raster::buffer(methods::as(data_pl, "Spatial"), width = as.numeric(method["width"]))
      result[[i]] <- data_pl
    }
    if (length(result) > 1) {
      result <- sapply(result, terra::vect)
      result <- do.call(terra::union, result)
    } else {
      result <- terra::vect(result[[1]])
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
