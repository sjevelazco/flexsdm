#' Perform geographical filtering on species occurrences
#'
#' @description This function perform geographical filtering of species occurrences based on
#' different approach to define the minimum nearest-neighbor distance between points.
#'
#'
#' @param data data.frame. Data.frame or tibble object with presences
#' (or presence-absence) records, and coordinates
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param env_layer SpatRaster. Raster variables that will be used to fit the model
#' @param method character. Method to perform geographical thinning. Pairs of points are filtered
#' based on a geographical distance criteria.The following methods are available:
#' \itemize{
#'   \item moran: records are filtered based on the smallest distance that reduces Moran's I to
#'   values lower than 0.1. Latlong = TRUE if occurrences are in a geographical projection.
#'   Usage method: method = c('moran').
#'   \item cellsize: records are filtered based on the resolution of the environmental variables which
#'   can be aggregated to coarser resolution defined by the factor.
#'   Usage method: method = c('cellsize', factor = '2').
#'   \item defined: records are filtered based on a  distance value (d) provided in km.
#'   Usage method: method = c('defined', d = 300).
#' }
#' @param prj character. Projection string (PROJ4) for occurrences. Not necessary if
#' the projection used is WGS84 ("+proj=longlat +datum=WGS84"). Default "+proj=longlat +datum=WGS84"
#'
#' @return
#' A tibble object with data filtered geographically
#'
#' @details In this function three alternatives are implemented to determine the
#' distance threshold between pair of points:
#' 1-"moran" determines the minimum nearest-neighbor distance that minimizes the spatial
#' autocorrelation in occurrence data, following a Moran's semivariogram. A Principal Component
#' Analysis with the environmental variables is performed and then the first Principal Component is used
#' to calculate the semivariograms. Because of this, this method only allow the use of continuous variables.
#' Sometimes, this method can (too) greatly reduce the number of presences.
#' 2-"cellsize" filters occurrences based on the predictors' resolution. This method will calculate
#' the distance between the first two cells of the environmental variable and use this distance
#' as minimum nearest-neighbor distance to filter occurrences.
#' The resolution of the raster is aggregated based on the values used in "factor". Thus, the distance
#' used for filtering can be adjusted to represent a larger grid size.
#' 3-"determined" this method uses any minimum nearest-neighbor distance specified in km.
#'
#'
#' For the third method the "thin" function from spThin package is used
#' (Aiello-Lammens et al., 2015) with the following argument settings reps = 20, write.files = FALSE,
#' locs.thinned.list.return = TRUE, and write.log.file = FALSE.
#'
#'
#' @references
#' \itemize{
#' \item Aiello-Lammens, M. E., Boria, R. A., Radosavljevic, A., Vilela, B., & Anderson,
#' R. P. (2015). spThin: An R package for spatial thinning of species occurrence records for use
#' in ecological niche models. Ecography, 38(5), 541-545. https://doi.org/10.1111/ecog.01132
#' }
#'
#' @export
#'
#' @importFrom dplyr tibble
#' @importFrom spThin thin
#' @importFrom stats complete.cases prcomp
#' @importFrom terra vect project geom extract as.data.frame predict distance xyFromCell
#' @importFrom utils capture.output
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#'
#' # Environmental variables
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' plot(somevar)
#'
#' # Species occurrences
#' data("spp")
#' spp
#' spp1 <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)
#'
#' somevar[[1]] %>% plot()
#' points(spp1 %>% select(x, y))
#'
#' # Using Moran method
#' filtered_1 <- occfilt_geo(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   env_layer = somevar,
#'   method = c("moran"),
#'   prj = crs(somevar)
#' )
#'
#' somevar[[1]] %>% plot(col = gray.colors(10))
#' points(spp1 %>% select(x, y)) # raw data
#' points(filtered_1 %>% select(x, y), pch = 19, col = "yellow") # filtered data
#'
#' # Using cellsize method
#' filtered_2 <- occfilt_geo(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   env_layer = somevar,
#'   method = c("cellsize", factor = "3"),
#'   prj = crs(somevar)
#' )
#'
#' somevar[[1]] %>% plot(col = gray.colors(10))
#' points(spp1 %>% select(x, y)) # raw data
#' points(filtered_2 %>% select(x, y), pch = 19, col = "yellow") # filtered data
#'
#'
#' # Using defined method
#' filtered_3 <- occfilt_geo(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   env_layer = somevar,
#'   method = c("defined", d = "30"),
#'   prj = crs(somevar)
#' )
#'
#' somevar[[1]] %>% plot(col = gray.colors(10))
#' points(spp1 %>% select(x, y)) # raw data
#' points(filtered_3 %>% select(x, y), pch = 19, col = "yellow") # filtered data
#' }
#'
#' @seealso \code{\link{occfilt_env}}
#'
occfilt_geo <- function(data, x, y, env_layer, method, prj = "+proj=longlat +datum=WGS84") {
  da <- data0 <- data
  rm(data)
  da <- da[c(x, y)]


  if (prj != "+proj=longlat +datum=WGS84") {
    da <- terra::vect(data, geom = c(x, y), prj)
    da <- terra::project(da, "+proj=longlat +datum=WGS84")
    da <- data.frame(terra::geom(da))
    names(da)[names(da) %in% c("x", "y")] <- c(x, y)
    da <- dplyr::tibble(cbind(da[c(x, y)]))
    if (!("moran" %in% method)) {
      env_layer <- env_layer[[1]]
    }
    env_layer <-
      terra::project(env_layer, "+proj=longlat +datum=WGS84")
  }

  # Remove NAs
  da <- da[c(x, y)]
  coord <- da[c(x, y)]

  message("Extracting values from raster ... ")
  env <- terra::extract(env_layer, coord)
  env$ID <- NULL

  filt <- stats::complete.cases(env)
  if (sum(!filt) > 0) {
    message(sum(!filt), " records were removed because they have NAs for some variables")
    da <- da[filt, ]
    coord <- coord[filt, ]
    env <- env[filt, ]
  }
  rm(filt)

  message("Number of unfiltered records: ", nrow(da))


  if ("moran" %in% method) {
    # Perform PCA
    p <- terra::as.data.frame(env_layer, xy = FALSE, na.rm = TRUE)

    p <- stats::prcomp(p,
      retx = TRUE,
      scale. = TRUE,
      center = TRUE,
      rank. = 1
    )

    env_layer <- terra::predict(env_layer, p)
    names(env_layer) <- "PC1"

    # Extract occurrence PC1 values
    vals <- terra::extract(env_layer, coord)$PC1
    # Trasform coord to vect
    coord <- vect(coord, geom = c(x, y), crs = "+proj=longlat +datum=WGS84")

    # Calculate distances & Define breaks
    dists <- terra::distance(coord, coord) / 1000 # dist. in km
    br <- max(dists) / 30
    br <- seq(from = br, to = max(dists), by = br)

    # Calculate Moran for breaks
    mor <- rep(NA, length(br))
    for (i in 1:length(br)) {
      dists2 <- dists
      dists2[dists2 > br[i]] <- 0
      dists2 <- as.matrix(dists2)
      diag(dists2) <- 0
      mor[i] <- abs(morani(x = vals, weight = dists2, scaled = TRUE))
    }

    if (any(mor <= 0.1)) {
      message("Threshold for Moran: ", 0.1)
      pos <- which(mor <= 0.1)[1]
      d <- br[pos]
    } else {
      message("Threshold for Moran: ", round(min(mor), 2))
      pos <- which(mor == min(mor))
      d <- br[pos]
    }

    options(warn = 1)

    # Thinning
    da$.spp <- "sp"
    invisible(utils::capture.output(
      occT <-
        spThin::thin(
          loc.data = da,
          lat.col = y,
          long.col = x,
          spec.col = ".spp",
          thin.par = d,
          reps = 20,
          write.files = FALSE,
          locs.thinned.list.return = TRUE,
          write.log.file = FALSE
        )
    ))
    occT <-
      occT[[which(sapply(occT, function(x) {
        nrow(x)
      }) == max(sapply(occT, function(x) {
        nrow(x)
      })))[1]]]
    occT <- as.integer(row.names(occT))

    # Select Thinned Occurrences
    coord_filter <- data0[occT, ]

    # Results
    message("Distance threshold(km): ", round(br[pos], 3))
    message("Number of filtered records: ", nrow(coord_filter))
  }

  if ("cellsize" %in% method) {
    # Haversine Transformation & Distance Threshold
    factor <- as.numeric(method["factor"])
    distance <-
      terra::xyFromCell(env_layer[[1]], 1:2)
    distance <- as.numeric(abs(terra::distance(distance, distance, lonlat = TRUE)[2] / 1000 * factor))

    # Thinning
    da$.spp <- "sp"
    invisible(utils::capture.output(
      occT <-
        spThin::thin(
          loc.data = da,
          lat.col = y,
          long.col = x,
          spec.col = ".spp",
          thin.par = distance,
          reps = 20,
          write.files = FALSE,
          locs.thinned.list.return = TRUE,
          write.log.file = FALSE
        )
    ))
    occT <-
      occT[[which(sapply(occT, function(x) {
        nrow(x)
      }) == max(sapply(occT, function(x) {
        nrow(x)
      })))[1]]]
    occT <- as.integer(row.names(occT))

    # Select Thinned Occurrences
    coord_filter <- data0[occT, ]
    message("Distance threshold(km): ", round(distance, 3))
    message("Number of filtered records: ", nrow(coord_filter))
  }

  if ("defined" %in% method) {
    # Thinning
    da$.spp <- "sp"
    invisible(utils::capture.output(
      occT <-
        spThin::thin(
          loc.data = da,
          lat.col = y,
          long.col = x,
          spec.col = ".spp",
          thin.par = as.numeric(method["d"]),
          reps = 20,
          write.files = FALSE,
          locs.thinned.list.return = TRUE,
          write.log.file = FALSE
        )
    ))
    occT <-
      occT[[which(sapply(occT, function(x) {
        nrow(x)
      }) == max(sapply(occT, function(x) {
        nrow(x)
      })))[1]]]
    occT <- as.integer(row.names(occT))

    # Select Thinned Occurrences
    coord_filter <- data0[occT, ]
    message("Distance threshold(km): ", round(as.numeric(method["d"]), 3))
    message("Number of filtered records: ", nrow(coord_filter))
  }
  return(dplyr::tibble(coord_filter))
}
