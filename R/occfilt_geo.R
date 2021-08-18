#' Perform geographical filtering on species occurrences
#'
#' @param data data.frame. Data.frame or tibble object with presences
#' (or presence-absence) records, and coordinates
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param id character. Column names with rows id. It is important that each row has its own unique code.
#' @param env_layer SpatRaster. Raster variables that will be used to fit the model
#' @param method character. Method to perform geographical thinning. Pair of points are filtered based on a geographical distance criteria.The following methods are available:
#' \itemize{
#'   \item moran: Points are filetered based on the smallest distance that reduces Moran's I to values lower than 0.1. Latlong = T if occurrences are in a gegraphical projection.
#'   Usage method: method = c('moran').
#'   \item cellsize: Points are filetered based on the environmental variables resolution. The factor can be specified.
#'   Usage method: method = c('cellsize', factor = '2').
#'   \item defined: Points are filetered based on a provided distance value in km.
#'   Usage method: method = c('defined', d = 300).
#' }
#' @param prj character. Projection string (PROJ4) for occurrences. Not necessary if the projection is WGS84.
#'
#' @return
#' A tibble object with filtered data
#'
#' @export
#'
#' @importFrom ape Moran.I
#' @importFrom dplyr tibble
#' @importFrom spThin thin
#' @importFrom stats complete.cases prcomp dist
#' @importFrom terra vect project geom extract as.data.frame predict xyFromCell distance
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
#' # Moran
#' filtered_1 <- occfilt_geo(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   id = "species",
#'   env_layer = somevar,
#'   method = c("moran"),
#'   prj = crs(somevar, proj = T)
#' )
#'
#' # Cellsize
#' filtered_2 <- occfilt_geo(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   id = "species",
#'   env_layer = somevar,
#'   method = c("cellsize", factor = "3"),
#'   prj = crs(somevar, proj = T)
#' )
#'
#' # Defined
#' filtered_3 <- occfilt_geo(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   id = "species",
#'   env_layer = somevar,
#'   method = c("defined", d = "30"),
#'   prj = crs(somevar, proj = T)
#' )
#' }
occfilt_geo <- function(data, x, y, id, env_layer, method, prj = NULL) {
  data <- data[c(x, y, id)]

  if (!is.null(prj) || prj != "+proj=longlat +datum=WGS84") {
    da <- terra::vect(data, geom = c(x, y), prj)
    da <- terra::project(da, "+proj=longlat +datum=WGS84")
    da <- data.frame(terra::geom(da))
    da <- dplyr::tibble(cbind(da[c(x, y)], data[c(id)]))

    env_layer <- terra::project(env_layer, "+proj=longlat +datum=WGS84")
  } else {
    da <- data
  }

  # Remove NAs
  da <- da[c(x, y, id)]
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

  if ("moran" %in% method) {
    if (!all(grepl("PC", names(env_layer)))) {
      p <- terra::as.data.frame(env_layer, xy = FALSE, na.rm = TRUE)

      p <- stats::prcomp(p,
        retx = TRUE,
        scale. = TRUE,
        center = TRUE
      )

      env_layer <- terra::predict(env_layer, p, index = 1)
    } else {
      env_layer <- env_layer[[1]]
    }
    # Extract occurrence PC1 values
    vals <- terra::extract(env_layer, coord)$PC1

    # Convert from decimals to km
    lat2grd <- function(input) {
      toradians <- atan(1) / 45
      radiusearth <- 0.5 * (6378.2 + 6356.7)
      sine51 <- sin(51.5 * toradians)
      output <- data.frame(cbind(
        x = (input[, 1] * toradians) * radiusearth * sine51,
        y = (input[, 2] * toradians) * radiusearth
      ))
    }
    coord <- lat2grd(coord)

    # Calculate distances & Define breaks
    dists <- stats::dist(coord)
    br <- max(dists) / 20
    br <- seq(from = br, to = max(dists), by = br)

    # Calculate Moran for breaks
    mor <- rep(NA, length(br))
    for (i in 1:length(br)) {
      dists2 <- dists
      dists2[dists2 > br[i]] <- 0
      dists2 <- as.matrix(dists2)
      diag(dists2) <- 0
      mor[i] <- abs(ape::Moran.I(x = vals, weight = dists2)$observed)
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
    invisible(utils::capture.output(
      occT <-
        spThin::thin(
          loc.data = da,
          lat.col = y,
          long.col = x,
          spec.col = id,
          thin.par = d,
          reps = 20,
          write.files = F,
          locs.thinned.list.return = T,
          write.log.file = F
        )
    ))
    occT <-
      occT[[which(sapply(occT, function(x) {
        nrow(x)
      }) == max(sapply(occT, function(x) {
        nrow(x)
      })))[1]]]
    occPOS <- as.integer(row.names(occT))

    # Select Thinned Occurrences
    coord_filter <- data[occPOS, ]

    # Results
    message("Distance threshold(km): ", round(br[pos], 3))
    message("Number of filtered records: ", nrow(coord_filter))
    return(dplyr::tibble(coord_filter))
  }

  if ("cellsize" %in% method) {

    # Haversine Transformation & Distance Threshold
    factor <- as.numeric(method["factor"])
    distance <-
      terra::xyFromCell(env_layer[[1]], 1:2)
    distance <- as.numeric(abs(terra::distance(distance, lonlat = T) / 1000 * factor))

    # Thinning
    invisible(utils::capture.output(
      occT <-
        spThin::thin(
          loc.data = da,
          lat.col = y,
          long.col = x,
          spec.col = id,
          thin.par = distance,
          reps = 20,
          write.files = F,
          locs.thinned.list.return = T,
          write.log.file = F
        )
    ))
    occT <-
      occT[[which(sapply(occT, function(x) {
        nrow(x)
      }) == max(sapply(occT, function(x) {
        nrow(x)
      })))[1]]]
    occPOS <- as.integer(row.names(occT))

    # Select Thinned Occurrences
    coord_filter <- data[occPOS, ]
    message("Distance threshold(km): ", round(distance, 3))
    message("Number of filtered records: ", nrow(coord_filter))
    return(dplyr::tibble(coord_filter))
  }

  if ("defined" %in% method) {

    # Thinning
    invisible(utils::capture.output(
      occT <-
        spThin::thin(
          loc.data = da,
          lat.col = y,
          long.col = x,
          spec.col = id,
          thin.par = as.numeric(method["d"]),
          reps = 20,
          write.files = F,
          locs.thinned.list.return = T,
          write.log.file = F
        )
    ))
    occT <-
      occT[[which(sapply(occT, function(x) {
        nrow(x)
      }) == max(sapply(occT, function(x) {
        nrow(x)
      })))[1]]]
    occPOS <- as.integer(row.names(occT))

    # Select Thinned Occurrences
    coord_filter <- data[occPOS, ]
    message("Distance threshold(km): ", round(as.numeric(method["d"]), 3))
    message("Number of filtered records: ", nrow(coord_filter))
    return(dplyr::tibble(coord_filter))
  }
}
