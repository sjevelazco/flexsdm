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
#' based on a geographical distance criteria. For the three method, it is possible to use several values.
#' If several values are provided, the function will return a list with the results. The following methods are available:
#' \itemize{
#'   \item moran: records are filtered based on the smallest distance to the Moran's I value provided. If no Moran's I
#'   values is provided it will use 0.1. Usage method: method = c('moran') or  method = c('moran', val = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35)).
#'   \item cellsize: records are filtered based on the resolution of the environmental variables which
#'   can be aggregated to coarser resolution defined by the factor.
#'   Usage method: method = c('cellsize', factor = '2') or method = c('cellsize', factor = c(1, 4, 8)).
#'   \item defined: records are filtered based on a  distance value (d) provided in km.
#'   Usage method: method = c('defined', d = c(20, 40, 60, 80)).
#' }
#' @param prj character. Projection string (PROJ4) for occurrences. Not necessary if
#' the projection used is WGS84 ("+proj=longlat +datum=WGS84"). Default "+proj=longlat +datum=WGS84"
#' @param reps integer. Number of times to repeat the thinning process. Default 20
#'
#' @return
#' If one value is used to filter occurrence funtion will return a tibble object with filtered data. If several
#' values are used to filter occurrences, the function will return a list of tibbles with filtered data.
#'
#' @details In this function three alternatives are implemented to determine the
#' distance threshold between pair of points:
#' ' \itemize{
#'   \item "moran" determines the minimum nearest-neighbor distance that approximate to the spatial
#' autocorrelation in occurrence data, following a Moran's I. To do so, a Principal Component
#' Analysis with the environmental variables is performed and then the first Principal Component is used
#' to calculate the semivariograms. Because of this, this method only allow the use of continuous variables.
#' Sometimes, this method can (too) greatly reduce the number of presences.
#'   \item "cellsize" filters occurrences based on the predictors' resolution. This method will calculate
#' the distance between the first two cells of the environmental variable and use this distance
#' as minimum nearest-neighbor distance to filter occurrences.
#' The resolution of the raster is aggregated based on the values used in "factor". Thus, the distance
#' used for filtering can be adjusted to represent a larger grid size.
#'  \item "determined" this method uses any minimum nearest-neighbor distance specified in km.
#'   }
#'
#' The "thin" function from spThin package is used to filter data
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
#' @importFrom dplyr tibble as_tibble left_join
#' @importFrom spThin thin
#' @importFrom stats complete.cases prcomp
#' @importFrom terra vect project geom extract as.data.frame predict distance xyFromCell
#' @importFrom utils capture.output
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#' require(ggplot2)
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
#' ##%######################################################%##
#' ####                  Cellsize method                   ####
#' ##%######################################################%##
#'
#' # Using cellsize method
#' filtered_occ <- occfilt_geo(
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
#' points(filtered_occ %>% select(x, y), pch = 19, col = "yellow") # filtered data
#'
#'
#' # Using cellsize method with several values
#' filtered_occ <- occfilt_geo(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   env_layer = somevar,
#'   method = c("cellsize", factor = c(1, 4, 8, 12, 16, 20)),
#'   prj = crs(somevar)
#' )
#'
#' filtered_occ # Note that several values are provided for any filtering method
#' # fuction will return a list of tibbles with the results.
#' # So user must select the desired filtered dataset
#'
#' # Let's explore the results
#' bind_rows(filtered_occ, .id="cellSize") %>%
#'   dplyr::mutate(cellSize = as.numeric(cellSize)) %>%
#'   ggplot(aes(x, y)) +
#'   geom_point() +
#'   facet_wrap(~cellSize)
#'
#'
#' ##%######################################################%##
#' ####                   Defined method                   ####
#' ##%######################################################%##
#' # Using defined method
#' filtered_occ <- occfilt_geo(
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
#' points(filtered_occ %>% select(x, y), pch = 19, col = "yellow") # filtered data
#'
#' # Using defined method with several values
#' filtered_occ <- occfilt_geo(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   env_layer = somevar,
#'   method = c("defined", factor = c(5, 10, 15, 30, 35, 40)),
#'   prj = crs(somevar)
#' )
#'
#' bind_rows(filtered_occ, .id="cellSize") %>%
#'   dplyr::mutate(cellSize = as.numeric(cellSize)) %>%
#'   ggplot(aes(x, y)) +
#'   geom_point() +
#'   facet_wrap(~cellSize)
#'
#'
#' ##%######################################################%##
#' ####                  Moran's I method                  ####
#' ##%######################################################%##
#'
#' # Using Moran's I method
#' filtered_occ <- occfilt_geo(
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
#' points(filtered_occ %>% select(x, y), pch = 19, col = "yellow") # filtered data
#'
#' # Using Moran's I method with several values
#' filtered_occ <- occfilt_geo(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   env_layer = somevar,
#'   method = c("moran", c(0.05, 0.15, 0.2, 0.5, 0.7)),
#'   prj = crs(somevar)
#' )
#'
#' bind_rows(filtered_occ, .id="moran") %>%
#'   dplyr::mutate(moran = as.numeric(moran)) %>%
#'   ggplot(aes(x, y)) +
#'   geom_point() +
#'   facet_wrap(~moran)
#' }
#'
#' @seealso \code{\link{occfilt_env}}
#'
occfilt_geo <- function(data, x, y, env_layer, method, prj = "+proj=longlat +datum=WGS84", reps = 20) {

  da <- data0 <- data
  rm(data)
  da <- da[c(x, y)]


  if (prj != "+proj=longlat +datum=WGS84") {
    da <- terra::vect(da, geom = c(x, y), prj)
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
  env <- terra::extract(env_layer, da[c(x, y)])
  env$ID <- NULL

  filt <- stats::complete.cases(env)
  if (sum(!filt) > 0) {
    message(sum(!filt), " records were removed because they have NAs for some variables")
    da <- da[filt, ]
    data0 <- data0[filt, ]
    coord <- coord[filt, ]
  }
  rm(filt)
  rm(env)

  message("Number of unfiltered records: ", nrow(data), "\n")


  occfilt_geo_0 <- function(data, da, x, y, env_layer, method, reps) {
    #### Morans'I ####
    if ("moran" %in% method) {
      # Extract occurrence PC1 values
      vals <- terra::extract(env_layer, coord)$PC1

      # Trasform coord to vect
      coord <- vect(coord, geom = c(x, y), crs = "+proj=longlat +datum=WGS84")

      # Calculate distances & Define breaks
      dists <- terra::distance(coord, coord) / 1000 # dist. in km
      # br <- max(dists) / 2000
      dists[dists == 0] <- NA
      br <- min(dists, na.rm = TRUE) # min distance of 1 km
      br <- seq(from = br, to = max(dists, na.rm = TRUE), length.out = 1000)

      # Calculate Moran for breaks
      mor <- rep(NA, length(br))
      for (i in 1:length(br)) {
        dists2 <- dists
        dists2[dists2 > br[i]] <- 0
        dists2 <- as.matrix(dists2)
        diag(dists2) <- 0
        mor[i] <- abs(morani(x = vals, weight = dists2, scaled = TRUE))
      }

      # Select threshold
        if (any(mor <= method[2])) {
          mor_thr <- as.numeric(method[-1])
          pos <- which(mor <= mor_thr & mor > mor_thr - 0.005)
          # Filter distance and imoran values
          br <- br[pos]
          mor <- mor[pos]
          # Select threshold
          d <- min(br)
          mor <- mor[which.min(br)]
          message("Moran's I threshold closest to the supplied value: ", round(mor[br==d], 3))
        } else {
          message("No Moran's I <= than the supplied value were found")
          message("Moran's I threshold used: ", round(min(mor), 2))
          pos <- which(mor == min(mor))
          d <- br[pos]
        }

      options(warn = 1)

      # Thinning
      da$.spp <- "sp"
      set.seed(1)
      invisible(utils::capture.output(
        occT <-
          spThin::thin(
            loc.data = da,
            lat.col = y,
            long.col = x,
            spec.col = ".spp",
            thin.par = d,
            reps = reps,
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
      colnames(occT) <- c(x, y)
      occT <- as.integer(row.names(occT))

      # Select Thinned Occurrences
      coord_filter <- data0[occT, ]

      # Results
      message("Distance threshold (km) : ", round(d, 3))
      message("Number of filtered records: ", nrow(coord_filter), "\n")
    }

    #### cellsize ####
    if ("cellsize" %in% method) {
      # Haversine Transformation & Distance Threshold
      factor <- as.numeric(method["factor"])
      distance <-
        terra::xyFromCell(env_layer[[1]], 1:2)
      distance <- as.numeric(abs(terra::distance(distance, distance, lonlat = TRUE)[2] / 1000 * factor))

      # Thinning
      da$.spp <- "sp"
      set.seed(1)
      invisible(utils::capture.output(
        occT <-
          spThin::thin(
            loc.data = da,
            lat.col = y,
            long.col = x,
            spec.col = ".spp",
            thin.par = distance,
            reps = reps,
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
      message("Factor: ", paste0("x", factor))
      message("Distance threshold (km): ", round(distance, 3))
      message("Number of filtered records: ", nrow(coord_filter), "\n")
    }

    #### defined ####
    if ("defined" %in% method) {
      # Thinning
      da$.spp <- "sp"
      set.seed(1)
      invisible(utils::capture.output(
        occT <-
          spThin::thin(
            loc.data = da,
            lat.col = y,
            long.col = x,
            spec.col = ".spp",
            thin.par = as.numeric(method["d"]),
            reps = reps,
            write.files = FALSE,
            locs.thinned.list.return = TRUE,
            write.log.file = FALSE
          )
      ))
      occT <-
        occT[[which.max(sapply(occT, function(x) {
          nrow(x)
        }))[1]]]
      colnames(occT) <- c(x, y)

      # Select Thinned Occurrences
      occT <- as.integer(row.names(occT))

      # Select Thinned Occurrences
      coord_filter <- data0[occT, ] %>%
        dplyr::as_tibble()
      message("Distance threshold (km): ", round(as.numeric(method["d"]), 3))
      message("Number of filtered records: ", nrow(coord_filter), "\n")
    }
    return(dplyr::tibble(coord_filter))
  }


  ##%######################################################%##
  ####    Approach for testing several filter values      ####
  ##%######################################################%##

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
    env_layer <- env_layer[[1]]
    rm(p)

    # Set 0.1 if value was not supplied
    if (length(method[-1]) == 0) {
      method <- c("moran", 0.1)
    }

    # Run filter
    val <- as.numeric(method[-1])
    result <- list()
    for(ii in 1:length(val)){
      result[[ii]] <- occfilt_geo_0(
        data = data0, # data without NA
        da = da,
        x = x,
        y = y,
        env_layer = env_layer,
        method = c("moran", factor = val[ii]),
        reps = reps
      )

    }
    names(result) <- val
  }


  if ("cellsize" %in% method) {
    val <- as.numeric(method[-1])
    result <- list()
    for(ii in 1:length(val)){
      result[[ii]] <- occfilt_geo_0(
        data = data0, # data without NA
        da = da,
        x = x,
        y = y,
        env_layer = env_layer,
        method = c("cellsize", factor = val[ii]),
        reps = reps
      )

    }
    names(result) <- val
  }


  if ("defined" %in% method) {
    val <- as.numeric(method[-1])
    result <- list()
    for(ii in 1:length(val)){
      result[[ii]] <- occfilt_geo_0(
        data = data0,# data without NA
        da = da,
        x = x,
        y = y,
        env_layer = env_layer,
        method = c("defined", d = val[ii]),
        reps = reps
      )
    }
    names(result) <- val
  }

  if(length(val) == 1) {
    result <- result[[1]]
  }

return(result)

}

