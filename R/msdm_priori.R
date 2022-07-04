#' Create spatial predictor variables to reduce overprediction of species distribution models
#'
#' @description This function creates geographical predictor variables that, together with
#' environmental variables, can be used to construct constrained species distribution models.
#'
#' @param data tibble or data.frame. A database with geographical coordinates of species presences.
#' @param x character. Column name with spatial x coordinates.
#' @param y character. Column name with spatial y coordinates.
#' @param method character. A character string indicating which MSDM method will be used. The
#' following methods are available: 'xy', 'min', 'cml', and 'ker'. Usage method = 'cml'
#' @param env_layer A raster layer used to construct species distribution models. This object will
#' be used to create constraining variables with the same resolution, extent, and pattern of empty
#' cells as the environmental variables. It is advisable to use a raster of an environmental layer
#' that will be used to create the species distribution models to avoid mismatch (e.g. resolution,
#' extent, cells with NA) between environmental and constraining variables.
#'
#' @return This function returns a SpatRaster object. Such raster/s have to be used together with
#' environmental variables to construct species distribution models.
#' The 'xy' approach creates a single pair of raster layers that can be used for all species that
#' share the same study region. Otherwise,
#' 'cml', 'min', and 'ker' create a species-specific raster layer.
#'
#' @details
#' This function creates geographical predictor variables that, together with environmental
#' variables, can be used to construct constrained species distribution models.
#' It is recommended to use these approaches to create models that will only be projected for current
#' conditions and not for different time periods (past or future).
#'
#' Four methods are implemented:
#'
#' xy (Latlong method). This method assumes that spatial structure can partially explain species
#' distribution (Bahn & McGill, 2007). Therefore, two raster layers will be created, containing the
#' latitude and longitude of pixels, respectively. These raster layers should be included as
#' covariates with the environmental layers to construct species distribution models. This method
#' does not interact with species occurrence and is generic for a given study region; for this
#' reason, it is possible to use this method for all species set that share the same study region.
#'
#' min (Nearest neighbor distance method). Compiled and adapted from Allouche et al. (2008), this
#' method calculates for each cell the Euclidean geographic distance to the nearest presence point.
#'
#' cml (Cumulative distance method). Compiled and adapted from Allouche et al. (2008), this method
#' assumes that pixels closer to presences are likely included in species distributions. Therefore,
#' a raster layer will be created containing the sum of Euclidean geographic distances from each
#' pixel to all occurrences of a species. Obtained values are normalized to vary from zero to one.
#' This raster layer should be included with the environmental layers to construct species
#' distribution models.
#'
#' ker (Kernel method). Compiled and adapted from Allouche et al. (2008), this method, like cml,
#' assumes that pixels located in areas with a higher density of occurrences are likely included in
#' the actual species distribution. Thus, a raster layer will be created containing the Gaussian
#' values based on the density of occurrences of a species. The standard deviation of the Gaussian
#' distribution was the maximum value in a vector of minimum distances between pairs of occurrences
#' of a species. Gaussian values are normalized to vary from zero to one. This raster layer should
#' be included  with the environmental layers to construct species distribution models.
#'
#'
#' See Mendes et al. (2020) for further methodological and performance details.
#'
#' If used one these constraining method cite Mendes et al 2020.
#'
#' @references
#' \itemize{
#' \item Mendes, P.; Velazco S.J.E.; Andrade, A.F.A.; De Marco, P. (2020) Dealing with
#' overprediction in species distribution models: how adding distance constraints can improve model
#' accuracy, Ecological Modelling, in press. https://doi.org/10.1016/j.ecolmodel.2020.109180
#' \item Allouche, O.; Steinitz, O.; Rotem, D.; Rosenfeld, A.; Kadmon, R. (2008). Incorporating
#' distance constraints into species distribution models. Journal of Applied Ecology, 45(2),
#' 599-609. doi:10.1111/j.1365-2664.2007.01445.x
#' \item Bahn, V.; McGill, B. J. (2007). Can niche-based distribution models outperform spatial
#' interpolation? Global Ecology and Biogeography, 16(6), 733-742.
#' doi:10.1111/j.1466-8238.2007.00331.x
#' }
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' require(terra)
#'
#' data("spp")
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' # Select the presences of a species
#' occ <- spp %>%
#'   dplyr::filter(species == "sp3", pr_ab == 1)
#'
#' # Select a raster layer to be used as a basic raster
#' a_variable <- somevar[[1]]
#' plot(a_variable)
#' points(occ %>% dplyr::select(x, y))
#'
#' ### xy method
#' m_xy <- msdm_priori(
#'   data = occ,
#'   x = "x",
#'   y = "y",
#'   method = "xy",
#'   env_layer = a_variable
#' )
#'
#' plot(m_xy)
#'
#' ### min method
#' m_min <- msdm_priori(
#'   data = occ,
#'   x = "x",
#'   y = "y",
#'   method = "min",
#'   env_layer = a_variable
#' )
#'
#' plot(m_min)
#'
#' ### cml method
#' m_cml <- msdm_priori(
#'   data = occ,
#'   x = "x",
#'   y = "y",
#'   method = "cml",
#'   env_layer = a_variable
#' )
#'
#' plot(m_cml)
#'
#' ### ker method
#' m_ker <- msdm_priori(
#'   data = occ,
#'   x = "x",
#'   y = "y",
#'   method = "ker",
#'   env_layer = a_variable
#' )
#'
#' plot(m_ker)
#' }
#'
#' @seealso \code{\link{msdm_posteriori}}
#'
#' @importFrom dplyr tibble select all_of
#' @importFrom stats na.omit
#' @importFrom methods is
#' @importFrom terra rast as.data.frame cellFromXY xyFromCell
#'
#' @export
msdm_priori <- function(data,
                        x,
                        y,
                        method = c("xy", "min", "cml", "ker"),
                        env_layer) {
  if (is.null(env_layer)) {
    stop("Complete env_layer argument")
  }

  ### Transform inputs

  # raster
  if (!methods::is(env_layer, "SpatRaster")) {
    env_layer <- terra::rast(env_layer)
  }
  env_layer <- env_layer[[1]]

  # database
  data <- dplyr::tibble(data)
  records <-
    data %>% dplyr::select(dplyr::all_of(x), dplyr::all_of(y))
  names(records) <- c("x", "y")


  # Prepare data for calculation
  if (any(method %in% c("min", "cml", "ker"))) {
    spi <- terra::as.data.frame(env_layer,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
      dplyr::select(x, y)

    r <- terra::cellFromXY(env_layer, xy = as.matrix(records))
    r <- terra::xyFromCell(object = env_layer, cell = r) %>%
      data.frame() %>%
      stats::na.omit(r)
  }


  #### Method 1: xy ####
  if (method == "xy") {
    # Extract coordinates of raster layer
    df <-
      terra::as.data.frame(env_layer,
        xy = TRUE,
        na.rm = TRUE,
        cells = TRUE
      )[1:3]
    lon <- lat <- env_layer
    lon[df$cell] <- df$x
    lat[df$cell] <- df$y

    # Create LatLong Stack
    result <- terra::rast(list(lon, lat))
    names(result) <- c("msdm_lon", "msdm_lat")
    rm(lat, lon, df)

    return(result)
  }

  # Method 2- Minimum distance ----
  if (method == "min") {
    distr <- euc_dist(x = spi, y = r)
    distr <- apply(distr, 1, min)
    distr <- (distr - min(distr)) / (max(distr) - min(distr))
    result <- env_layer
    result[!is.na(result)] <- distr
    names(result) <- "msdm_min"
    rm(distr, spi, r)
    return(result)
  }

  # Method 3: Cummulative distance----
  if (method == "cml") {
    distr <- euc_dist(spi, r)
    distr <- distr + 1
    distr <- 1 / (1 / distr^2)
    distr <- apply(distr, 1, sum)
    distr <- (distr - min(distr)) / (max(distr) - min(distr))

    result <- env_layer
    result[!is.na(result)] <- distr
    names(result) <- "msdm_cml"
    rm(distr, spi, r)
    return(result)
  }

  # Method 4: Gaussian Kernel----
  if (method == "ker") {
    distp <- euc_dist(r, r)
    distp1 <- matrix(0, nrow(euc_dist(r, r)), 1)

    for (c in 1:nrow(distp)) {
      vec <- distp[c, ]
      distp1[c] <- min(vec[vec != min(vec)])
    }
    rm(vec)
    sd_graus <- max(distp1)
    distr <- euc_dist(spi, r)
    # distr2 <- distr
    distr <-
      (1 / sqrt(2 * pi * sd_graus) * exp(-1 * (distr / (2 * sd_graus^
        2))))
    distr <- apply(distr, 1, sum)
    distr <- (distr - min(distr)) / (max(distr) - min(distr))
    result <- env_layer
    result[!is.na(result)] <- distr
    rm(distr, spi, r)
    names(result) <- "msdm_ker"
    return(result)
  }
}
