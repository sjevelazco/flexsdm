#' Peseudo-absence allocation method
#'
#' @param data data.frame or tibble. Database with presences
#' (or presence-absence, o presences-pseudo-absence) records, and coordinates
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param n integer. Number of pseudo-absence to be sampled
#' @param method character. Pseudo-absence allocation method. It is necessary to
#' provide a vector for this argument. The next methods are implemented:
#' \itemize{
#' \item rnd: Random allocation of pseudo-absences throughout the area used for model fitting. Usage method='rnd'.
#' \item env_const: Pseudo-absences are environmentally constrained to regions with lower suitability values predicted by a Bioclim model. For this method, it is necessary to provide a raster stack or brick object with environmental variables Usage method=c(method='env_const', env = somevar).
#' \item geo_const: Pseudo-absences are allocated far from occurrences based on a geographical buffer. For this method, it is necessary to provide a value of the buffer width in m if raster (used in rlayer) has a longitude/latitude CRS, or map units in other cases. Usage method=c('geo_const', width='50000').
#' \item geo_env_const: Pseudo-absences are constrained environmentally (based on Bioclim model) and distributed geographically far from occurrences based on a geographical buffer. For this method, it is necessary to provide a raster with environmental variables stored as SpatRaster object. Also it is necessary provide a value of the buffer width in m if raster (used in rlayer) has a longitude/latitude CRS, or map units in other cases. Usage method=c('geo_env_const', width='50000', env = somevar)
#' \item geo_env_km_const: Pseudo-absences are constrained on a three-level procedure; it is similar to the geo_env_const with an additional step which distributes the pseudo-absences in the environmental space using k-means cluster analysis. For this method, it is necessary to provide a raster stack or brick object with environmental variables and a value of the buffer width in m if raster (used in rlayer) has a longitude/latitude CRS, or map units in other cases. Usage method=c('geo_env_km_const', width='50000', env = somevar)
#' }
#'
#' @param rlayer SpatRaster. A raster layer used for sampling pseudo-absence
#' It is recommended to use a layer with the same resolution and extent that environmental variables that will be used for modeling. In the case use maskval argument, this raster layer must contain the values to sampling constraint
#' @param maskval integer or numeric. Values of the raster layer used for constraining the pseudo-absence sampling
#' @param calibarea SpatVector A SpatVector which delimit the calibration area used for a given species (see \code{\link{calib_area}} function).
#'
#' @importFrom dplyr %>% select mutate
#' @importFrom stats na.exclude kmeans
#' @importFrom terra mask ext extract match as.data.frame
#'
#' @return
#' @export
#'
#'
#' @seealso \code{\link{sample_background}} and \code{\link{calib_area}}.
#'
#' @examples
#' \dontrun{
#' data("spp")
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' regions <- system.file("external/regions.tif", package = "flexsdm")
#' regions <- terra::rast(regions)
#'
#' plot(regions)
#'
#'
#' single_spp <-
#'   spp %>%
#'   dplyr::filter(species == "sp3") %>%
#'   dplyr::filter(pr_ab == 1) %>%
#'   dplyr::select(-pr_ab)
#'
#'
#' # Pseudo-absences randomly sampled throughout study area
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 10,
#'     method = "rnd",
#'     rlayer = regions,
#'     maskval = NULL
#'   )
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19) # presences
#' points(ps1, col = "red", cex = 0.7, pch = 19) # absences
#'
#'
#' # Pseudo-absences randomly sampled within a regions where a species occurs
#' ## Regions where this species occurrs
#' samp_here <- terra::extract(regions, single_spp[2:3])[, 2] %>%
#'   unique() %>%
#'   na.exclude()
#'
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 10,
#'     method = "rnd",
#'     rlayer = regions,
#'     maskval = samp_here
#'   )
#'
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#'
#'
#' # Pseudo-absences sampled with geographical constraint
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 10,
#'     method = c("geo_const", width = "30000"),
#'     rlayer = regions,
#'     maskval = samp_here
#'   )
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#'
#' # Pseudo-absences sampled with environmental constraint
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 10,
#'     method = c("env_const", env = somevar),
#'     rlayer = regions,
#'     maskval = samp_here
#'   )
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#'
#' # Pseudo-absences sampled with environmental and geographical contraint
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 10,
#'     method = c("geo_env_const", width = "50000", env = somevar),
#'     rlayer = regions,
#'     maskval = samp_here
#'   )
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#'
#' # Pseudo-absences sampled with environmental and geographical contraint and with k-mean
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 10,
#'     method = c("geo_env_km_const", width = "50000", env = somevar),
#'     rlayer = regions,
#'     maskval = samp_here
#'   )
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#'
#' # Sampling pseudo-absence using a calibration area
#' ca_ps1 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = c("buffer", width = 50000),
#' )
#' plot(regions, col = gray.colors(9))
#' plot(ca_ps1, add = T)
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#'
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 50,
#'     method = "rnd",
#'     rlayer = regions,
#'     maskval = NULL,
#'     calibarea = ca_ps1
#'   )
#' plot(regions, col = gray.colors(9))
#' plot(ca_ps1, add = T)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#'
#'
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 50,
#'     method = "rnd",
#'     rlayer = regions,
#'     maskval = samp_here,
#'     calibarea = ca_ps1
#'   )
#' plot(regions, col = gray.colors(9))
#' plot(ca_ps1, add = T)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' }
sample_pseudoabs <- function(data, x, y, n, method, rlayer, maskval = NULL, calibarea = NULL) {
  ID <- NULL

  if (!any(c(
    "rnd",
    "env_const",
    "geo_const",
    "geo_env_const",
    "geo_env_km_const"
  ) %in% method)) {
    stop(
      "argument 'method' was misused, available methods rnd, env_const, geo_const, geo_env_const, and geo_env_km_const"
    )
  }

  rlayer <- rlayer[[1]]
  data <- data[, c(x, y)]

  if (!is.null(calibarea)) {
    rlayer <- terra::mask(rlayer, calibarea)
  }

  # Random method
  if (any(method %in% "rnd")) {
    cell_samp <- sample_background(n = n, rlayer = rlayer, maskval = maskval)
  }

  # env_const method
  if (any(method == "env_const")) {
    if (is.null(method["env"])) {
      stop("Provide a environmental stack/brick variables for env_const method, \ne.g. method = c('env_const', env=somevar)")
    }

    env <- method[["env"]]

    # Restriction for a given region
    envp <- inv_bio(e = env, p = data[, c(x, y)])
    if (terra::ext(env) != terra::ext(rlayer)) {
      rlayer2 <- rlayer
      rlayer2[] <- terra::extract(envp, coordinates(rlayer2), method = "simple", xy = TRUE) %>%
        dplyr::select(-c(ID, x, y))
      envp <- rlayer2
    }

    envp <- terra::mask(rlayer, envp)
    cell_samp <- sample_background(n = n, rlayer = envp, maskval = maskval)
  }

  # geo_const method
  if (any(method == "geo_const")) {
    if (!"width" %in% names(method)) {
      stop("Provide a width value for 'geo_const' method, \ne.g. method=c('geo_const', width='50000')")
    }

    # Restriction for a given region
    envp <- inv_geo(e = rlayer, p = data[, c(x, y)], d = as.numeric(method["width"]))

    cell_samp <- sample_background(n = n, rlayer = envp, maskval = maskval)
  }

  # geo_env_const method
  if (any(method == "geo_env_const")) {
    if (!all(c("env", "width") %in% names(method))) {
      stop("Provide a width value and environmental stack/brick variables for 'geo_env_const' method, \ne.g. method=c('geo_env_const', width='50000', env=somevar)")
    }

    env <- method[["env"]]

    # Restriction for a given region
    envp <- inv_geo(e = rlayer, p = data[, c(x, y)], d = as.numeric(method["width"]))
    envp2 <- inv_bio(e = env, p = data[, c(x, y)])
    if (terra::ext(env) != terra::ext(rlayer)) {
      rlayer2 <- rlayer
      rlayer2[] <- terra::extract(envp2, coordinates(rlayer2), method = "simple", xy = TRUE) %>%
        dplyr::select(-c(ID, x, y))
      envp2 <- rlayer2
    }
    envp <- (envp2 + envp)
    rm(envp2)
    envp <- terra::mask(rlayer, envp)
    cell_samp <- sample_background(n = n, rlayer = envp, maskval = maskval)
  }

  # geo_env_km_const
  if (any(method == "geo_env_km_const")) {
    if (!all(c("env", "width") %in% names(method))) {
      stop("Provide a width value and environmental stack/brick variables for 'geo_env_km_const' method, \ne.g. method=c('geo_env_const', width='50000', env=somevar)")
    }

    env <- method[["env"]]

    # Restriction for a given region
    envp <- inv_geo(e = rlayer, p = data[, c(x, y)], d = as.numeric(method["width"]))
    envp2 <- inv_bio(e = env, p = data[, c(x, y)])
    if (terra::ext(env) != terra::ext(rlayer)) {
      rlayer2 <- rlayer
      rlayer2[] <- terra::extract(envp2, coordinates(rlayer2), method = "simple", xy = TRUE) %>%
        dplyr::select(-c(ID, x, y))
      envp2 <- rlayer2
    }
    envp <- (envp2 + envp)
    envp <- terra::mask(rlayer, envp)
    rm(envp2)

    if (!is.null(maskval)) {
      if (is.factor(maskval)) {
        maskval <-
          which(levels(maskval)[-1] %in% as.character(maskval))
        rlayer <- rlayer * 1
      }
      filt <- terra::match(rlayer, maskval)
      rlayer <- terra::mask(rlayer, filt)
      rm(filt)
    }

    envp <- terra::mask(rlayer, envp)

    # K-mean procedure
    if (terra::ext(env) != terra::ext(envp)) {
      env <- crop(env, envp)
    }

    env_changed <- terra::mask(env, envp)
    env_changed <- terra::as.data.frame(env_changed, xy = TRUE)
    env_changed <- stats::na.exclude(env_changed)

    suppressWarnings(km <- stats::kmeans(env_changed, centers = n))
    cell_samp <- km$centers[, 1:2] %>% data.frame()
    val <- terra::extract(envp, cell_samp, method = "simple", xy = TRUE) %>%
      dplyr::select(-c(ID, x, y))
    cell_samp <-
      cell_samp %>% dplyr::mutate(val = val[, 1])
    cell_samp <- cell_samp[!is.na(cell_samp$val), -3]
  }
  colnames(cell_samp) <- c(x, y)
  return(cell_samp)
}
