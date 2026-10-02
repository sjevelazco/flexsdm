#' Sample pseudo-absences
#'
#' @description This function provide several methods for sampling pseudo-absences, for instance
#' totally random sampling method, with the options of using different environmental and or geographical constraints.
#'
#' @param data data.frame or tibble. Database with presences
#' (or presence-absence, or presences-pseudo-absence) records, and coordinates
#' @param x character. Column name with spatial x coordinates
#' @param y character. Column name with spatial y coordinates
#' @param n integer. Number of pseudo-absences to be sampled
#' @param method character. Pseudo-absence allocation method. It is necessary to
#' provide a vector for this argument. The methods implemented are:
#' \itemize{
#' \item random: Random allocation of pseudo-absences throughout the area used for model fitting. Usage method='random'.
#' \item kmeans: Pseudo-absences are sampled randomly within environmental clusters defined by K-means
#' cluster analysis (the number of clusters will equal than number of pseudo-absences defined in 'n').
#' For this method, it is necessary to provide a raster object with environmental
#' variables. Usage method = c('kmeans', env = somevar).
#' throughout the area used for model fitting. Usage method='random'.
#' \item env_const: Pseudo-absences are environmentally constrained to regions with lower suitability values predicted by a Bioclim model. Areas with suitability values < 0.1 are selected as candidate regions to sample pseudo-absences. For this method, it is necessary to provide a raster object with environmental variables Usage method=c(method='env_const', env = somevar).
#' \item env_const_kmeans: A combination of 'env_const' and 'kmeans' methods.
#' \item geo_const: Pseudo-absences are allocated far from occurrences based on a geographical buffer. A value of the buffer width in m must be provided if raster (used in rlayer) has a longitude/latitude CRS, or in map units in other cases. Usage method=c('geo_const', width='50000').
#' \item geo_const_kmeans: A combination of 'geo_const' and 'kmeans' methods.
#' \item geoenv_const: Pseudo-absences are constrained environmentally (based on Bioclim model) and distributed geographically far from occurrences based on a geographical buffer. For this method, a raster with environmental variables stored as SpatRaster object should be provided. A value of the buffer width in m must be provided if raster (used in rlayer) has a longitude/latitude CRS, or in map units in other cases. Usage method=c('geoenv_const', width='50000', env = somevar).
#' \item geoenv_const_kmeans: Pseudo-absences are constrained using a three-level procedure; it is similar to
#' the geoenv_const with an additional step which distributes the pseudo-absences in the environmental space
#' using K-means cluster analysis. For this method, it is necessary to provide a raster object with
#' environmental variables and a value of the buffer width in m if raster (used in rlayer) has a
#' longitude/latitude CRS, or map units in other cases. Usage method=c('geoenv_const_kmeans',
#' width='50000', env = somevar).
#' }
#'
#' @param rlayer SpatRaster. A raster layer used for sampling pseudo-absence
#' A layer with the same resolution and extent that environmental variables that will be used for modeling is recommended. In the case use maskval argument, this raster layer must contain the values used to constrain sampling
#' @param maskval integer, character, or factor. Values of the raster layer used for constraining the pseudo-absence sampling
#' @param calibarea SpatVector A SpatVector which delimit the calibration area used for a given species (see \code{\link{calib_area}} function).
#' @param sp_name character. Species name for which the output will be used.
#' If this argument is used, the first output column will have the species name. Default NULL.
#'
#' @importFrom dplyr %>% select mutate
#' @importFrom stats na.exclude kmeans
#' @importFrom terra mask ext extract match as.data.frame scale
#'
#' @return A tibble object with x y coordinates of sampled pseudo-absence points
#'
#' @export
#'
#' @seealso \code{\link{sample_background}} and \code{\link{calib_area}}.
#'
#' @examples
#' \donttest{
#' require(terra)
#' require(dplyr)
#' data("spp")
#'
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
#'     method = "random",
#'     rlayer = regions,
#'     maskval = NULL,
#'     sp_name = "sp3"
#'   )
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19) # presences
#' points(ps1[-1], col = "red", cex = 0.7, pch = 19) # absences
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
#'     method = "random",
#'     rlayer = regions,
#'     maskval = samp_here
#'   )
#'
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#'
#' # Pseudo-absences sampled with K-means approach
#' set.seed(123)
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 10,
#'     method = c(method = "kmeans", env = somevar),
#'     rlayer = regions
#'   )
#'
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
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
#' # Pseudo-absences sampled with environmental and geographical constraint
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 10,
#'     method = c("geoenv_const", width = "50000", env = somevar),
#'     rlayer = regions,
#'     maskval = samp_here
#'   )
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#'
#' # Pseudo-absences sampled with environmental and geographical constraint and with k-mean clustering
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 10,
#'     method = c("geoenv_const_kmeans", width = "50000", env = somevar),
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
#'   crs = crs(somevar)
#' )
#' plot(regions, col = gray.colors(9))
#' plot(ca_ps1, add = TRUE)
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#'
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 50,
#'     method = "random",
#'     rlayer = regions,
#'     maskval = NULL,
#'     calibarea = ca_ps1
#'   )
#' plot(regions, col = gray.colors(9))
#' plot(ca_ps1, add = TRUE)
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
#'     method = "random",
#'     rlayer = regions,
#'     maskval = samp_here,
#'     calibarea = ca_ps1
#'   )
#' plot(regions, col = gray.colors(9))
#' plot(ca_ps1, add = TRUE)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' }
sample_pseudoabs <- function(data, x, y, n, method, rlayer, maskval = NULL, calibarea = NULL, sp_name = NULL) {
  . <- ID <- NULL

  if (!any(c(
    "random",
    "kmeans",
    "env_const",
    "env_const_kmeans",
    "geo_const",
    "geo_const_kmeans",
    "geoenv_const",
    "geoenv_const_kmeans"
  ) %in% method)) {
    stop(
      "argument 'method' was misused, available methods are 'random', 'kmeans', 'env_const', 'env_const_kmeans', 'geo_const', 'geo_const_kmeans', 'geoenv_const', and 'geoenv_const_kmeans'"
    )
  }

  rlayer <- rlayer[[1]]
  data <- data[, c(x, y)]

  if (!is.null(calibarea)) {
    rlayer <- rlayer %>%
      terra::crop(., calibarea) %>%
      terra::mask(., calibarea)
  }

  #### Random method ####
  if (any(method %in% "random")) {
    cell_samp <- sample_background(data = data, x = x, y = y, method = "random", n = n, rlayer = rlayer, maskval = maskval)
  }

  #### K-means method ####
  if (any(method == "kmeans")) {
    if (is.na(method["env"])) {
      stop("Provide a environmental stack/brick variables for env_const method, \ne.g. method = c('kmeans', env=somevar)")
    }

    env <- method[["env"]]
    # Test extent
    if (!all(as.vector(terra::ext(env)) %in% as.vector(terra::ext(rlayer)))) {
      message("Extents do not match, raster layers used were croped to minimum extent")
      df_ext <- data.frame(as.vector(terra::ext(env)), as.vector(terra::ext(rlayer)))

      e <- terra::ext(apply(df_ext, 1, function(x) x[which.min(abs(x))]))
      env <- crop(env, e)
      rlayer <- crop(rlayer, e)
    }

    # Restriction for a given region
    if (!is.null(maskval)) {
      if (is.factor(rlayer)) {
        lv <- terra::levels(rlayer)[[1]]
        maskval <- lv[, 1][lv[, 2] %in% as.character(maskval)]
        rlayer <- rlayer * 1
      }
      filt <- terra::match(rlayer, maskval)
      rlayer <- terra::mask(rlayer, filt)
      rm(filt)
      env <- terra::mask(env, rlayer)
    }

    # K-mean procedure
    env <- terra::scale(env)
    env_changed <- terra::as.data.frame(env, na.rm = TRUE, cells = TRUE, xy = TRUE)
    cell_samp <- kf(df = env_changed, n)
  }

  #### env_const method ####
  if (any(method == "env_const")) {
    if (is.na(method["env"])) {
      stop("Provide an environmental stack/brick variables for env_const method, \ne.g. method = c('env_const', env=somevar)")
    }

    env <- method[["env"]]
    # Test extent
    if (!all(as.vector(terra::ext(env)) %in% as.vector(terra::ext(rlayer)))) {
      message("Extents do not match, raster layers used were croped to minimum extent")
      df_ext <- data.frame(as.vector(terra::ext(env)), as.vector(terra::ext(rlayer)))

      e <- terra::ext(apply(df_ext, 1, function(x) x[which.min(abs(x))]))
      env <- crop(env, e)
      rlayer <- crop(rlayer, e)
    }

    # Restriction for a given region
    envp <- inv_bio(e = env, p = data[, c(x, y)])
    envp <- terra::mask(rlayer, envp)
    cell_samp <- sample_background(data = data, x = x, y = y, method = "random", n = n, rlayer = envp, maskval = maskval)
  }

  #### env_const_kmeans method ####
  if (any(method == "env_const_kmeans")) {
    if (is.na(method["env"])) {
      stop("Provide an environmental stack/brick variables for env_const_kmeans method, \ne.g. method = c('env_const_kmeans', env=somevar)")
    }

    env <- method[["env"]]
    # Test extent
    if (!all(as.vector(terra::ext(env)) %in% as.vector(terra::ext(rlayer)))) {
      message("Extents do not match, raster layers used were croped to minimum extent")
      df_ext <- data.frame(as.vector(terra::ext(env)), as.vector(terra::ext(rlayer)))

      e <- terra::ext(apply(df_ext, 1, function(x) x[which.min(abs(x))]))
      env <- crop(env, e)
      rlayer <- crop(rlayer, e)
    }

    # Restriction for a given region
    envp <- inv_bio(e = env, p = data[, c(x, y)])
    envp <- terra::mask(rlayer, envp)

    if (!is.null(maskval)) {
      if (is.factor(rlayer)) {
        lv <- terra::levels(rlayer)[[1]]
        maskval <- lv[, 1][lv[, 2] %in% as.character(maskval)]
        rlayer <- rlayer * 1
      }
      filt <- terra::match(rlayer, maskval)
      envp <- terra::mask(envp, filt)
      rm(filt)
    }

    env_changed <- terra::mask(env, envp)
    env_changed <- terra::scale(env_changed)
    env_changed <- terra::as.data.frame(env_changed, na.rm = TRUE, cells = TRUE, xy = TRUE)
    cell_samp <- kf(df = env_changed, n)
  }

  #### geo_const method ####
  if (any(method == "geo_const")) {
    if (!"width" %in% names(method)) {
      stop("Provide a width value for 'geo_const' method, \ne.g. method=c('geo_const', width='50000')")
    }

    # Restriction for a given region
    envp <- inv_geo(e = rlayer, p = data[, c(x, y)], d = as.numeric(method["width"]))
    cell_samp <- sample_background(data = data, x = x, y = y, method = "random", n = n, rlayer = envp, maskval = maskval)
  }

  #### geo_const_kmeans method ####
  if (any(method == "geo_const_kmeans")) {
    if (!all(c("env", "width") %in% names(method))) {
      stop("Provide a width value and environmental stack/brick variables for 'geo_const_kmeans' method, \ne.g. method=c('geo_const_kmeans', width='50000', env=somevar)")
    }

    env <- method[["env"]]
    # Test extent
    if (!all(as.vector(terra::ext(env)) %in% as.vector(terra::ext(rlayer)))) {
      message("Extents do not match, raster layers used were croped to minimum extent")
      df_ext <- data.frame(as.vector(terra::ext(env)), as.vector(terra::ext(rlayer)))

      e <- terra::ext(apply(df_ext, 1, function(x) x[which.min(abs(x))]))
      env <- crop(env, e)
      rlayer <- crop(rlayer, e)
    }

    # Restriction for a given region
    envp <- inv_geo(e = rlayer, p = data[, c(x, y)], d = as.numeric(method["width"]))

    if (!is.null(maskval)) {
      if (is.factor(rlayer)) {
        lv <- terra::levels(rlayer)[[1]]
        maskval <- lv[, 1][lv[, 2] %in% as.character(maskval)]
        rlayer <- rlayer * 1
      }
      filt <- terra::match(rlayer, maskval)
      envp <- terra::mask(envp, filt)
      rm(filt)
    }

    env_changed <- terra::mask(env, envp)
    env_changed <- terra::scale(env_changed)
    env_changed <- terra::as.data.frame(env_changed, na.rm = TRUE, cells = TRUE, xy = TRUE)
    cell_samp <- kf(df = env_changed, n)
  }

  #### geoenv_const method ####
  if (any(method == "geoenv_const")) {
    if (!all(c("env", "width") %in% names(method))) {
      stop("Provide a width value and environmental stack/brick variables for 'geoenv_const' method, \ne.g. method=c('geoenv_const', width='50000', env=somevar)")
    }

    env <- method[["env"]]

    # Test extent
    if (!all(as.vector(terra::ext(env)) %in% as.vector(terra::ext(rlayer)))) {
      message("Extents do not match, raster layers used were croped to minimum extent")
      df_ext <- data.frame(as.vector(terra::ext(env)), as.vector(terra::ext(rlayer)))

      e <- terra::ext(apply(df_ext, 1, function(x) x[which.min(abs(x))]))
      env <- crop(env, e)
      rlayer <- crop(rlayer, e)
    }

    # Restriction for a given region
    envp <- inv_geo(e = rlayer, p = data[, c(x, y)], d = as.numeric(method["width"]))
    envp2 <- inv_bio(e = env, p = data[, c(x, y)])

    envp <- (envp2 + envp)
    rm(envp2)
    envp <- terra::mask(rlayer, envp)
    cell_samp <- sample_background(data = data, x = x, y = y, method = "random", n = n, rlayer = envp, maskval = maskval)
  }

  #### geoenv_const_kmeans ####
  if (any(method == "geoenv_const_kmeans")) {
    if (!all(c("env", "width") %in% names(method))) {
      stop("Provide a width value and environmental stack/brick variables for 'geoenv_const_kmeans' method, \ne.g. method=c('geoenv_const_kmeans', width='50000', env=somevar)")
    }

    env <- method[["env"]]

    # Test extent
    if (!all(as.vector(terra::ext(env)) %in% as.vector(terra::ext(rlayer)))) {
      message("Extents do not match, raster layers used were croped to minimum extent")
      df_ext <- data.frame(as.vector(terra::ext(env)), as.vector(terra::ext(rlayer)))

      e <- terra::ext(apply(df_ext, 1, function(x) x[which.min(abs(x))]))
      env <- crop(env, e)
      rlayer <- crop(rlayer, e)
    }

    # Restriction for a given region
    envp <- inv_geo(e = rlayer, p = data[, c(x, y)], d = as.numeric(method["width"]))
    envp2 <- inv_bio(e = env, p = data[, c(x, y)])
    envp <- (envp2 + envp)
    envp <- terra::mask(rlayer, envp)
    rm(envp2)

    if (!is.null(maskval)) {
      if (is.factor(rlayer)) {
        lv <- terra::levels(rlayer)[[1]]
        maskval <- lv[, 1][lv[, 2] %in% as.character(maskval)]
        rlayer <- rlayer * 1
      }
      filt <- terra::match(rlayer, maskval)
      rlayer <- terra::mask(rlayer, filt)
      rm(filt)
    }

    envp <- terra::mask(rlayer, envp)

    # K-mean procedure
    env_changed <- terra::mask(env, envp)
    env_changed <- terra::scale(env_changed)
    env_changed <- terra::as.data.frame(env_changed, na.rm = TRUE, cells = TRUE, xy = TRUE)
    env_changed <- stats::na.exclude(env_changed)
    cell_samp <- kf(df = env_changed, n)
  }

  colnames(cell_samp) <- c(x, y, "pr_ab")
  if (!is.null(sp_name)) {
    cell_samp <- tibble(sp = sp_name, cell_samp)
  }
  return(cell_samp)
}
