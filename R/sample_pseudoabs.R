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
#' \item RND: Random allocation of pseudo-absences throughout the area used for model fitting. Usage method='RND'.
#' \item ENV_CONST: Pseudo-absences are environmentally constrained to regions with lower suitability values predicted by a Bioclim model. For this method, it is necessary to provide a raster stack or brick object with environmental variables Usage method=c(method='ENV_CONST', env = somevar).
#' \item GEO_CONST: Pseudo-absences are allocated far from occurrences based on a geographical buffer. For this method, it is necessary to provide a value of the buffer width in m if raster (used in rlayer) has a longitude/latitude CRS, or map units in other cases. Usage method=c('GEO_CONST', width='50000').
#' \item GEO_ENV_CONST: Pseudo-absences are constrained environmentally (based on Bioclim model) and distributed geographically far from occurrences based on a geographical buffer. For this method, it is necessary to provide a raster stack or brick object with environmental variables and a value of the buffer width in m if raster (used in rlayer) has a longitude/latitude CRS, or map units in other cases. Usage method=c('GEO_ENV_CONST', width='50000', env = somevar)
#' \item GEO_ENV_KM_CONST: Pseudo-absences are constrained on a three-level procedure; it is similar to the GEO_ENV_CONST with an additional step which distributes the pseudo-absences in the environmental space using k-means cluster analysis. For this method, it is necessary to provide a raster stack or brick object with environmental variables and a value of the buffer width in m if raster (used in rlayer) has a longitude/latitude CRS, or map units in other cases. Usage method=c('GEO_ENV_KM_CONST', width='50000', env = somevar)
#' }
#'
#' @param rlayer raster. A raster layer used for sampling pseudo-absence
#' It is recommended to use a layer with the same resolution and extent that environmental variables that will be used for modeling. In the case use maskval argument, this raster layer must contain the values to sampling constraint
#' @param maskval integer or numeric. Values of the raster layer used for constraining the pseudo-absence sampling
#' @param calibarea raster or shapefile. Raster or shapefile with delimited calibration area used for a given species (see calibarea function).
#'
#' @return
#' @export
#'
#' @importFrom dplyr mutate
#' @importFrom raster extent extract cellStats match mask coordinates
#' @importFrom stats na.exclude kmeans
#'
#' @examples
#' \dontrun{
#' data("spp")
#' data("somevar")
#' data("regions")
#'
#' plot(regions)
#'
#' single_spp <-
#'   spp %>%
#'   dplyr::filter(species == "sp3") %>%
#'   dplyr::filter(pr_ab == 1) %>%
#'   dplyr::select(-pr_ab)
#' names(somevar) <- paste("var", 1:4)
#'
#'
#' # Pseudo-absences randomly sampled throughout study area
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 10,
#'     method = "RND",
#'     rlayer = regions,
#'     maskval = NULL
#'   )
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#'
#'
#' # Pseudo-absences randomly sampled within a regions where a species occurs
#'
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 10,
#'     method = "RND",
#'     rlayer = regions,
#'     maskval = samp_here
#'   )
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#'
#'
#' # Pseudo-absences sampled with environmental and geographical constraint
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = nrow(single_spp) * 10,
#'     method = c(method = "GEO_ENV_CONST", env = somevar, width = 50000),
#'     rlayer = regions,
#'     maskval = NULL
#'   )
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#'
#' # Pseudo-absences randomly sampled for regions where species occurs with environmental and geographical constraint
#' samp_here <- raster::extract(regions, single_spp[, c("x", "y")]) %>%
#'   unique() %>%
#'   na.exclude()
#'
#' ps1 <-
#'   sample_pseudoabs(
#'     data = single_spp,
#'     x = "x",
#'     y = "y",
#'     n = 1000,
#'     method = c(method = "GEO_ENV_CONST", env = somevar, width = 50000),
#'     rlayer = regions,
#'     maskval = samp_here
#'   )
#' plot(regions, col = gray.colors(9))
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' points(ps1, col = "red", cex = 0.7, pch = 19)
#' }
sample_pseudoabs <- function(data, x, y, n, method, rlayer, maskval = NULL, calibarea = NULL) {

  if(!any(c('RND', 'ENV_CONST', 'GEO_CONST', 'GEO_ENV_CONST', 'GEO_ENV_KM_CONST')%in%method)){
    stop("argument 'method' was misused, available methods RND, ENV_CONST, GEO_CONST, GEO_ENV_CONST, and GEO_ENV_KM_CONST")
  }

  rlayer <- rlayer[[1]]

  if (any(method %in% "RND")) {
    cell_samp <- sample_background(n = n, rlayer = rlayer, maskval = maskval)
  }

  if (any(method == "ENV_CONST")) {
    if (is.null(method["env"])) {
      stop("Provide a environmental stack/brick variables for ENV_CONST method, \ne.g. method = c('ENV_CONST', env=somevar)")
    }

    env <- method[["env"]]

    # Restriction for a given region
    envp <- inv_bio(e = env, p = data[, c(x, y)])
    if (raster::extent(env) != raster::extent(rlayer)) {
      rlayer2 <- rlayer
      rlayer2[] <- raster::extract(envp, coordinates(rlayer2))
      envp <- rlayer2
    }

    if (!is.null(maskval)) {
      rvalues <- raster::cellStats(rlayer, unique) %>% stats::na.exclude()
      filt <- raster::match(rlayer, maskval)
      filt[filt[] == 0] <- NA
      envp <- raster::mask(envp, filt)
    }

    cell_samp <- sample_background(n = n, rlayer = envp)
  }

  if (any(method == "GEO_CONST")) {
    if (!"width" %in% names(method)) {
      stop("Provide a width value for 'GEO_CONST' method, \ne.g. method=c('GEO_CONST', width='50000')")
    }

    # Restriction for a given region
    envp <- inv_geo(e = rlayer, p = data[, c(x, y)], d = as.numeric(method["width"]))

    if (!is.null(maskval)) {
      rvalues <- raster::cellStats(rlayer, unique) %>% stats::na.exclude()
      filt <- raster::match(rlayer, maskval)
      filt[filt[] == 0] <- NA
      envp <- raster::mask(envp, filt)
    }

    cell_samp <- sample_background(n = n, rlayer = envp)
  }

  if (any(method == "GEO_ENV_CONST")) {
    if (!all(c("env", "width") %in% names(method))) {
      stop("Provide a width value and environmental stack/brick variables for 'GEO_ENV_CONST' method, \ne.g. method=c('GEO_ENV_CONST', width='50000', env=somevar)")
    }

    # Restriction for a given region
    envp <- inv_geo(e = rlayer, p = data[, c(x, y)], d = as.numeric(method["width"]))
    envp2 <- inv_bio(e = env, p = data[, c(x, y)])
    if (raster::extent(env) != raster::extent(rlayer)) {
      rlayer2 <- rlayer
      rlayer2[] <- raster::extract(envp2, coordinates(rlayer2))
      envp2 <- rlayer2
    }
    envp <- (envp2 + envp)
    rm(envp2)

    if (!is.null(maskval)) {
      rvalues <- raster::cellStats(rlayer, unique) %>% stats::na.exclude()
      filt <- raster::match(rlayer, maskval)
      filt[filt[] == 0] <- NA
      envp <- raster::mask(envp, filt)
    }

    cell_samp <- sample_background(n = n, rlayer = envp)
  }


  if (any(method == "GEO_ENV_KM_CONST")) {
    if (!all(c("env", "width") %in% names(method))) {
      stop("Provide a width value and environmental stack/brick variables for 'GEO_ENV_KM_CONST' method, \ne.g. method=c('GEO_ENV_CONST', width='50000', env=somevar)")
    }

    # Restriction for a given region
    envp <- inv_geo(e = rlayer, p = data[, c(x, y)], d = as.numeric(method["width"]))
    envp2 <- inv_bio(e = env, p = data[, c(x, y)])
    if (raster::extent(env) != raster::extent(rlayer)) {
      rlayer2 <- rlayer
      rlayer2[] <- raster::extract(envp2, coordinates(rlayer2))
      envp2 <- rlayer2
    }
    envp <- (envp2 + envp)
    rm(envp2)

    if (!is.null(maskval)) {
      rvalues <- raster::cellStats(rlayer, unique) %>% stats::na.exclude()
      filt <- raster::match(rlayer, maskval)
      filt[filt[] == 0] <- NA
      envp <- raster::mask(envp, filt)
    }

    # K-mean procedure
    if (raster::extent(env) != raster::extent(envp)) {
      env_changed <- crop(env, envp)
    }
    env_changed <- mask(env_changed, envp)

    env_changed <- data.frame(raster::coordinates(env_changed), val = values(env_changed))
    env_changed <- stats::na.exclude(env_changed)

    suppressWarnings(km <- stats::kmeans(env_changed, centers = n))
    cell_samp <- km$centers[, 1:2] %>%data.frame()
    cell_samp <- cell_samp %>% dplyr::mutate(val = raster::extract(envp, cell_samp))
    cell_samp <- cell_samp[!is.na(cell_samp$val),-3]
  }

  return(cell_samp)
}
