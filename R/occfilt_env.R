#' Perform environmental filtering on species occurrences
#'
#' @description This function perform filtering on species occurrences based on their environmental conditions.
#'
#' @param data data.frame. Data.frame or tibble object with presences
#' (or presence-absence) records, and coordinates
#' @param x character. Column name with spatial x coordinates
#' @param y character. Column name with spatial y coordinates
#' @param id character. Column names with rows id. It is important that each row has its own unique code.
#' @param env_layer SpatRaster. Rasters with environmental conditions
#' @param nbins integer. A number of classes used to split each environmental condition
#'
#' @return
#' A tibble object with data environmentally filtered
#'
#' @details This function uses an approach adapted from the approach proposed by Varela et al. (2014).
#'  It consists of filtering occurrences in  environmental space. First, a regular
#'  multidimensional grid is created in  environmental space. This multidimensional
#'  grid is determined by the environmental variables (always use continuous variables) the
#'  grid cell size is defined by the number of bins, used for dividing variable range into
#'  interval classes (Varela et al. 2014; Castellanos et al., 2019). The number of bins is set in
#'  the "nbins" argument. Then, a single occurrence is randomly selected within each cell of the
#'  multivariate grid. Consider that there is a trade-off between the number of bins and the number
#'  of filtered records because as the number of bins decreases, the cell size of the grids
#'  increases, and the number of filtered records decreases (Castellanos et al., 2019).
#'  occfilt_env works for any number of dimensions (variables) and with the original variables without performing
#'  a PCA beforehand.
#'
#'  The greater the number of predictor variables (i.e., the number of dimensions of the
#'  multidimensional environmental grid) and the greater the number of bins, the greater the time processing
#'  and the computer memory used. Therefore, it is recommended to use a small number of bins
#'  between 2-5 if more than ten variables are used.
#'
#' @references
#' \itemize{
#' \item Castellanos, A. A., Huntley, J. W., Voelker, G., & Lawing, A. M. (2019). Environmental
#' filtering improves ecological niche models across multiple scales. Methods in Ecology and
#' Evolution, 10(4), 481-492. https://doi.org/10.1111/2041-210X.13142
#'
#' \item Varela, S., Anderson, R. P., Garcia-Valdes, R., & Fernandez-Gonzalez, F. (2014).
#' Environmental filters reduce the effects of sampling bias and improve predictions of
#' ecological niche models. Ecography, 37, 1084-1091.
#' https://doi.org/10.1111/j.1600-0587.2013.00441.x
#' }
#'
#' @export
#'
#' @importFrom dplyr %>% mutate select starts_with pull tibble
#' @importFrom stats complete.cases
#' @importFrom terra extract nlyr levels classify
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
#' spp1$idd <- 1:nrow(spp1)
#'
#'
#' # split environmental variables into 5 bins
#' filtered_1 <- occfilt_env(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   id = "idd",
#'   env_layer = somevar,
#'   nbins = 5
#' )
#'
#' # split into 8 bins
#' filtered_2 <- occfilt_env(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   id = "idd",
#'   env_layer = somevar,
#'   nbins = 8
#' )
#'
#' # split into 12 bins
#' filtered_3 <- occfilt_env(
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   id = "idd",
#'   env_layer = somevar,
#'   nbins = 12
#' )
#' # note that the higher the nbins parameter the more
#' # classes must be processed (4 variables, 30 bins = 923521 classes)
#'
#' # While the greater the number of bins the fewer records retained
#' }
#'
#' @seealso \code{\link{occfilt_geo}}
#'
occfilt_env <- function(data, x, y, id, env_layer, nbins) {
  s <- . <- l <- NULL

  da <- data[c(x, y, id)]
  coord <- data[c(x, y)]

  message("Extracting values from raster ...")
  env <- env_layer
  env_layer <- terra::extract(env_layer, coord)
  env_layer$ID <- NULL

  filt <- stats::complete.cases(env_layer)
  if (sum(!filt) > 0) {
    message(sum(!filt), " records were removed because they have NAs for some variables")
    da <- da[filt, ]
    coord <- coord[filt, ]
    env_layer <- env_layer[filt, ]
  }
  rm(filt)

  n <- ncol(env_layer)
  res <- (apply(env_layer, 2, max) - apply(env_layer, 2, min)) / nbins

  classes <- list()
  for (i in 1:n) {
    ext1 <- range(env_layer[, i])
    ext1[2] <- ext1[2] + 0.05
    ext1 <- c(ext1[1] - 0.000001, ext1[2] + 0.000001)
    classes[[i]] <- seq(ext1[1], ext1[2], by = res[i])
    classes[[i]][length(classes[[i]])] <- ext1[2]
  }

  for (i in 1:(terra::nlyr(env))) {
    env[[i]] <- terra::classify(env[[i]], classes[[i]], include.lowest = TRUE)
    levels(env[[i]]) <- as.numeric(as.factor(terra::levels(env[[i]])[[1]]))
  }
  real_p <- terra::extract(env, coord)[-1]
  real_p$groupID <- apply(real_p, 1, function(x) paste(x, collapse = "."))

  message("Number of unfiltered records: ", nrow(da))

  if (any(is.na(real_p$groupID))) {
    nas <- da[is.na(real_p$groupID), c(id, x, y)]
    no_nas <- da[!duplicated(real_p$groupID) & !is.na(real_p$groupID), c(id, x, y)]
    coord_filter <- unique(dplyr::bind_rows(no_nas, nas))
  } else {
    coord_filter <- da[!duplicated(real_p$groupID), c(id, x, y)]
  }

  message("Number of filtered records: ", nrow(coord_filter))
  return(dplyr::tibble(coord_filter))
}
