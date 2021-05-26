#' Sample background points
#'
#' @param n integer. Number of background point to be sampled
#' @param rlayer raster. A raster layer used for sampling background-point.
#' It is recommended to use a layer with the same resolution and extent that environmental variables that will be used for modeling. In the case use maskval argument, this raster layer must contain the values to sampling constraint
#' @param maskval integer or numeric. Values of the raster layer used for constraining the background points sampling
#' @param calibarea shapefile. A SpatialPolygon or SpatialPolygonDataFrame which delimit the calibration area used for a given species (see calib_area function).
#'
#' @return
#' @export
#'
#' @seealso \code{\link{sample_pseudoabs}} and \code{\link{calib_area}}.
#'
#' @importFrom dplyr %>% tibble
#' @importFrom stats na.exclude
#' @importFrom terra mask freq match values ncell xyFromCell
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' data(spp)
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' # Lest practice with a single species
#' single_spp <- spp %>% dplyr::filter(species == "sp3")
#'
#' part <- block_partition(
#'   env_layer = somevar,
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   min_res_mult = 100,
#'   max_res_mult = 500,
#'   num_grids = 10,
#'   n_part = 2
#' )
#'
#' grid_env <- get_block(env_layer = somevar, best_grid = part$grid)
#' plot(grid_env)
#'
#' # Sample background points throughout study area
#' bg <- sample_background(n = 1000, rlayer = grid_env)
#' plot(grid_env)
#' points(bg)
#'
#' # Sample background points constrained to a region regions with a give set of values
#' plot(grid_env)
#' sample_background(n = 1000, rlayer = grid_env, maskval = 1) %>% points()
#' plot(grid_env)
#' sample_background(n = 1000, rlayer = grid_env, maskval = 2) %>% points()
#' plot(grid_env)
#' sample_background(n = 1000, rlayer = grid_env, maskval = c(1, 2)) %>% points()
#'
#' # Sample background within a calibration area and constrained to a region regions
#' # Delimit a calibration area with calib_area
#' ca_ps1 <- calib_area(
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   method = c("buffer", width = 50000),
#' )
#' plot(grid_env)
#' plot(ca_ps1, add = T)
#' points(single_spp[-1], col = "blue", cex = 0.7, pch = 19)
#' sample_background(n = 1000, rlayer = grid_env, maskval = 1, calibarea = ca_ps1) %>%
#'   points(col = "red")
#' }
#'
sample_background <- function(n, rlayer, maskval = NULL, calibarea = NULL) {
  rlayer <- rlayer[[1]]

  if (!is.null(calibarea)) {
    rlayer <- terra::mask(rlayer, calibarea)
  }

  if (!is.null(maskval)) {
    if (is.factor(maskval)) {
      maskval <-
        which(levels(maskval)[-1] %in% as.character(maskval))
      rlayer <- rlayer * 1
    }
    filt <- terra::match(rlayer, maskval)
    rlayer <- terra::mask(rlayer, filt)
  }

  ncellr <- sum(!is.na(terra::values(rlayer)))

  if (ncellr < n) {
    message(
      "Number of background-points exceeds number of cell will be returned ",
      ncellr,
      " background-points"
    )
    cell_samp <- 1:terra::ncell(rlayer)
    cell_samp <- cell_samp[!is.na(terra::values(rlayer))]
    cell_samp <- terra::xyFromCell(rlayer, cell_samp) %>%
      data.frame() %>%
      dplyr::tibble()
  } else {
    cell_samp <- 1:terra::ncell(rlayer)
    cell_samp <- cell_samp[!is.na(terra::values(rlayer))]
    cell_samp <-
      sample(cell_samp,
        size = n,
        replace = FALSE,
        prob = NULL
      )
    cell_samp <- terra::xyFromCell(rlayer, cell_samp) %>%
      data.frame() %>%
      dplyr::tibble()
  }
  colnames(cell_samp) <- c("x", "y")
  return(cell_samp)
}
