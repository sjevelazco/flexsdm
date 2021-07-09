#' Sample background points
#'
#' @description Sampling background points with the possibility of
#' using different geographical restrictions and sampling method.
#'
#' @param data data.frame or tibble. Database with presences records, and coordinates
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param n integer. Number of background point to be sampled
#' @param method character. Background allocation method. The methods implemented are:
#' \itemize{
#' \item random: Random allocation of background points. Usage method = 'random'.
#' \item thickening: Thickening background point based on Vollering et al. (2019) method. For method, it is necessary to define a buffer width that will be used around presences points. It possible to define a buffer using the argument as method = c("thickening", width = 20000). Buffer width must be in m if raster (used in rlayer) has a longitude/latitude CRS, or map units in other cases. If a width buffer is not provided function will use a width value equal to the mean distance of the pair-wise presences distance. In case of a width value is not provided, the argument must be used as method='thickening'.
#' }
#' Usage method='thickening' or method = c("thickening", width = 20000). Default 'random'
#' @param rlayer SpatRaster used for sampling background points.
#' It is recommended to use a layer with the same resolution and extent that environmental variables that will be used for modeling. In the case use maskval argument, this raster layer must contain the values to sampling constraint
#' @param maskval integer or numeric. Values of the raster layer used for constraining the background points sampling
#' @param calibarea SpatVect that delimits the calibration area used for a given species (see calib_area function).
#'
#' @references
#' \itemize{
#' \item Vollering, J., Halvorsen, R., Auestad, I., & Rydgren, K. (2019). Bunching up the background betters bias in species distribution models. Ecography, 42(10), 1717-1727. https://doi.org/10.1111/ecog.04503
#' }
#'
#' @return
#' A tibble object with x y coordinates of sampled background points
#'
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
#' spp_pa <- spp %>% dplyr::filter(species == "sp3")
#'
#' #
#' part <- part_sblock(
#'   env_layer = somevar,
#'   data = spp_pa,
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
#'
#' ##%######################################################%##
#' #                                                          #
#' ####             Random background method               ####
#' #                                                          #
#' ##%######################################################%##
#'
#' # Sample background points throughout study area with random method
#' spp_p <- spp_pa %>% dplyr::filter(pr_ab==1)
#' bg <-
#'   sample_background(
#'     data = spp_p,
#'     x = 'x',
#'     y = 'y',
#'     n = 1000,
#'     method = "random",
#'     rlayer = grid_env
#'   )
#' plot(grid_env)
#' points(bg)
#'
#' # Sample random background points constrained to a region with a give set of values
#' plot(grid_env)
#' sample_background(
#'   data = spp_p,
#'   x = 'x',
#'   y = 'y',
#'   n = 1000,
#'   method = "random",
#'   rlayer = grid_env,
#'   maskval = 1
#' ) %>% points()
#'
#' plot(grid_env)
#' sample_background(
#'   data = spp_p,
#'   x = 'x',
#'   y = 'y',
#'   n = 1000,
#'   method = "random",
#'   rlayer = grid_env,
#'   maskval = 2
#' ) %>% points()
#'
#' plot(grid_env)
#' sample_background(
#'   data = spp_p,
#'   x = 'x',
#'   y = 'y',
#'   n = 1000,
#'   method = "random",
#'   rlayer = grid_env,
#'   maskval = c(1, 2)
#' ) %>% points()
#'
#' # Sample random background within a calibration area and constrained to a region
#' ca_ps1 <- calib_area(
#'   data = spp_pa,
#'   x = "x",
#'   y = "y",
#'   method = c("buffer", width = 50000),
#' )
#' plot(grid_env)
#' plot(ca_ps1, add = T)
#' points(spp_pa[-1], col = "blue", cex = 0.7, pch = 19)
#' sample_background(
#'   data = spp_p,
#'   x = 'x',
#'   y = 'y',
#'   n = 1000,
#'   method = "random",
#'   rlayer = grid_env,
#'   maskval = 1,
#'   calibarea = ca_ps1
#' ) %>%
#'   points(col = "red")
#'
#'
#'
#' ##%######################################################%##
#' #                                                          #
#' ####            Thickening background method            ####
#' #                                                          #
#' ##%######################################################%##
#'
#' # Thickening background without containing them
#' spp_p # presences database of a species
#' grid_env # The raster layer used for sampling background
#' bg <- sample_background(
#'   data = spp_p,
#'   x = 'x',
#'   y = 'y',
#'   n = 5000,
#'   method = "thickening",
#'   rlayer = grid_env,
#' )
#'
#' plot(grid_env)
#' bg %>%
#'   points(col = "red")
#'
#'
#' # Thickening background without usin a given buffer width
#' spp_p # presences database of a species
#' grid_env # The reaster layer used for sampling background
#' bg <- sample_background(
#'   data = spp_p,
#'   x = 'x',
#'   y = 'y',
#'   n = 5000,
#'   method = c("thickening", width = 150000),
#'   rlayer = grid_env
#' )
#'
#' plot(grid_env)
#' bg %>%
#'   points(col = "red")
#'
#' # Sample thickening background within a calibration area and constrained to a region
#' bg <- sample_background(
#'   data = spp_p,
#'   x = 'x',
#'   y = 'y',
#'   n = 3000,
#'   method = "thickening",
#'   rlayer = grid_env,
#'   maskval = 2,
#'   calibarea = ca_ps1
#' )
#'
#' plot(grid_env)
#' plot(ca_ps1, add=T)
#' bg %>%
#'   points(col = "red", cex=0.3)
#' points(spp_p[c('x', 'y')], pch=19)
#' }
#'
sample_background <-
  function(data,
           x,
           y,
           n,
           method = "random",
           rlayer,
           maskval = NULL,
           calibarea =  NULL) {
  if (!method[1] %in% c("random", "thickening")) {
    stop("argument 'method' was misused, available methods 'random', 'thickening'")
  }

  # Prepare datasets
  if (class(rlayer) != "SpatRaster") {
    rlayer <- terra::rast(rlayer)
  }

  rlayer <- rlayer[[1]]
  data <- data[, c(x, y)]

  # Remove cell with presences
  rlayer[terra::cellFromXY(rlayer, as.matrix(data))] <- NA

  # Mask to calibration area
  if (!is.null(calibarea)) {
    rlayer <- terra::mask(rlayer, calibarea)
  }

  # Mask to maksvalue
  if (!is.null(maskval)) {
    if (is.factor(maskval)) {
      maskval <-
        which(levels(maskval)[-1] %in% as.character(maskval))
      rlayer <- rlayer * 1
    }
    filt <- terra::match(rlayer, maskval)
    rlayer <- terra::mask(rlayer, filt)
  }

  ncellr <- terra::global(!is.na(rlayer), sum)

  # Create density of buffers sum
  if (any(method == "thickening")) {
    data2 <- terra::vect(data, geom = c(x, y), crs = crs(rlayer))
    if (is.na(method["width"])) {
      buf_with <- mean(terra::distance(data2))
    } else {
      buf_with <- as.numeric(method["width"])
    }
    buf <- terra::buffer(data2, buf_with, quadsegs = 10)
    buf_r <- !is.na(rasterize(buf[1], rlayer))
    for (i in 2:nrow(buf)) {
      buf_r <- buf_r + !is.na(rasterize(buf[i], rlayer))
    }
    buf_r <- terra::mask(buf_r, rlayer)
  }

  if (ncellr < n) {
    message(
      "Number of background-points exceeds number of cell will be returned ",
      ncellr,
      " background-points"
    )
    cell_samp <- terra::as.data.frame(rlayer, na.rm = TRUE, cells = TRUE)[, "cell"]
    cell_samp <- terra::xyFromCell(rlayer, cell_samp) %>%
      data.frame() %>%
      dplyr::tibble()
  } else {
    cell_samp <- terra::as.data.frame(rlayer, na.rm = TRUE, cells = TRUE)[, "cell"]

    if (any(method == "random")) {
      cell_samp <-
        sample(cell_samp,
          size = n,
          replace = FALSE,
          prob = NULL
        )
    } else {
      cell_samp <-
        sample(cell_samp,
          size = n,
          replace = FALSE,
          prob = buf_r[cell_samp][, 1]
        )
    }

    cell_samp <- terra::xyFromCell(rlayer, cell_samp) %>%
      data.frame() %>%
      dplyr::tibble()
  }
  colnames(cell_samp) <- c("x", "y")
  return(cell_samp)
}
