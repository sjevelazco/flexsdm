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
#' \item thickening: Thickening background points is based on Vollering et al. (2019) method. For this method, it is necessary to define a buffer width that will be used around presences points. It possible to define a buffer using the argument as method = c("thickening", width = 20000). Buffer width must be in m if raster (used in rlayer) has a longitude/latitude CRS, or map units in other cases. If a width buffer is not provided function will use a width value equal to the mean distance of the pair-wise presences distance. In case of a width value is not provided, the argument must be used as method='thickening'.
#' }
#' Usage method='thickening' or method = c("thickening", width = 20000). Default 'random'
#' @param rlayer SpatRaster used for sampling background points.
#' It is recommended to use a layer with the same resolution and extent that environmental variables that will be used for modeling. In the case use maskval argument, this raster layer must contain the values to sampling constraint
#' @param maskval integer or numeric. Values of the raster layer used for constraining the background points sampling
#' @param calibarea SpatVect that delimits the calibration area used for a given species (see calib_area function).
#' @param rbias SpatRaster used for bias background points. When use 'biased' method must be provided the raster with
#' with bias data. It is recommended that rbias match resolution and extent of rlayer.
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
#' ## %######################################################%##
#' #                                                          #
#' ####             Random background method               ####
#' #                                                          #
#' ## %######################################################%##
#'
#' # Sample background points throughout study area with random method
#' spp_p <- spp_pa %>% dplyr::filter(pr_ab == 1)
#' bg <-
#'   sample_background(
#'     data = spp_p,
#'     x = "x",
#'     y = "y",
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
#'   x = "x",
#'   y = "y",
#'   n = 1000,
#'   method = "random",
#'   rlayer = grid_env,
#'   maskval = 1
#' ) %>% points()
#'
#' plot(grid_env)
#' sample_background(
#'   data = spp_p,
#'   x = "x",
#'   y = "y",
#'   n = 1000,
#'   method = "random",
#'   rlayer = grid_env,
#'   maskval = 2
#' ) %>% points()
#'
#' plot(grid_env)
#' sample_background(
#'   data = spp_p,
#'   x = "x",
#'   y = "y",
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
#'   x = "x",
#'   y = "y",
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
#' ## %######################################################%##
#' #                                                          #
#' ####            Thickening background method            ####
#' #                                                          #
#' ## %######################################################%##
#'
#' # Thickening background without containing them
#' spp_p # presences database of a species
#' grid_env # The raster layer used for sampling background
#' bg <- sample_background(
#'   data = spp_p,
#'   x = "x",
#'   y = "y",
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
#' # Thickening background without using a given buffer width
#' spp_p # presences database of a species
#' grid_env # The raster layer used for sampling background
#' bg <- sample_background(
#'   data = spp_p,
#'   x = "x",
#'   y = "y",
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
#'   x = "x",
#'   y = "y",
#'   n = 3000,
#'   method = "thickening",
#'   rlayer = grid_env,
#'   maskval = 2,
#'   calibarea = ca_ps1
#' )
#'
#' plot(grid_env)
#' plot(ca_ps1, add = T)
#' bg %>%
#'   points(col = "red", cex = 0.3)
#' points(spp_p[c("x", "y")], pch = 19)
#'
#' ## %######################################################%##
#' #                                                          #
#' ####             Biased background method               ####
#' #                                                          #
#' ## %######################################################%##
#' require(dplyr)
#' require(terra)
#' data(spp)
#'
#' # Lets select the presences of a species
#' spp_p <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)
#'
#' # Raster layer with density of poinst to obtain a biased sampling background
#' occ_density <- system.file("external/occ_density.tif", package = "flexsdm")
#' occ_density <- terra::rast(occ_density)
#' plot(occ_density)
#' points(spp_p %>% dplyr::select(x, y), cex = 0.5)
#'
#' # A layer with region used to contraint background
#' regions <- system.file("external/regions.tif", package = "flexsdm")
#' regions <- terra::rast(regions)
#' plot(regions)
#' points(spp_p %>% dplyr::select(x, y), cex = 0.5)
#'
#'
#' # Biased background points
#' spp_p # presences database of a species
#' bg <- sample_background(
#'   data = spp_p,
#'   x = "x",
#'   y = "y",
#'   n = 3000,
#'   method = "biased",
#'   rlayer = regions,
#'   rbias = occ_density
#' )
#'
#' plot(occ_density)
#' bg %>%
#'   points(col = "red", cex = 0.1)
#' spp_p %>%
#'   dplyr::select(x, y) %>%
#'   points(., col = "black", pch = 19, cex = 0.5)
#'
#'
#' # Biased background points constrained in a region
#' # It will be selected region 6
#' plot(regions)
#' plot(regions %in% c(1, 6))
#'
#' bg <- sample_background(
#'   data = spp_p,
#'   x = "x",
#'   y = "y",
#'   n = 500,
#'   method = "biased",
#'   rlayer = regions,
#'   rbias = occ_density,
#'   maskval = c(1, 2)
#' )
#'
#' plot(occ_density)
#' bg %>%
#'   points(col = "red", cex = 0.5)
#' spp_p %>%
#'   dplyr::select(x, y) %>%
#'   points(., col = "black", pch = 19, cex = 0.5)
#' }
sample_background <-
  function(data,
           x,
           y,
           n,
           method = "random",
           rlayer,
           maskval = NULL,
           calibarea = NULL,
           rbias = NULL) {
    if (!method[1] %in% c("random", "thickening", "biased")) {
      stop("argument 'method' was misused, available methods 'random', 'thickening'")
    }

    if (method[1] %in% c("biased") & is.null(rbias)) {
      stop("for using 'random' method a raster layer with biases data must be provided in 'rbias' argument")
    }

    # Prepare datasets
    if (class(rlayer)[1] != "SpatRaster") {
      rlayer <- terra::rast(rlayer)
    }

    rlayer <- rlayer[[1]]
    data <- data[, c(x, y)]

    # Remove cell with presences
    rlayer[na.omit(terra::cellFromXY(rlayer, as.matrix(data)))] <- NA

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

    # Correct rbias data in case it don't match resolution or extent of rlayer
    if (method[1] %in% c("biased")) {
      if (any(!(ext(rlayer)[1:4] %in% ext(rbias)[1:4])) | all(!res(rlayer) %in% res(rbias))) {
        rbias2 <- rlayer
        terra::values(rbias2) <- NA
        df <- terra::as.data.frame(rlayer, xy = TRUE)
        rbias2[as.numeric(rownames(df))] <-
          terra::extract(rbias, df[c("x", "y")])[, 2]
        rbias <- rbias2
        rm(rbias2, df)
      }
    }

    if (method[1] %in% c("biased")) {
      rlayer <- mask(rlayer, rbias)
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
      } else if (any(method == "thickening")) {
        cell_samp <-
          sample(cell_samp,
            size = n,
            replace = FALSE,
            prob = buf_r[cell_samp][, 1]
          )
      } else if (any(method == "biased")) {
        cell_samp <-
          sample(cell_samp,
            size = n,
            replace = FALSE,
            prob = rbias[cell_samp][, 1]
          )
      }

      cell_samp <- terra::xyFromCell(rlayer, cell_samp) %>%
        data.frame() %>%
        dplyr::tibble()
    }
    colnames(cell_samp) <- c("x", "y")
    return(cell_samp)
  }
