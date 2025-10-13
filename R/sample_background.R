#' Sample background points
#'
#' @description Sampling background points with the options of
#' using different geographical restrictions and sampling methods.
#'
#' @param data data.frame or tibble. Database with presences records, and coordinates
#' @param x character. Column name with spatial x coordinates
#' @param y character. Column name with spatial y coordinates
#' @param n integer. Number of background point to be sampled
#' @param method character. Background allocation method. The methods implemented are:
#' \itemize{
#' \item random: Random allocation of background points. Usage method = 'random'
#' \item thickening: Thickening background points is based on the Vollering et al. (2019) method.
#' For this method, a buffer width must be defined that will be used around presences points.
#' A buffer can be defined using the argument as method = c("thickening", width = 20000).
#' Buffer width must be in m if raster (used in rlayer) has a longitude/latitude CRS, or map
#' units in other cases. If a buffer width is not provided the function will use a width value
#' equal to the mean of the pair-wise presence distances. If a width value is not provided, the
#' argument must be used as method = 'thickening'.
#' \item biased: This method, similar to "thickening", sample background biased with the same
#' bias of presences. However, the background points are sampled used a presences probability
#' throughout the entire study area, and not restricting such bias within buffers as in the
#' “thickening” approach. For using this method, it is necessary to provide a layer with presences
#' bias in "rbias" argument (Phillips et al., 2009).
#' }
#' Usage method='thickening' or method = c("thickening", width = 20000). Default 'random'
#' @param rlayer SpatRaster used for sampling background points.
#' It is best to use a layer with the same resolution and extent that environmental variables that
#' will be used for modeling. If using maskval argument, this raster layer must contain the values
#'  to constrain sampling
#' @param maskval integer, character, or factor. Values of the raster layer used for constraining sampling of
#' background points
#' @param calibarea SpatVect that delimits the calibration area used for a given species
#' (see calib_area function).
#' @param rbias SpatRaster used for choosing background points using the bias method. A raster
#' with bias data must be provided. It is recommended that rbias match resolution and extent
#' of rlayer.
#' @param sp_name character. Species name for which the output will be used.
#' If this argument is used, the first output column will have the species name. Default NULL.
#'
#' @references
#' \itemize{
#' \item Phillips, S. J., Dudík, M., Elith, J., Graham, C. H., Lehmann, A., Leathwick, J., &
#' Ferrier, S. (2009). Sample selection bias and presence-only distribution models:
#' Implications for background and pseudo-absence data. Ecological Applications, 19(1), 181-197.
#'
#' \item Vollering, J., Halvorsen, R., Auestad, I., & Rydgren, K. (2019). Bunching up the
#' background betters bias in species distribution models. Ecography, 42(10), 1717-1727.
#' https://doi.org/10.1111/ecog.04503
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
#' require(terra)
#' require(dplyr)
#' data(spp)
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' # Example for a single species
#' spp_pa <- spp %>% dplyr::filter(species == "sp3")
#'
#' # Spatially structured partition
#' part <- part_sblock(
#'   env_layer = somevar,
#'   data = spp_pa,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   min_res_mult = 100,
#'   max_res_mult = 500,
#'   num_grids = 30,
#'   min_occ = 5,
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
#' # Sample background points throughout study area with random sampling method
#' spp_p <- spp_pa %>% dplyr::filter(pr_ab == 1)
#' bg <-
#'   sample_background(
#'     data = spp_p,
#'     x = "x",
#'     y = "y",
#'     n = 1000,
#'     method = "random",
#'     rlayer = grid_env,
#'     sp_name = "sp3"
#'   )
#'
#' bg
#' plot(grid_env)
#' points(bg[-1])
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
#'   crs = crs(somevar)
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
#' ## %######################################################%##
#' #                                                          #
#' ####            Thickening background method            ####
#' #                                                          #
#' ## %######################################################%##
#'
#' # Thickening background without constraining them
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
#' # Thickening background
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
#' # Select the presences of a species
#' spp_p <- spp %>% dplyr::filter(species == "sp1", pr_ab == 1)
#'
#' # Raster layer with density of points to obtain a biased sampling background
#' occ_density <- system.file("external/occ_density.tif", package = "flexsdm")
#' occ_density <- terra::rast(occ_density)
#' plot(occ_density)
#' points(spp_p %>% dplyr::select(x, y), cex = 0.5)
#'
#' # A layer with region used to contrain background sampling area
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
#' # Biased background points constrained to a region
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
           rbias = NULL,
           sp_name = NULL) {
    . <- NULL
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
    if (!is.null(rbias)) {
      if (class(rbias)[1] != "SpatRaster") {
        rbias <- terra::rast(rbias)
      }
    }

    rlayer <- rlayer[[1]]
    data <- data[, c(x, y)]

    # Remove cell with presences
    rlayer[na.omit(terra::cellFromXY(rlayer, as.matrix(data)))] <- NA

    # Mask to calibration area
    if (!is.null(calibarea)) {
      rlayer <- rlayer %>%
        terra::crop(., calibarea) %>%
        terra::mask(., calibarea)
    }

    # Mask to maksvalue
    if (!is.null(maskval)) {
      if (is.factor(rlayer)) {
        maskval <-
          which(levels(rlayer)[[1]][, 2] %in% as.character(maskval))
        rlayer <- rlayer * 1
      }
      filt <- terra::match(rlayer, maskval)
      rlayer <- terra::mask(rlayer, filt)
    }

    # Correct rbias data in case it don't match resolution or extent of rlayer
    if (method[1] %in% c("biased")) {
      if (any(!(ext(rlayer)[1:4] %in% ext(rbias)[1:4])) | all(!res(rlayer) %in% res(rbias))) {
        if (!all(res(rlayer) %in% res(rbias))) {
          rbias <- terra::resample(rbias, rlayer, method = "bilinear")
        }
        rbias2 <- rbias %>%
          terra::crop(., rlayer) %>%
          terra::mask(., rlayer)
        # df <- terra::as.data.frame(rlayer, xy = TRUE)
        # rbias2[] <- NA
        # rbias2[as.numeric(rownames(df))] <-
        #   terra::extract(rbias, df[c("x", "y")])[, 2]
        # rn(df)
        rbias <- rbias2
        rm(rbias2)
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
    colnames(cell_samp) <- c(x, y)
    cell_samp$pr_ab <- 0
    if (!is.null(sp_name)) {
      cell_samp <- tibble(sp = sp_name, cell_samp)
    }
    return(cell_samp)
  }
