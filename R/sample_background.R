#' Sample background points
#'
#' @param n integer. Number of background point to be sampled
#' @param rlayer raster. A raster layer used for sampling background-point.
#' It is recommended to use a layer with the same resolution and extent that environmental variables that will be used for modeling. In the case use maskval argument, this raster layer must contain the values to sampling constraint
#' @param maskval integer or numeric. Values of the raster layer used for constraining the background points sampling
#'
#' @return
#' @export
#'
#' @importFrom dplyr tibble
#' @importFrom raster values mask ncell xyFromCell
#' @importFrom stats na.exclude
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' data(spp)
#' data(somevar)
#'
#' # Lest practice with a single species
#' single_spp <- spp %>% dplyr::filter(species == "sp3")
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
#' grid_env <- get_block(env_layer = somevar, best_grid = part$Grid)
#' plot(grid_env)
#'
#'
#' # Sample background points throughout study area
#' bg <- sample_background(n = 1000, rlayer = grid_env)
#' bg %>% plot()
#'
#' # Sample background points constrained to a region regions with a give set of values
#' sample_background(n = 1000, rlayer = grid_env, maskval = 1) %>% plot()
#' sample_background(n = 1000, rlayer = grid_env, maskval = 2) %>% plot()
#' sample_background(n = 1000, rlayer = grid_env, maskval = c(1, 2)) %>% plot()
#' }
#'
#' @importFrom dplyr tibble
#' @importFrom raster cellStats match mask ncell xyFromCell
#' @importFrom stats na.exclude
#'
#' @examples
sample_background <- function(n, rlayer, maskval = NULL) {
  rlayer <- rlayer[[1]]

  if (!is.null(maskval)) {
    rvalues <- raster::cellStats(rlayer, unique) %>% stats::na.exclude()
    filt <- raster::match(rlayer, maskval)
    filt[filt[] == 0] <- NA
    rlayer <- raster::mask(rlayer, filt)
  }

  ncellr <- sum(!is.na(rlayer[]))

  if (ncellr < n) {
    message(
      "Number of background-points exceeds number of cell will be returned ",
      ncellr,
      " background-points"
    )
    cell_samp <- 1:raster::ncell(rlayer)
    cell_samp <- cell_samp[!is.na(rlayer[])]
    cell_samp <- raster::xyFromCell(rlayer, cell_samp) %>%
      data.frame() %>%
      dplyr::tibble()
  } else {
    cell_samp <- 1:raster::ncell(rlayer)
    cell_samp <- cell_samp[!is.na(rlayer[])]
    cell_samp <-
      sample(cell_samp,
        size = n,
        replace = FALSE,
        prob = NULL
      )
    cell_samp <- raster::xyFromCell(rlayer, cell_samp) %>%
      data.frame() %>%
      dplyr::tibble()
  }
  return(cell_samp)
}
