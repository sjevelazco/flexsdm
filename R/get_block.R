#' Transform a spatial partition layer to the same spatial properties as environmental variables
#'
#' @details Transform a layer originating from the function block_partition or band_partition to the same spatial properties as the environmental variables
#'
#' @param env_layer SpatRaster object with some environmental variables used in the block_partition or band_partition function. Function always will select the first layer
#' @param best_grid SpatRaster object returned by block_partition or band_partition
#'
#' @return A SpatRaster layer with the same resolution and extent as the environmental variables
#' @export
#'
#' @importFrom dplyr %>% select
#' @importFrom terra as.data.frame extract
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' require(terra)
#' data(spp)
#' f <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(f)
#'
#' # Example for a single species
#' single_spp <- spp %>% dplyr::filter(species == "sp3")
#'
#' part <- part_sblock(
#'   env_layer = somevar,
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   min_res_mult = 100,
#'   max_res_mult = 500,
#'   num_grids = 10,
#'   min_occ = 5,
#'   n_part = 2
#' )
#'
#' grid_env <- get_block(env_layer = somevar, best_grid = part$grid)
#' grid_env
#' part$grid
#'
#' plot(part$grid)
#' plot(grid_env)
#' }
get_block <- function(env_layer, best_grid) {
  maskr <- env_layer[[!is.factor(env_layer)]][[1]]
  rdf <- terra::as.data.frame(maskr, xy = TRUE)
  val <- terra::extract(best_grid, rdf[1:2])[, 2]
  maskr[as.numeric(rownames(rdf))] <- val
  names(maskr) <- ".part"
  return(maskr)
}
