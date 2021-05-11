#' Transform a spatial partition layer to the same spatial properties of environmental variables
#'
#' @details Transform a layer originated from the function block_partition or band_partition to the same properties of environmental variables
#'
#' @param env_layer raster. A raster, stack, or brick object with some environmental variables used in the block_partition or band_partition function. Function always will select the first layer
#' @param best_grid raster. Raster object returned by block_partition or band_partition
#'
#' @return
#' @export
#'
#' @importFrom dplyr select
#' @importFrom raster coordinates ncell values extract
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' data(spp)
#' f <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(f)
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
#' grid_env
#' part$grid
#'
#' plot(part$grid)
#' plot(grid_env)
#' }
#'
get_block <- function(env_layer, best_grid) {
  maskr <- env_layer[[1]]
  rdf <- terra::as.data.frame(maskr, xy=TRUE)
  val <- terra::extract(best_grid, rdf[1:2])[,2]
  maskr[as.numeric(rownames(rdf))] <- val
  return(maskr)
}
