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
#' # See examples in block_partition
#' }
#'
get_block <- function(env_layer, best_grid) {
  maskr <- env_layer[[1]]
  rdf <- data.frame(raster::coordinates(maskr),
    ncell = 1:raster::ncell(maskr),
    val = raster::values(maskr)
  ) %>%
    stats::na.exclude() %>%
    dplyr::select(-val)
  val <- raster::extract(best_grid, rdf[1:2])
  maskr[rdf$ncell] <- c(val)
  return(maskr)
}
