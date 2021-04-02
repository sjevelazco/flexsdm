#' Plot the maximum resolution to be used in block_partition_pa  
#'
#' @param r raster. A raster layer. Preferably a layer of environmental variables to be used  
#' @param max_res_mult numeric. Maximum resolution to be tested.
#'
#' @importFrom raster res plot rasterToPolygons
#'
#' @return A plot with raster and grid with the resolution setted 
#' @export
#'
#' @examples
#' \dontrun{
#' f <- system.file("external/test.grd", package="raster")
#' f
#' r <- raster(f)
#' plot_max_res(r, max_res_mult = 2)
#' plot_max_res(r, max_res_mult = 10)
#' plot_max_res(r, max_res_mult = 20)
#' }
plot_max_res <- function(r, max_res_mult){
  r0 <- r
  r0[!is.na(values(r0))] <- 1
  
  r[] <- 0
  raster::res(r) <- raster::res(r)*max_res_mult 
  raster::plot(r0)
  r <- raster::rasterToPolygons(r)
  raster::plot(r, add=TRUE)
}
