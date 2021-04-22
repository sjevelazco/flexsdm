#' Plot different resolutions to be used in block_partition
#'
#' @description This function is useful to display the maximum and minimum resolution that you want to test with the block_partition function. Note that if the resolution to be tested is very fine, the plot display may take a long time.
#'
#' @param r raster. A raster layer. Preferably a layer of environmental variables to be used
#' @param res_mult numeric. Maximum or minimum resolution to be tested.
#'
#' @importFrom raster res plot rasterToPolygons
#'
#' @return A plot with the original raster overlapped by a grid with the resolution used
#' @export
#'
#' @importFrom raster res plot rasterToPolygons
#'
#' @examples
#' \dontrun{
#' f <- system.file("external/test.grd", package = "raster")
#' f
#' r <- raster(f)
#' plot_res(r, max_res_mult = 2)
#' plot_res(r, max_res_mult = 10)
#' plot_res(r, max_res_mult = 20)
#' }
plot_res <- function(r, max_res_mult) {
  r0 <- r
  r0[!is.na(values(r0))] <- 1
  r[] <- 0
  raster::res(r) <- raster::res(r) * max_res_mult
  raster::plot(r0)
  r <- raster::rasterToPolygons(r)
  raster::plot(r, add = TRUE)
}
