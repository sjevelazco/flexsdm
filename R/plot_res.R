#' Plot different resolutions to be used in part_sblock
#'
#' @description This function is useful to display the maximum and minimum resolution that you want to test with the block_partition function. Note that if the resolution to be tested is very fine, the plot display may take a long time.
#'
#' @param r SpatRaster. A raster layer, preferably a layer of environmental variables to be used
#' @param res_mult numeric. Maximum or minimum resolution to be tested.
#'
#' @return A plot with the original raster overlapped by a grid with the resolution used
#' @export
#'
#' @importFrom terra res plot as.polygons
#'
#' @examples
#' \dontrun{
#' f <- system.file("external/somevar.tif", package = "flexsdm")
#' r <- terra::rast(f)
#' r <- r$CFP_1
#' plot_res(r, res_mult = 100)
#' plot_res(r, res_mult = 200)
#' }
plot_res <- function(r, res_mult) {
  r0 <- r
  r0[!is.na(r0)] <- 1
  r[] <- 0
  terra::res(r) <- terra::res(r) * res_mult
  terra::plot(r0)
  r <- terra::as.polygons(r)
  terra::plot(r, add = TRUE)
}
