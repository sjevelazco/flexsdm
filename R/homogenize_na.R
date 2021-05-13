#' Homogenize cells with NAs across all layers
#'
#' @param x A SpatRaster.
#'
#' @return a SpatRaster
#' @export
#'
#' @importFrom terra mask
#'
#' @examples
#' \dontrun{
#' #' require(terra)
#'
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' somevar2 <- homogenize_na(somevar)
#' par(mfrow = c(2, 1))
#' plot(somevar$CFP_4)
#' plot(somevar2$CFP_4)
#' par(mfrow = c(1, 1))
#'
#' # In somevar2 all layers have the same cells with NAs
#' }
homogenize_na <- function(x) {
  x2 <- sum(is.na(x)) == 0
  x <- terra::mask(x, x2, maskvalues = 0)
  rm(x2)
  return(x)
}
