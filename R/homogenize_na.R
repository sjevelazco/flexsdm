#' Title
#'
#' @param x A stack or brick object.
#'
#' @return
#' @export
#'
#' @examples
homogenize_na <- function(x) {
  if (raster::canProcessInMemory(x, n = 2))
  {
    val <- raster::getValues(x)
    NA.pos <- unique(which(is.na(val), arr.ind = T)[, 1])
    val[NA.pos,] <- NA
    x <- raster::setValues(x, val)
    return(x)
  } else
  {
    x <- raster::mask(x, calc(x, fun = sum))
    return(x)
  }
}

