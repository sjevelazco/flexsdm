#' Correction of suitability based on extrapolation
#'
#' @description Constraint suitability values under a given extrapolation value
#'
#' @param suit SpatRaster with suitability values
#' @param extra SpatRaster with extrapolation values measured in percentage (output from extra_eval function)
#' @param threshold numeric. vector with one or more values used for correct extrapolation. Default 50\%
#'
#' @returns
#' A SpatRaster object with corrected suitability values
#'
#' @seealso \code{\link{extra_eval}}
#' @export
#'
#' @importFrom terra rast
#'
#' @examples
#' \dontrun{
#' # see examples in extra_eval function
#' }
extra_correct <- function(suit, extra, threshold = 50) {
  names(suit) <- "suit"
  l <- as.list(threshold)
  for (i in 1:length(threshold)) {
    l[[i]] <- suit
    l[[i]][extra > threshold[i]] <- 0
  }
  names(l) <- threshold
  l <- terra::rast(l)
  return(l)
}
