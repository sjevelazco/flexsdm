#' Truncate suitability predictions based on an extrapolation value
#'
#' @description Exclusion of suitability predictions in environmental conditions assumed as with high extrapolation value (EXPERIMENTAL)
#'
#' @param suit SpatRaster with suitability values
#' @param extra SpatRaster with extrapolation values preferable measured with extra_eval function
#' @param threshold numeric. Vector with one or more extrapolation values used for truncate suitability Default 50\%
#' @param trunc_value numeric. Numeric value to be used to those cells assumed to be extrapolative
#'
#' @returns
#' A SpatRaster object with truncated suitability values
#'
#' @details This function truncates suitability predictions assigning a given value, generally 0 or NA. Usage trunc_value = NA. Default 0.
#' to those cells assumed to be extrapolative (i.e., higher than a given threshold) based on an extrapolation metric like SHAPE methods (calculated with extra_eval).
#' We recommend using the function p_pdp and p_extra.
#'
#' @seealso \code{\link{extra_eval}}, \code{\link{p_extra}}, \code{\link{p_pdp}}
#'
#' @export
#'
#' @importFrom terra rast
#'
#' @examples
#' \dontrun{
#' # see examples in extra_eval function
#' }
extra_truncate <- function(suit, extra, threshold = 50, trunc_value = 0) {
  # names(suit) <- "suit"
  l <- as.list(threshold)
  for (i in 1:length(threshold)) {
    l[[i]] <- suit
    for (ii in 1:terra::nlyr(l[[i]])) {
      l[[i]][[ii]][extra > threshold[i]] <- trunc_value
    }
  }
  names(l) <- threshold
  return(l)
}
