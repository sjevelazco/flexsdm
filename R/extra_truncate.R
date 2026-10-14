#' Truncate suitability predictions based on an extrapolation value
#'
#' @description Exclusion of suitability predictions in environmental conditions assumed to be
#' extrapolative.
#'
#'
#' @param suit SpatRaster with suitability values
#' @param extra SpatRaster with extrapolation values preferable measured with extra_eval function
#' @param threshold numeric. Vector with one or more extrapolation values used for truncate suitability Default 50\%
#' @param trunc_value numeric. Numeric value to be used to those cells assumed to be extrapolative
#'
#' @returns
#' A SpatRaster object with truncated suitability values
#'
#' @details Exclusion of suitability predictions in environmental conditions assumed to be
#' extrapolative. In this function it is possible to use any metric measuring degree of
#' extrapolation (e.g., MESS-Multivariate Environmental Similarity Surfaces, EO-Environmental
#' Overlap, MOP-Mobility-Oriented Parity, EXDET-Extrapolation Detection, or AOA-Area of
#' Applicability). However, we recommend to use Shape approach (see \code{\link{extra_eval}},
#' and Velazco et al., 2023; \doi{10.1111/ecog.06992}).
#'
#' This function truncates suitability predictions assigning a given value, generally 0 or NA.
#'  Usage trunc_value = NA. Default 0.
#'
#' To those cells assumed to be extrapolative, i.e., higher than a given threshold of a given
#' extrapolation metric.
#'
#' See this \href{https://sjevelazco.github.io/flexsdm/articles/v06_Extrapolation_example.html}{vignette at flexsdm website}
#' for further details about Shape metric, model truncation, and tools to explore model extrapolation.
#'
#' @seealso \code{\link{extra_eval}}, \code{\link{p_extra}}, \code{\link{p_pdp}}, \code{\link{p_bpdp}}
#'
#' @export
#'
#' @importFrom terra rast
#'
#' @examples
#' # Toy data: small raster and records of one species to keep this example fast
#' data("spp")
#' somevar <- terra::rast(system.file("external/somevar.tif", package = "flexsdm"))
#' somevar <- terra::aggregate(somevar, fact = 6) # coarser resolution
#' sp3 <- spp[spp$species == "sp3", ]
#' sp3 <- sdm_extract(data = sp3, x = "x", y = "y", env_layer = somevar)
#'
#' # Extrapolation (Shape metric) in the projection area
#' extr <- extra_eval(
#'   training_data = sp3,
#'   pr_ab = "pr_ab",
#'   projection_data = somevar,
#'   metric = "mahalanobis"
#' )
#'
#' # Any suitability layer can be truncated; here, a toy layer rescaled to 0-1
#' suit <- somevar[[1]]
#' suit <- suit / terra::global(suit, "max", na.rm = TRUE)[[1]]
#' names(suit) <- "suit"
#'
#' # Suitability is set to 0 where extrapolation is higher than the threshold
#' suit_trunc <- extra_truncate(
#'   suit = suit,
#'   extra = extr,
#'   threshold = 50,
#'   trunc_value = 0
#' )
#' suit_trunc # one raster for each threshold value
#'
#' # see more examples in extra_eval function
extra_truncate <- function(suit, extra, threshold = 50, trunc_value = 0) {
  # names(suit) <- "suit"
  l <- as.list(threshold)
  for (i in 1:length(threshold)) {
    l[[i]] <- suit
    for (ii in 1:terra::nlyr(l[[i]])) {
      l[[i]][[ii]][extra[[1]] > threshold[i]] <- trunc_value
    }
  }
  names(l) <- threshold
  return(l)
}
