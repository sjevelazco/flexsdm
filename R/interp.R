#' Raster interpolation between two time periods
#'
#' @description This function interpolates values for each year between two specified years with simple interpolation using two raster objects containing e.g. habitat suitability values predicted using a species distribution model.
#'
#' @param r1 SpatRaster. Raster object for the initial year
#' @param r2 SpatRaster. Raster object for the final year
#' @param y1 numeric. Initial year
#' @param y2 numeric. Final year
#' @param rastername character. Word used as prefix in raster file name. Default NULL
#' @param dir_save character. Directory path and name of the folder in which
#' the raster files will be saved. If NULL, function will return a SpatRaster object, else, it
#' will save raster in a given directory. Default NULL
#'
#' @importFrom terra rast nlyr writeRaster tmpFiles
#'
#' @return If dir_save is NULL, the function returns a SpatRaster with suitability interpolation for each year.
#' If dir_save is used, function outputs are saved in the directory
#' specified in dir_save.
#'
#' @details This function interpolates suitability values assuming that annual changes in
#' suitability are linear. This function could be useful for linking SDM output based on averaged
#' climate data and climate change scenarios to other models that require suitability values
#' disaggregated in time periods, such as population dynamics
#' (Keith et al., 2008; Conlisk et al., 2013; Syphard et al., 2013).
#'
#' @references
#' \itemize{
#' \item Keith, D.A., Akçakaya, H.R., Thuiller, W., Midgley, G.F., Pearson, R.G., Phillips, S.J.,
#'  Regan, H.M., Araujo, M.B. & Rebelo, T.G. (2008) Predicting extinction risks under climate
#'  change: coupling stochastic population models with dynamic bioclimatic habitat models.
#'  Biology Letters, 4, 560-563.
#'
#' \item Conlisk, E., Syphard, A.D., Franklin, J., Flint, L., Flint, A. & Regan, H.M. (2013)
#'  Management implications of uncertainty in assessing impacts of multiple landscape-scale
#'  threats to species persistence using a linked modeling approach. Global Change Biology
#'   3, 858-869.
#'
#' \item Syphard, A.D., Regan, H.M., Franklin, J. & Swab, R. (2013) Does functional type
#' vulnerability to multiple threats depend on spatial context in Mediterranean-climate regions?
#' Diversity and Distributions, 19, 1263-1274.
#'
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#'
#' f <- system.file("external/suit_time_step.tif", package = "flexsdm")
#' abma <- terra::rast(f)
#' plot(abma)
#'
#' int <- interp(
#'   r1 = abma[[1]],
#'   r2 = abma[[2]],
#'   y1 = 2010,
#'   y2 = 2020,
#'   rastername = "Abies",
#'   dir_save = NULL
#' )
#'
#' int
#' }
interp <- function(r1, r2, y1, y2, rastername = NULL, dir_save = NULL) {
  annual <- (r1 - r2) / (y2 - y1)

  rlist <- list()
  for (i in 1:(y2 - y1)) {
    rlist[[i]] <- (r1 - (annual * (i - 1)))
  }
  i <- length(rlist)
  rlist[[i + 1]] <- (r1 - (annual * (i)))

  if (is.null(rastername)) {
    rastername <- (y1:y2)
  } else {
    rastername <- paste(rastername, (y1:y2), sep = "_")
  }
  rlist <- terra::rast(rlist)
  names(rlist) <- rastername
  if (!is.null(dir_save)) {
    message("saving raster...")
    for (i in 1:terra::nlyr(rlist)) {
      terra::writeRaster(
        x = rlist,
        filename = paste0(file.path(dir_save, names(rlist)), ".tif"),
        overwrite = TRUE
      )
    }
    message(paste0("rasters were saved in: ", dir_save))
  } else {
    return(rlist)
  }
  terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
}
