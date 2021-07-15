#' Raster interpolation between two time periods
#'
#' @description This function calculates suitability values between two time periods with simple interpolation using two raster objects with suitability values.
#'
#' @param r1 SpatRaster. Raster object for the initial year
#' @param r2 SpatRaster. Raster object for the final year
#' @param y1 numeric. Initial year
#' @param y2 numeric. Final year
#' @param rastername character. Word used as prefix in raster file name. Default NULL
#' @param dir_save character. Directory path and name of the folder in which you want to
#' save the raster files. If NULL, function will return a SpatRaster object, on the contrary, it
#' will save raster in a given directory. Default NULL
#'
#' @importFrom terra rast nlyr writeRaster
#'
#' @return This function returns a SpatRaster if dir_save is used as NULL, If dirsave is used,
#' function will return anythink
#' @export
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#'
#' f <- system.file("external/suit_time_steps.tif", package = "flexsdm")
#' abma <- terra::rast(f)
#' plot(abma)
#'
#' int <- inter(
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
inter <- function(r1, r2, y1, y2, rastername = NULL, dir_save = NULL) {
  # dir_save: character. Directory path and folder name where you want to save raster

  annual <- (r1 - r2) / (y2 - y1)

  rlist <- list()
  for(i in 1:(y2 - y1)){
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
    for(i in 1:terra::nlyr(rlist)){
      terra::writeRaster(
        x = rlist[[i]],
        filename = paste0(file.path(dir_save, names(rlist)[i]), ".tif"),
        overwrite = TRUE
      )
    }
    message(paste0("rasters were saved in: ", dir_save))
  } else {
    return(rlist)
  }
}
