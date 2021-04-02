#' Raster interpolation
#'
#' @param r1 raster. Raster object for the initial year
#' @param r2 raster. Raster object for the final year
#' @param y1 numeric. Initial year
#' @param y2 numeric. Final year
#' @param rastername character. word used as prefix in raster file name 
#' @param dir_save character. Directory path and name of the folder in which you want to save the raster files
#' @param n_cores numeric. Number of cores use for parallelization
#'
#' @importFrom raster brick writeRaster
#' 
#' @return
#' @export
#'
#' @examples
inter <- function(r1, r2, y1, y2, rastername = NULL, dir_save = NULL, n_cores=1) {
  # dir_save: character. Directory path and folder name where you want to save raster
  
  annual <- (r1 - r2) / (y2 - y1)
  
  message('Number of cores: ', parallel::detectCores())
  message('Cores used: ', n_cores)
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  rlist <- foreach::foreach (i = 1:(y2 - y1), .export = 'raster') %dopar%{
    (r1 - (annual * (i - 1)))
  }
  i <- length(rlist)
  rlist[[i + 1]] <- (r1 - (annual * (i)))
  # unique(values(r2-rlist[[length(rlist)]])) #test
  
  if (is.null(rastername)) {
    rastername <- (y1:y2)
  } else{
    rastername <- paste(rastername, (y1:y2), sep = '_')
  }
  names(rlist) <- rastername
  rlist <- raster::brick(rlist)
  if(!is.null(dir_save)){
    message('saving raster...')
    foreach::foreach (i = 1:raster::nlayers(rlist), .export = 'raster') %dopar%{
    raster::writeRaster(
      x = rlist[[i]],
      filename = file.path(dir_save, names(rlist)[i]),
      bylayer = TRUE,
      format = 'GTiff',
      overwrite = TRUE
    )
    NULL
    }
    parallel::stopCluster(cl)
    message(paste0('rasters were saved in: ', dir_save))
  }else{
    return(rlist)
  }
}
