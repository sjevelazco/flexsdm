##%######################################################%##
#                                                          #
####                Raster interpolation                ####
#                                                          #
##%######################################################%##

inter <- function(r1, r2, y1, y2, rastername = NULL, dir_save = NULL, n_cores=1) {
  # r1: raster. Raster for the initial year
  # r2: raster. Raster object for the final year
  # y1: numeric. Initial year
  # y2: numeric. Final year
  # dir_save: character. Directory path and folder name where you want to save rasters
  require(raster)
  require(doParallel)
  require(parallel)
  
  annual <- (r1 - r2) / (y2 - y1)
  
  message('Number of cores: ', detectCores())
  message('Cores used: ', n_cores)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  rlist <- foreach (i = 1:(y2 - y1), .export = 'raster') %dopar%{
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
    foreach (i = 1:nlayers(rlist), .export = 'raster') %dopar%{
    raster::writeRaster(
      x = rlist[[i]],
      filename = file.path(dir_save, names(rlist)[i]),
      bylayer = TRUE,
      format = 'GTiff',
      overwrite = TRUE
    )
    NULL
    }
    stopCluster(cl)
    message(paste0('rasters were saved in: ', dir_save))
  }else{
    return(rlist)
  }
}
