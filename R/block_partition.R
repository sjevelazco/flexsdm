## Written by Santiago Velazco

#' Title
#'
#' @param env_layer 
#' @param occ_data 
#' @param sp 
#' @param x 
#' @param y 
#' @param pr_ab 
#' @param dir_save 
#' @param cores 
#' @param max_res_mult 
#' @param num_grids 
#' @param n_part 
#' @param save_part_raster 
#'
#' @return
#' @export
#'
#' @examples
block_partition_pa <- function(env_layer,
                               occ_data,
                               sp,
                               x,
                               y,
                               pr_ab,
                               max_res_mult=200,
                               num_grids=30,
                               n_part = 2, 
                               cores = 1,
                               save_part_raster=FALSE,
                               dir_save) {
  
  require(foreach)
  require(doParallel)
  require(parallel)
  require(flexclust)
  require(raster)
  require(dplyr)
  require(ape)
  
  # max_res_mult:numeric. Maximum value will multiply raster resolution and will define the coarsest resolution to be tested, default 50.
  # num_grids: numeric. Number of grid to be tested between 2x(raster resolution) and max_res_mult*(raster resolution), default 30
  
  # occ_data: matrix or data frame with presences records
  # n_part: 2 (dafault). integer  Number of group for data  partitioning
  # env_layer: raster. Raster stack or brick with environmental variable. This will be used to evaluate spatial autocorrelation and environmental similarity between training and testing partition
  
  dir.create(file.path(dir_save, 'block'))
  
  # Transform occ_data to data.frame and list
  occ_data <- occ_data[,c(sp, pr_ab, x, y)]
  colnames(occ_data) <- c('sp', 'pr_ab', 'x', 'y')
  occ_data <- data.frame(occ_data)
  occ_data <- split(occ_data[,-1], occ_data[,'sp'])
  
  #Vector with grid cell-size used 
  cellSize = seq(res(env_layer[[1]])[1] * 2, 
                 res(env_layer[[1]])[1] * max_res_mult, 
                 length.out = num_grids)
  
  message("The following grid cell sizes will be tested:\n", 
          paste(round(cellSize, 2), collapse = " | "), "\n")
  
  # Mask
  message("Creating basic raster mask...\n")
  
  mask <- env_layer[[1]]
  if (class(mask) != "brick") {
    mask <- raster::brick(mask)
  }
  names(mask) <- "group"
  mask[!is.na(mask[, ])] <- 1
  
  #Extent
  e <- raster::extent(mask)
  
  # Loop for each species-----
  SpNames <- names(occ_data)
  
  #Start Cluster
  if (Sys.getenv("RSTUDIO") == "1" &&
      !nzchar(Sys.getenv("RSTUDIO_TERM")) &&
      Sys.info()["sysname"] == "Darwin" &&
      as.numeric(gsub('[.]', '', getRversion())) >= 360) {
    cl <-
      parallel::makeCluster(cores, outfile = "", setup_strategy = "sequential")
  } else{
    cl <- parallel::makeCluster(cores, outfile = "")
  }
  doParallel::registerDoParallel(cl)
  
  # LOOP----
  
  message("Searching for the optimal grid size...\n")
  
  results <-
    foreach(
      s = 1:length(occ_data),
      .packages = c("raster", "ape", "dismo", "flexclust", "sp")
    ) %dopar% {
      
      # Extract coordinates----
      mask2 <- mask
      mask2[] <- 0
      
      # Eliminate any recrods wity NA
      presences2 <- occ_data[[s]]
      filt <- raster::extract(env_layer, presences2[,c(x, y)])
      filt <- complete.cases(filt)
      presences2 <- presences2[filt,]
      presences2 <- occ_data[[s]]
      
      # Transform the presences points in a DataFrameSpatialPoints
      sp::coordinates(presences2) = presences2[, c("x", "y")]
      raster::crs(presences2) <- raster::projection(mask)
      
      #### Data partitioning using a grid approach ####
      
      # Create a list of grids based on different raster resolution
      grid <- list() #List of grids
      
      # raster resolution
      DIM <-
        matrix(0, length(cellSize), 2) # the number of rows and columns of each grid
      colnames(DIM) <- c("R", "C")
      
      for (i in 1:length(cellSize)) {
        mask3 <- mask2
        res(mask3) <- cellSize[i]
        DIM[i,] <- dim(mask3)[1:2]
        raster::values(mask3) <- 1 # Add values to cells /
        NAS <-
          c(raster::extract(mask3, presences2)) # Extract values to test if exist NAs
        if (any(is.na(NAS))) {
          while (any(is.na(NAS))) {
            raster::extent(mask3) <- raster::extent(mask3) + cellSize[i]
            raster::res(mask3) <- cellSize[i] # Give to cells a size
            DIM[i,] <- dim(mask3)[1:2]
            raster::values(mask3) <- 1
            NAS <- raster::extract(mask3, presences2)
          }
        }
        grid[[i]] <- mask3
      }
      rm(mask3)
      rm(mask2)
      
      # In this section is assigned the group of each cell
      for (i in 1:length(grid)) {
        if (any(n_part == c(2, 4, 6, 8, 10))) {
          # odds number of partition
          
          group <- c(
            rep(1:n_part, DIM[i, 2])[1:DIM[i, 2]],
            rep(c((n_part / 2 + 1):n_part, 1:(n_part / 2)), DIM[i, 2])[1:DIM[i, 2]]
            )
          
          values(grid[[i]]) <- rep(group, length.out=raster::ncell(grid[[i]]))
        }
      }
      
      # Matrix within each columns represent the partitions of points
      # for each grid resolution
      part <- data.frame(matrix(0, nrow(presences2@data), length(grid)))
      for (i in 1:length(grid)) {
        part[, i] <- raster::extract(grid[[i]], presences2)
      }
      
      ### Remove problematic grids based on presences
      # Grids that assigned partitions less than the number of groups will be removed 
      pa <- presences2@data[, 1] # Vector with presences and absences
      pp <- sapply(part[pa==1, ], function(x)
        length(unique(range(x))))
      pp <- ifelse(pp == n_part, TRUE, FALSE)
      # Elimination of those partition that have one record in some group
      pf <- sapply(part[pa==1, ], table)
      if (is.list(pf) == TRUE) {
        pf <- which(sapply(pf, min) <= 1)
      } else{
        pf <- which(apply(pf, 2, min) <= 1)
      }
      pp[pf] <- FALSE
      grid <- grid[pp]
      part <- data.frame(part[, pp])
      names(part) <- names(which(pp == TRUE))
      
      ### Remove problematic grids based on presences
      # Grids that assigned partitions less than the number of groups will be removed 
      pa <- presences2@data[, 1] # Vector with presences and absences
      pp <- sapply(part[pa==0, ], function(x)
        length(unique(range(x))))
      pp <- ifelse(pp == n_part, TRUE, FALSE)
      # Elimination of those partition that have one record in some group
      pf <- sapply(part[pa==0, ], table)
      if (is.list(pf) == TRUE) {
        pf <- which(sapply(pf, min) <= 1)
      } else{
        pf <- which(apply(pf, 2, min) <= 1)
      }
      pp[pf] <- FALSE
      grid <- grid[pp]
      part <- data.frame(part[, pp])
      names(part) <- names(which(pp == TRUE))
      
      
      
      # Performance of cells ----
      # SD of number of records per cell size-----
      
      Sd.Grid.P <- rep(NA, length(grid))
      Sd.Grid.A <- rep(NA, length(grid))
      
      for (i in 1:ncol(part)) {
        Sd.Grid.A[i] <- stats::sd(table(part[pa == 0, i])) /
          mean(table(part[pa == 0, i]))
        Sd.Grid.P[i] <- stats::sd(table(part[pa == 1, i])) /
          mean(table(part[pa == 1, i]))
      }
      
      # Environmental similarity between train and test based on euclidean  -----
      EnvirDist.Grid <- rep(NA, length(grid))
      Env.P <- raster::extract(env_layer, presences2)
      for (i in 1:ncol(part)) {
        Env.P1 <- cbind(part[i], Env.P)
        Env.P2 <- split(Env.P1[, -1], Env.P1[, 1])
        euq1 <- flexclust::dist2(Env.P2[[1]], Env.P2[[2]])
        EnvirDist.Grid[i] <- mean(euq1)
        rm(Env.P1)
        rm(Env.P2)
      }
      
      # I moran-----
      Imoran.Grid <- rep(NA, length(grid))
      species2 <-
        cbind(presences2@data, raster::extract(env_layer, presences2))
      for (p in 1:length(grid)) {
        part3 <- part[,p]
          odd <- which((part3 == 1))
          even <- which((part3 == 2))
          dist <- as.matrix(dist(presences2@data[,c('x', 'y')]))
          dist <- 1 / dist
          diag(dist) <- 0
          dist[which(dist == Inf)] <- 0
          dist[odd, odd] <- 0
          dist[even, even] <- 0
          mins <- apply(dist, 2, max)
          for (i in 1:length(mins)) {
            dist[, i] <- ifelse(dist[, i] == mins[i], mins[i], 0)
          }
        
        if (nrow(species2) < 3) {
          Imoran.Grid[p] <- NA
        } else{
          im <- sapply(species2[, names(env_layer)],
                       function(x)
                         ape::Moran.I(x,
                                      dist,
                                      na.rm = T,
                                      scaled = T)$observed)
          Imoran.Grid[p] <- mean(im)
        }
      }
      
      Imoran.Grid <-
        abs(Imoran.Grid)
      N.grid <- 1:length(cellSize[pp])
      
      Opt <-
        data.frame(N.grid, cellSize=cellSize[pp], round(data.frame(
          Imoran.Grid, EnvirDist.Grid, Sd.Grid.P, Sd.Grid.A
        ), 3))
      # Cleaning those variances based in data divided in a number of partition less than
      # the number of groups
      
      # SELLECTION OF THE BEST CELL SIZE----
      Opt2 <- Opt
      Dup <-
        !duplicated(Opt2[c("Imoran.Grid", "EnvirDist.Grid", "Sd.Grid.P", "Sd.Grid.A")])
      Opt2 <- Opt2[Dup, ]
      
      while (nrow(Opt2) > 1) {
        # I MORAN
        if (nrow(Opt2) == 1)
          break
        Opt2 <-
          Opt2[which(Opt2$Imoran.Grid <= summary(Opt2$Imoran.Grid)[2]), ]
        if (nrow(Opt2) == 1)
          break
        # Euclidean
        Opt2 <-
          Opt2[which(Opt2$EnvirDist.Grid >= summary(Opt2$EnvirDist.Grid)[5]), ]
        if (nrow(Opt2) == 1)
          break
        # SD
        Opt2 <-
          Opt2[which(Opt2$Sd.Grid.P <= summary(Opt2$Sd.Grid.P)[2]), ]
        if (nrow(Opt2) == 2)
          break
        # SD
        Opt2 <-
          Opt2[which(Opt2$Sd.Grid.A <= summary(Opt2$Sd.Grid.A)[2]), ]
        if (nrow(Opt2) == 2)
          break
        
        if (unique(Opt2$Imoran.Grid) &&
            unique(Opt2$EnvirDist.Grid) && unique(Opt2$Sd.Grid.P)) {
          Opt2 <- Opt2[nrow(Opt2), ]
        }
      }
      
      if (nrow(Opt2) > 1) {
        Opt2 <- Opt2[nrow(Opt2), ]
      }
      
      # Optimum size for presences
      print(Opt2)
      Optimum.Grid <- grid[[Opt2$N.grid]]
      
      # Final data.frame result----
      result <- data.frame(SpNames[s], presences2@data, partition = c(part[, Opt2$N.grid]))
      colnames(result) <- c("sp", "pr_ab", "x", "y", "partition")
      result <- result[c("sp", "x", "y", "pr_ab", "partition")]
      
      #Save blocks raster
      if(save_part_raster){
        dir.create(file.path(dir_save, 'block', 'block_layer'))
        dir_save_r <- file.path(dir_save, 'block', 'block_layer')
        pseudo.mask <- mask
        RtoP <- data.frame(raster::rasterToPoints(mask)[, -3])
        sp::coordinates(RtoP) = c("x", "y")
        raster::crs(RtoP) <- raster::projection(mask)
        Ncell <- raster::cellFromXY(mask, RtoP)
        RtoP <- raster::extract(Optimum.Grid, RtoP)
        pseudo.mask[Ncell] <- RtoP
        rm(RtoP)
        rm(Ncell)
        raster::writeRaster(
          pseudo.mask,
          paste(dir_save_r, paste(SpNames[s], '.tif', sep = ""), sep = '/'),
          format = 'GTiff',
          overwrite = TRUE
        )

      }
      
      
      Opt2 <- data.frame(Sp = SpNames[s], Opt2)
      
      # Final data.frame result2----
      out <- list(ResultList = result,
                  BestGridList = Opt2)
      return(out)
    }
  parallel::stopCluster(cl)
  
  message('Saving results...')
  
  FinalResult <- dplyr::bind_rows(lapply(results, function(x) x[[1]]))
  FinalInfoGrid <- dplyr::bind_rows(lapply(results, function(x) x[[2]]))
  rm(results)
  utils::write.table(
    FinalResult,
    file.path(dir_save, 'block',"OccBlocks.txt"),
    sep = "\t",
    row.names = F
  )
  utils::write.table(
    FinalInfoGrid,
    file.path(dir_save, 'block', "BestPartitions.txt"),
    sep = "\t",
    col.names = T,
    row.names = F
  )
  return(FinalResult)
}
