## Written by Santiago Velazco
utils::globalVariables("s")

env_layer = r_base
occ_data = sp_db
sp='species'
x='x'
y='y'
pr_ab='pr_ab'
max_res_mult=200
num_grids=30

block_partition_pa <- function(env_layer = NULL,
                               occ_data = NULL,
                               sp,
                               x,
                               y,
                               pr_ab,
                               dir_save = NULL,
                               cores = NULL,
                               max_res_mult=100,
                               num_grids=30,
                               type = NULL,
                               pseudoabsencesMethod = NULL,
                               N = 2,
                               PrAbRatio = NULL,
                               DirM = NULL,
                               MRst = NULL,
                               Geo_Buf = NULL) {
  
  # max_res_mult:numeric. Maximum value will multiply raster resolution and will define the coarsest resolution to be tested, default 50.
  # num_grids: numeric. Number of grid to be tested between 2x(raster resolution) and max_res_mult*(raster resolution), defaul 30
  
  # occ_data: matrix or data frame with presences records
  # N: 2 (dafault). interger  Number of group for data  paritioning
  # pseudoabsences: logical, TRUE (dafault).
  # mask: Raster object. Preferible on wiht the same resolution and extent than
  #       variable used in the future
  # pseudoabsencesMethod: a character string indicating which pseudo-absences
  #                       method is to be computed
  # PrAbRatio: numeric. value of PrAbRatio to be computed.
  # env_layer: Raster object. Variable set to be used in pseusoabsences
  # cellSize: numeric vector. a vector of values with different cell grid sizes
  
  # Transform occ_data to data.frame and list
  occ_data <- occ_data[,c(sp, pr_ab, x, y)]
  colnames(occ_data) <- c('sp', 'pr_ab', 'x', 'y')
  occ_data <- data.frame(occ_data)
  occ_data <- split(occ_data[,-1], occ_data[,'sp'])
  
  #Cellsize
  cellSize = seq(res(env_layer[[1]])[1] * 2, 
                 res(env_layer[[1]])[1] * max_res_mult, 
                 length.out = num_grids)
  
  
  # Mask
  mask <- env_layer[[1]]
  if (class(mask) != "brick") {
    mask <- raster::brick(mask)
  }
  names(mask) <- "Group"
  mask[!is.na(mask[, ])] <- 1
  
  #Extent
  e <- raster::extent(mask)
  
  # Loop for each species-----
  # ResultList <- rep(list(NULL),length(occ_data))
  SpNames <- names(occ_data)
  # BestGridList <- rep(list(NULL),length(occ_data))
  
  #Start Cluster
  if (Sys.getenv("RSTUDIO") == "1" &&
      !nzchar(Sys.getenv("RSTUDIO_TERM")) &&
      Sys.info()["sysname"] == "Darwin" &&
      as.numeric(gsub('[.]', '', getRversion())) >= 360) {
    cl <- parallel::makeCluster(cores,outfile="", setup_strategy = "sequential")
  }else{
    cl <- parallel::makeCluster(cores,outfile="")
  }
  doParallel::registerDoParallel(cl)
  
  # LOOP----
  results <-
    foreach(
      s = 1:length(occ_data),
      .packages = c("raster", "ape", "dismo"),
      .export = c("inv_bio","MESS","inv_geo","KM_BLOCK","OptimRandomPoints")
    ) %dopar% {
      
      print(paste(s, SpNames[s]))
      # Extract coordinates----
      mask2 <- mask
      mask2[] <- 0
      presences2 <- occ_data[[s]]
      presences <- occ_data[[s]]
      
      # Transform the presences points in a DataFrameSpatialPoints
      sp::coordinates(presences2) = presences2[, c("x", "y")]
      raster::crs(presences2) <- raster::projection(mask)
      presences2@data <- presences2@data['pr_ab']
      
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
        # grid[[i]] <- rasterToPolygons(mask3)
      }
      rm(mask3)
      
      # In this section is assigned the group of each cell
      for (i in 1:length(grid)) {
        if (any(N == c(2, 4, 6, 8, 10))) {
          # odds number of partition
          
          group <- c(
            rep(1:N, DIM[i, 2])[1:DIM[i, 2]],
            rep(c((N / 2 + 1):N, 1:(N / 2)), DIM[i, 2])[1:DIM[i, 2]]
            )
          
          values(grid[[i]]) <- rep(group, length.out=raster::ncell(grid[[i]]))
        }
      }
      plot(grid[[5]])
      plot(grid[[10]])
      plot(grid[[20]])
      plot(grid[[30]])
      
      # Matrix within each columns represent the partitions of points
      # for each grid resolution
      part <- data.frame(matrix(0, nrow(presences2@data), length(grid)))
      for (i in 1:length(grid)) {
        part[, i] <- raster::extract(grid[[i]], presences2)
      }
      
      part2 <- list()
      for (i in 1:length(grid)) {
        part2[[i]] <- data.frame(raster::extract(grid[[i]], presences2), presences)
      }
      
      # Here will be deleted grids that assigned partitions less than the number
      # of groups
      pp <- sapply(part[1:nrow(presences), ], function(x)
        length(unique(range(x))))
      pp <- ifelse(pp == N, TRUE, FALSE)
      # Elimination of those partition that have one record in some group
      pf <- sapply(part[1:nrow(presences), ], table)
      if (is.list(pf) == TRUE) {
        pf <- which(sapply(pf, min) == 1)
      } else{
        pf <- which(apply(pf, 2, min) == 1)
      }
      pp[pf] <- FALSE
      grid <- grid[pp]
      part <- data.frame(part[, pp])
      names(part) <- names(which(pp == T))
      part2 <- part2[pp]
      
      # Performace of cells ----
      # SD of number of records per cell size-----
      pa <- presences2@data[, 1] # Vector with presences and absences
      Sd.Grid.P <- rep(NA, length(grid))
      Sd.Grid.A <- rep(NA, length(grid))
      
      for (i in 1:ncol(part)) {
        Sd.Grid.A[i] <- stats::sd(table(part[pa == 0, i])) /
          mean(table(part[pa == 0, i]))
      }
      for (i in 1:ncol(part)) {
        Sd.Grid.P[i] <- stats::sd(table(part[pa == 1, i])) /
          mean(table(part[pa == 1, i]))
      }
      
      # MESS -----
      Mess.Grid.P <- rep(NA, length(grid))
      Env.P <- raster::extract(env_layer, presences)
      for (i in 1:ncol(part)) {
        Env.P1 <- cbind(part[i], Env.P)
        Env.P2 <- split(Env.P1[, -1], Env.P1[, 1])
        mess1 <- MESS(Env.P2[[1]], Env.P2[[2]])
        Mess.Grid.P[i] <- mean(mess1$TOTAL, na.rm = TRUE)
        rm(Env.P1)
      }
      
      # Imoran-----
      pc1 <- RStoolbox::rasterPCA(env_layer, spca = T,nComp = 1)$map
      Imoran.Grid.P <- rep(NA, length(grid))
      for (p in 1:length(grid)) {
        part3 <- part2[[p]]
        # rownames(part3) <-
        #   paste(part3$group, part3$C, part3$R, part3$lon, part3$lat)
        if (type == "nearest") {
          mineucli <- list()
          for (r in 1:DIM[p, 'R']) {
            part4 <- part3[part3$R == r,]
            if (nrow(part4) >= 2) {
              for (c in 1:DIM[p, 'C']) {
                A <- part4[part4$C == c, ]
                B <- part4[part4$C == (c + 1), ]
                euclilist <- matrix(0, (nrow(A) * nrow(B)), 3)
                if ((nrow(A) >= 1 & nrow(B) >= 1)) {
                  coord <- data.frame(
                    expand.grid(rownames(A), rownames(B)),
                    cbind(expand.grid(A[, 5], B[, 5]),
                          expand.grid(A[, 4], B[, 4]))
                  )
                  colnames(coord) <- c("A", "B", "xa", "xb", "ya", "yb")
                  coord$eucli <-
                    sqrt((coord$xa - coord$xb) ^ 2 + (coord$ya - coord$yb) ^ 2)
                  mineucli[[length(mineucli) + 1]] <-
                    coord[which.min(coord$eucli),]
                }
              }
            }
          }
          for (r in 1:DIM[p, 'C']) {
            part4 <- part3[part3$C == r,]
            if (nrow(part4) >= 2) {
              for (c in 1:DIM[p, 'R']) {
                A <- part4[part4$R == c, ]
                B <- part4[part4$R == (c + 1), ]
                if ((nrow(A) >= 1 & nrow(B) >= 1)) {
                  coord <- data.frame(
                    expand.grid(rownames(A), rownames(B)),
                    cbind(expand.grid(A[, 5], B[, 5]),
                          expand.grid(A[, 4], B[, 4]))
                  )
                  colnames(coord) <- c("A", "B", "xa", "xb", "ya", "yb")
                  coord$eucli <-
                    sqrt((coord$xa - coord$xb) ^ 2 + (coord$ya - coord$yb) ^ 2)
                  mineucli[[length(mineucli) + 1]] <-
                    coord[which.min(coord$eucli),]
                }
              }
            }
          }
          
          euclida <- plyr::ldply(mineucli, data.frame)
          euclida$A <- as.character(euclida$A)
          euclida$B <- as.character(euclida$B)
          species2 <- presences[c(euclida$A, euclida$B), ]
          dist <- as.matrix(dist(species2))
          dist <- 1 / dist
          diag(dist) <- 0
          dist[which(dist == Inf)] <- 0
          species2$pc1 <- raster::extract(env_layer[[1]], species2)
        }
        
        if (type == "all") {
          odd <- which((part3$group == 1))
          even <- which((part3$group == 2))
          dist <- as.matrix(dist(presences))
          dist <- 1 / dist
          diag(dist) <- 0
          dist[which(dist == Inf)] <- 0
          dist[odd, odd] <- 0
          dist[even, even] <- 0
          mins <- apply(dist, 2, max)
          for (i in 1:length(mins)) {
            dist[, i] <- ifelse(dist[, i] == mins[i], mins[i], 0)
          }
          species2 <-
            cbind(presences, pc1 = raster::extract(pc1, presences))
        }
        
        if (nrow(species2) < 3) {
          Imoran.Grid.P[p] <- NA
        } else{
          Imoran.Grid.P[p] <-
            Moran.I(species2$pc1,
                    dist,
                    na.rm = T,
                    scaled = T)$observed
        }
      }
      
      Imoran.Grid.P <-
        abs(Imoran.Grid.P) # OJO estamos dejando todos los valores positivos
      N.grid <- 1:length(cellSize[pp])
      
      Opt <-
        data.frame(N.grid, cellSize[pp], round(data.frame(
          Imoran.Grid.P, Mess.Grid.P, Sd.Grid.P
        ), 3))
      # Cleaning those variances based in data divided in a number of partition less than
      # the number of groups
      
      # SELLECTION OF THE BEST CELL SIZE----
      Opt2 <- Opt
      Dup <-
        !duplicated(Opt2[c("Imoran.Grid.P", "Mess.Grid.P", "Sd.Grid.P")])
      Opt2 <- Opt2[Dup, ]
      
      while (nrow(Opt2) > 1) {
        # I MORAN
        if (nrow(Opt2) == 1)
          break
        Opt2 <-
          Opt2[which(Opt2$Imoran.Grid.P <= summary(Opt2$Imoran.Grid.P)[2]), ]
        if (nrow(Opt2) == 1)
          break
        # MESS
        Opt2 <-
          Opt2[which(Opt2$Mess.Grid.P >= summary(Opt2$Mess.Grid.P)[5]), ]
        if (nrow(Opt2) == 1)
          break
        # SD
        Opt2 <-
          Opt2[which(Opt2$Sd.Grid.P <= summary(Opt2$Sd.Grid.P)[2]), ]
        if (nrow(Opt2) == 2)
          break
        
        if (unique(Opt2$Imoran.Grid.P) &&
            unique(Opt2$Mess.Grid.P) && unique(Opt2$Sd.Grid.P)) {
          Opt2 <- Opt2[nrow(Opt2), ]
        }
      }
      
      if (nrow(Opt2) > 1) {
        Opt2 <- Opt2[nrow(Opt2), ]
      }
      
      # Optimum size for presences
      print(Opt2)
      Optimum.Grid <- grid[[Opt2$N.grid]]
      Optimum.Grid@data[, c("C", "R")] <- NULL
      presences <- data.frame(Partition = part[, Opt2$N.grid], presences)
      
      #Save blocks raster
      pseudo.mask <- mask
      pseudo.mask2 <- list()
      RtoP <- data.frame(raster::rasterToPoints(mask)[, -3])
      sp::coordinates(RtoP) = c("x", "y")
      raster::crs(RtoP) <- raster::projection(mask)
      FILTER <- sp::over(RtoP, Optimum.Grid)
      pseudo.mask[which(pseudo.mask[] == 1)] <- as.matrix(FILTER)
      # writeRaster(pseudo.mask, paste(dir_save, paste(SpNames[s],'.tif',sep=""),sep='/'),
      #             format = 'GTiff', NAflag = -9999, overwrite = TRUE)
      #
      for (i in 1:N) {
        mask3 <- pseudo.mask
        mask3[!mask3[] == i] <- 0
        pseudo.mask2[[i]] <- mask3
      }
      #
      pseudo.mask <- raster::brick(pseudo.mask2)
      # rm(pseudo.mask2)
      
      ##%######################################################%##
      #                                                          #
      ####             Pseudoabsences allocation              ####
      #                                                          #
      ##%######################################################%##
      
      
      # Pseudo-Absences with Random allocation-----
      if (pseudoabsencesMethod == "RND") {
        pseudo.mask_p <- pseudo.mask
        pseudo.mask <- sum(pseudo.mask_p)
        pseudo.mask_p[pseudo.mask_p==0] <- NA
        
        raster::writeRaster(
          pseudo.mask,
          paste(dir_save, paste(SpNames[s], '.tif', sep = ""), sep = '/'),
          format = 'GTiff',
          NAflag = -9999,
          overwrite = TRUE
        )
        
        # Random allocation of Pseudo-Absences
        absences <- list()
        for (i in 1:N) {
          set.seed(s)
          if (!is.null(MRst)) {
            SpMask <- raster::raster(file.path(DirM, paste0(SpNames[s], ".tif")))
            pseudo.mask_p[[i]] <- pseudo.mask_p[[i]] * SpMask
            # if(sum(is.na(SpMask[])==F)<(PrAbRatio*nrow(occ_data[[s]]))){
            # warning("The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences")
            # stop("Please try again with another restriction type or without restricting the extent")
            # }
          }
          # absences.0 <-
          #   dismo::randomPoints(
          #     pseudo.mask_p[[i]],
          #     (1 / PrAbRatio) * sum(presences[, 1] == i),
          #     ext = e,
          #     prob = FALSE
          #   )
          
          absences.0 <- OptimRandomPoints(r=pseudo.mask_p[[i]], n=(1 / PrAbRatio)*sum((presences[, 1] == i)),p=presences[presences[, 1] == i, 2:3] )
          colnames(absences.0) <- c("lon", "lat")
          absences[[i]] <- as.data.frame(absences.0)
        }
        names(absences) <- 1:N
      }
      
      # Pseudo-Absences allocation with Environmental constrain ----
      if (pseudoabsencesMethod == "ENV_CONST") {
        pseudo.mask_p <- inv_bio(env_layer, presences[, -1])
        
        # Split the raster of environmental layer with grids
        pseudo.mask_p <- raster::mask(pseudo.mask, pseudo.mask_p)
        pseudo.mask <- sum(pseudo.mask_p)
        pseudo.mask_p[pseudo.mask_p==0] <- NA
        
        raster::writeRaster(
          pseudo.mask,
          paste(dir_save, paste(SpNames[s], '.tif', sep = ""), sep = '/'),
          format = 'GTiff',
          NAflag = -9999,
          overwrite = TRUE
        )
        
        absences <- list()
        for (i in 1:N) {
          set.seed(s)
          if (!is.null(MRst)) {
            SpMask <-
              raster::raster(file.path(DirM, paste0(SpNames[s], ".tif")))
            pseudo.mask_p[[i]] <- pseudo.mask_p[[i]] * SpMask
            if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(occ_data[[s]]))) {
              warning(
                "The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences"
              )
              stop(
                "Please try again with another restriction type or without restricting the extent"
              )
            }
          }
          # absences.0 <-
          #   dismo::randomPoints(
          #     pseudo.mask_p[[i]],
          #     (1 / PrAbRatio) * sum(presences[, 1] == i),
          #     ext = e,
          #     prob = FALSE
          #   )
          absences.0 <- OptimRandomPoints(r=pseudo.mask_p[[i]], n=(1 / PrAbRatio)*sum((presences[, 1] == i)),p=presences[presences[, 1] == i, 2:3] )
          
          colnames(absences.0) <- c("lon", "lat")
          absences[[i]] <- as.data.frame(absences.0)
        }
        
        names(absences) <- 1:N
      }
      
      # Pseudo-Absences allocation with Geographical constrain-----
      if (pseudoabsencesMethod == "GEO_CONST") {
        pseudo.mask_p <-
          inv_geo(e = env_layer, p = presences[, -1], d = Geo_Buf)
        
        # Split the raster of environmental layer with grids
        pseudo.mask_p <- raster::mask(pseudo.mask, pseudo.mask_p)
        pseudo.mask <- sum(pseudo.mask_p)
        pseudo.mask_p[pseudo.mask_p==0] <- NA
        
        raster::writeRaster(
          pseudo.mask,
          paste(dir_save, paste(SpNames[s], '.tif', sep = ""), sep = '/'),
          format = 'GTiff',
          NAflag = -9999,
          overwrite = TRUE
        )
        
        absences <- list()
        for (i in 1:N) {
          set.seed(s)
          if (!is.null(MRst)) {
            SpMask <-
              raster::raster(file.path(DirM, paste0(SpNames[s], ".tif")))
            pseudo.mask_p[[i]] <- pseudo.mask_p[[i]] * SpMask
            if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(occ_data[[s]]))) {
              warning(
                "The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences"
              )
              stop(
                "Please try again with a smaller geographical buffer or without restricting the accessible area"
              )
            }
          }
          # absences.0 <-
          #   dismo::randomPoints(
          #     pseudo.mask_p[[i]],
          #     (1 / PrAbRatio) * sum(presences[, 1] == i),
          #     ext = e,
          #     prob = FALSE
          #   )
          absences.0 <- OptimRandomPoints(r=pseudo.mask_p[[i]], n=(1 / PrAbRatio)*sum((presences[, 1] == i)),p=presences[presences[, 1] == i, 2:3] )
          
          colnames(absences.0) <- c("lon", "lat")
          absences[[i]] <- as.data.frame(absences.0)
        }
        
        names(absences) <- 1:N
      }
      
      # Pseudo-Absences allocation with Environmentla and Geographical  constrain-----
      if (pseudoabsencesMethod == "GEO_ENV_CONST") {
        pseudo.mask_p <- inv_bio(env_layer, presences[, -1])
        pseudo.mask_pg <-
          inv_geo(e = env_layer, p = presences[, -1], d = Geo_Buf)
        pseudo.mask_p <- pseudo.mask_p * pseudo.mask_pg
        
        # Split the raster of environmental layer with grids
        pseudo.mask_p <- raster::mask(pseudo.mask, pseudo.mask_p)
        pseudo.mask <- sum(pseudo.mask_p)
        pseudo.mask_p[pseudo.mask_p==0] <- NA
        
        raster::writeRaster(
          pseudo.mask,
          paste(dir_save, paste(SpNames[s], '.tif', sep = ""), sep = '/'),
          format = 'GTiff',
          NAflag = -9999,
          overwrite = TRUE
        )
        
        absences <- list()
        for (i in 1:N) {
          set.seed(s)
          if (!is.null(MRst)) {
            SpMask <-
              raster::raster(file.path(DirM, paste0(SpNames[s], ".tif")))
            pseudo.mask_p[[i]] <- pseudo.mask_p[[i]] * SpMask
            if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(occ_data[[s]]))) {
              warning(
                "The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences"
              )
              stop(
                "Please try again with another restriction type or without restricting the extent"
              )
            }
          }
          # absences.0 <-
          #   dismo::randomPoints(pseudo.mask_p[[i]],
          #                       (1 / PrAbRatio) * sum(presences[, 1] == i),
          #                       prob = FALSE)
          
          absences.0 <- OptimRandomPoints(r=pseudo.mask_p[[i]], n=(1 / PrAbRatio)*sum((presences[, 1] == i)),p=presences[presences[, 1] == i, 2:3] )
          colnames(absences.0) <- c("lon", "lat")
          absences[[i]] <- as.data.frame(absences.0)
        }
        
        names(absences) <- 1:N
      }
      
      # Pseudo-Absences allocation with Environmentla and Geographical and k-mean constrain-----
      if (pseudoabsencesMethod == "GEO_ENV_KM_CONST") {
        pseudo.mask_p <- inv_bio(env_layer, presences[, -1])
        pseudo.mask_pg <-
          inv_geo(e = env_layer, p = presences[, -1], d = Geo_Buf)
        pseudo.mask_p <- pseudo.mask_p * pseudo.mask_pg
        
        # Split the raster of environmental layer with grids
        pseudo.mask_p <- raster::mask(pseudo.mask, pseudo.mask_p)
        pseudo.mask <- sum(pseudo.mask_p)
        pseudo.mask_p[pseudo.mask_p==0] <- NA
        
        raster::writeRaster(
          pseudo.mask,
          paste(dir_save, paste(SpNames[s], '.tif', sep = ""), sep = '/'),
          format = 'GTiff',
          NAflag = -9999,
          overwrite = TRUE
        )
        
        absences <- list()
        for (i in 1:N) {
          set.seed(s)
          if (!is.null(MRst)) {
            SpMask <- raster::raster(file.path(DirM, paste0(SpNames[s], ".tif")))
            pseudo.mask_p[[i]] <- pseudo.mask_p[[i]] * SpMask
            if (sum(is.na(SpMask[]) == F) < (PrAbRatio * nrow(occ_data[[s]]))) {
              warning(
                "The ammount of cells in the M restriction is insuficient to generate a 1:1 number of pseudo-absences"
              )
              stop(
                "Please try again with another restriction type or without restricting the extent"
              )
            }
          }
          
          absences.0 <-
            KM_BLOCK(
              raster::rasterToPoints(pseudo.mask_p[[i]])[, -3],
              raster::mask(env_layer, pseudo.mask_p[[i]]),
              (1 / PrAbRatio) * sum(presences[, 1] == i)
            )
          colnames(absences.0) <- c("lon", "lat")
          
          absences[[i]] <- as.data.frame(absences.0)
        }
        names(absences) <- 1:N
      }
      
      absences <- plyr::ldply(absences, data.frame)
      names(absences) <- c("Partition", "x", "y")
      absences[, c("x", "y")] <- round(absences[, c("x", "y")], 4)
      colnames(absences) <- colnames(presences)
      # Final data.frame result----
      PresAbse <-
        rep(c(1, 0), sapply(list(presences, absences), nrow))
      result <-
        data.frame(
          Sp = SpNames[s],
          PresAbse,
          rbind(presences, absences),
          stringsAsFactors = F
        )
      result <- result[, c("Sp", "x", "y", "Partition", "PresAbse")]
      
      Opt2 <- data.frame(Sp = SpNames[s], Opt2)
      
      # Final data.frame result2----
      out <- list(ResultList = result,
                  BestGridList = Opt2)
      # utils::write.table(result,paste(dir_save, paste0(SpNames[s],".txt"), sep="\\"), sep="\t",row.names=F)
      return(out)
    }
  
  parallel::stopCluster(cl)
  FinalResult <- dplyr::bind_rows(lapply(results, function(x) x[[1]]))
  FinalInfoGrid <- dplyr::bind_rows(lapply(results, function(x) x[[2]]))
  
  colnames(FinalResult) <- c("sp", "x", "y", "Partition", "PresAbse")
  utils::write.table(
    FinalResult,
    paste(dir_save, "OccBlocks.txt", sep = "\\"),
    sep = "\t",
    row.names = F
  )
  utils::write.table(
    FinalInfoGrid,
    paste(dir_save, "BestPartitions.txt", sep = '/'),
    sep = "\t",
    col.names = T,
    row.names = F
  )
  
  return(FinalResult)
}