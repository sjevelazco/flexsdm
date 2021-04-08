#' Spatial block cross validation
#'
#' @param env_layer raster. Raster stack or brick with environmental variable. This will be used to evaluate spatial autocorrelation and environmental similarity between training and testing partition
#' @param data data.frame frame with presences records
#' @param x
#' @param y
#' @param pr_ab
#' @param min_res_mult
#' @param max_res_mult numeric. Maximum value will multiply raster resolution and will define the coarsest resolution to be tested, default 50.
#' @param num_grids numeric. Number of grid to be tested between 2x(raster resolution) and max_res_mult*(raster resolution), default 30
#' @param n_part
#'
#' @return
#' @export
#'
#' @importFrom ape Moran.I
#' @importFrom dplyr group_by pull
#' @importFrom flexclust dist2
#' @importFrom raster brick extent extract crs projection values res ncell cellFromXY
#' @importFrom sp coordinates
#' @importFrom stats sd
#'
#' @examples
block_partition_pa <- function(env_layer,
                               data,
                               x,
                               y,
                               pr_ab,
                               n_part = 2,
                               min_res_mult = 2,
                               max_res_mult = 200,
                               num_grids = 30) {
  if (n_part != 2) {
    stop("The use of n_part values other than 2 has not yet been implemented.")
  }

  # Transform data to data.frame and list
  data <- data.frame(data)
  data <- data[, c(pr_ab, x, y)]
  colnames(data) <- c("pr_ab", "x", "y")
  data <- data.frame(data)

  # Test some about 0 and 1 (TRY TO ADAPT THIS FUNCTION FOR WORKING ONLY FOR PRESENCE)
  if (any(!unique(data[, "pr_ab"]) %in% c(0, 1))) {
    stop(
      "values in pr_ab column did not match with 0 and 1:
unique list values in pr_ab column are: ",
      paste(unique(data[, "pr_ab"]), collapse = " ")
    )
  }

  # Vector with grid cell-size used
  cellSize <- seq(raster::res(env_layer[[1]])[1] * min_res_mult,
    raster::res(env_layer[[1]])[1] * max_res_mult,
    length.out = num_grids
  )

  message(
    "The following grid cell sizes will be tested:\n",
    paste(round(cellSize, 2), collapse = " | "),
    "\n"
  )

  # Mask
  message("Creating basic raster mask...\n")

  mask <- env_layer[[1]]
  if (class(mask) != "brick") {
    mask <- raster::brick(mask)
  }
  names(mask) <- "group"
  mask[!is.na(mask[, ])] <- 1

  # Extent
  e <- raster::extent(mask)

  # Start Cluster
  message("Searching for the optimal grid size...\n")

  # Extract coordinates----
  mask2 <- mask
  mask2[] <- 0

  # Eliminate any recrods wity NA
  filt <- raster::extract(env_layer, data[, c("x", "y")]) %>%
    stats::complete.cases()
  presences2 <- data[filt, ]

  # Transform the presences points in a DataFrameSpatialPoints
  sp::coordinates(presences2) <- presences2[, c("x", "y")]
  raster::crs(presences2) <- raster::projection(mask)

  #### Data partitioning using a grid approach ####

  # Create a list of grids based on different raster resolution
  grid <- list() # List of grids

  # raster resolution
  DIM <-
    matrix(0, length(cellSize), 2) # the number of rows and columns of each grid
  colnames(DIM) <- c("R", "C")

  for (i in 1:length(cellSize)) {
    mask3 <- mask2
    raster::res(mask3) <- cellSize[i]
    DIM[i, ] <- dim(mask3)[1:2]
    raster::values(mask3) <- 1 # Add values to cells /
    NAS <-
      c(raster::extract(mask3, presences2)) # Extract values to test if exist NAs
    if (any(is.na(NAS))) {
      while (any(is.na(NAS))) {
        raster::extent(mask3) <- raster::extent(mask3) + cellSize[i]
        raster::res(mask3) <- cellSize[i] # Give to cells a size
        DIM[i, ] <- dim(mask3)[1:2]
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

      raster::values(grid[[i]]) <- rep(group, length.out = raster::ncell(grid[[i]]))
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
  pp <- sapply(part[pa == 1, ], function(x) {
    length(unique(range(x)))
  })
  pp <- ifelse(pp == n_part, TRUE, FALSE)
  # Elimination of those partition that have one record in some group
  pf <- sapply(part[pa == 1, ], table)
  if (is.list(pf) == TRUE) {
    pf <- which(sapply(pf, min) <= 1)
  } else {
    pf <- which(apply(pf, 2, min) <= 1)
  }
  pp[pf] <- FALSE
  grid <- grid[pp]
  part <- data.frame(part[, pp])
  names(part) <- names(which(pp == TRUE))

  ### Remove problematic grids based on presences
  # Grids that assigned partitions less than the number of groups will be removed
  pa <- presences2@data[, 1] # Vector with presences and absences
  pp <- sapply(part[pa == 0, ], function(x) {
    length(unique(range(x)))
  })
  pp <- ifelse(pp == n_part, TRUE, FALSE)
  # Elimination of those partition that have one record in some group
  pf <- sapply(part[pa == 0, ], table)
  if (is.list(pf) == TRUE) {
    pf <- which(sapply(pf, min) <= 1)
  } else {
    pf <- which(apply(pf, 2, min) <= 1)
  }
  pp[pf] <- FALSE
  grid <- grid[pp]
  part <- data.frame(part[, pp])
  names(part) <- names(which(pp == TRUE))


  # Ncell
  ncell <- data.frame(matrix(
    0, nrow(presences2@data),
    length(grid)
  ))
  for (i in 1:length(grid)) {
    ncell[, i] <- raster::cellFromXY(grid[[i]], presences2)
  }

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

  dist <- flexclust::dist2(
    presences2@data[, c("x", "y")],
    presences2@data[, c("x", "y")]
  )
  dist <- 1 / dist
  diag(dist) <- 0
  dist[which(dist == Inf)] <- 0

  species2 <-
    cbind(presences2@data, raster::extract(env_layer, presences2))

  for (p in 1:length(grid)) {
    ncell3 <- ncell[, p]
    part3 <- c(part[, p])
    filt <- data.frame(
      nrow = 1:length(ncell3),
      ncell = ncell3,
      group = part3,
      pr_ab = presences2@data[c("pr_ab")]
    ) %>%
      dplyr::group_by(ncell, group, pr_ab) %>%
      dplyr::slice_sample(n = 1) %>%
      dplyr::pull(nrow) %>%
      sort()
    odd <- which((part3[filt] == 1))
    even <- which((part3[filt] == 2))
    dist2 <- dist[filt, filt]
    dist2[odd, odd] <- 0
    dist2[even, even] <- 0

    mins <- apply(dist2, 2, function(x) max(x, na.rm = TRUE))
    for (i in 1:length(mins)) {
      dist2[, i] <- ifelse(dist2[, i] == mins[i], mins[i], 0)
    }

    if (nrow(species2) < 3) {
      Imoran.Grid[p] <- NA
    } else {
      im <- sapply(
        species2[filt, names(env_layer)],
        function(x) {
          ape::Moran.I(x,
            dist2,
            na.rm = TRUE,
            scaled = TRUE
          )$observed
        }
      )
      Imoran.Grid[p] <- mean(im)
    }
  }

  Imoran.Grid <-
    abs(Imoran.Grid)
  N.grid <- 1:length(cellSize[pp])

  Opt <-
    data.frame(N.grid, cellSize = cellSize[pp], round(
      data.frame(Imoran.Grid, EnvirDist.Grid, Sd.Grid.P, Sd.Grid.A),
      3
    ))
  # Cleaning those variances based in data divided in a number of partition less than
  # the number of groups

  # SELLECTION OF THE BEST CELL SIZE----
  Opt2 <- Opt
  Dup <-
    !duplicated(Opt2[c("Imoran.Grid", "EnvirDist.Grid", "Sd.Grid.P", "Sd.Grid.A")])
  Opt2 <- Opt2[Dup, ]

  while (nrow(Opt2) > 1) {
    # I MORAN
    if (nrow(Opt2) == 1) {
      break
    }
    Opt2 <-
      Opt2[which(Opt2$Imoran.Grid <= summary(Opt2$Imoran.Grid)[2]), ]
    if (nrow(Opt2) == 1) {
      break
    }
    # Euclidean
    Opt2 <-
      Opt2[which(Opt2$EnvirDist.Grid >= summary(Opt2$EnvirDist.Grid)[5]), ]
    if (nrow(Opt2) == 1) {
      break
    }
    # SD
    Opt2 <-
      Opt2[which(Opt2$Sd.Grid.P <= summary(Opt2$Sd.Grid.P)[2]), ]
    if (nrow(Opt2) == 2) {
      break
    }
    # SD
    Opt2 <-
      Opt2[which(Opt2$Sd.Grid.A <= summary(Opt2$Sd.Grid.A)[2]), ]
    if (nrow(Opt2) == 2) {
      break
    }

    if (unique(Opt2$Imoran.Grid) &&
      unique(Opt2$EnvirDist.Grid) && unique(Opt2$Sd.Grid.P)) {
      Opt2 <- Opt2[nrow(Opt2), ]
    }
  }

  if (nrow(Opt2) > 1) {
    Opt2 <- Opt2[nrow(Opt2), ]
  }

  # Optimum size for presences
  Optimum.Grid <- grid[[Opt2$N.grid]]

  # Final data.frame result----
  result <- data.frame(presences2@data, partition = c(part[, Opt2$N.grid]))
  colnames(result) <- c("pr_ab", "x", "y", ".part")
  result <- result[c("x", "y", "pr_ab", ".part")]

  # Final data.frame result2----
  out <- list(
    ResultList = result,
    BestGridList = Opt2
  )
  return(out)
}
