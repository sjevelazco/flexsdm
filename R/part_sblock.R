#' Spatial block cross validation
#'
#' @description This function explores spatial blocks with different cell sizes and returns the
#' best one suited for a given presence or presence-absences database. The selection of the best
#' grid size is performed automatically considering spatial autocorrelation, environmental
#' similarity, and the number of presence and absence records in each partition.
#'
#' @param env_layer SpatRaster. Raster with environmental
#' variable. This will be used to evaluate spatial autocorrelation and
#' environmental similarity between training and testing partition. Because this function
#' calculate dissimilarity based on euclidean distances, it can only handle continuous
#' layers, do not use categorical layers as inputs
#' @param data data.frame. Data.frame or tibble object with presences
#' (or presence-absence, or presences-pseudo-absence) records, and coordinates
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param pr_ab character. Column with presences, presence-absence,
#' or pseudo-absence. Presences must be represented by 1 and absences by 0
#' @param min_res_mult integer. Minimum value used for multiplying
#' raster resolution and define the finest resolution to be tested, default 3.
#' @param max_res_mult integer. Maximum value used for multiplying
#' raster resolution and define the coarsest resolution to be tested, default 200.
#' @param num_grids integer. Number of grid to be tested between
#' min_res_mult X (raster resolution) and max_res_mult X (raster resolution), default 30
#' @param n_part  integer. Number of partition. Default 2, values other than
#' 2 has not yet been implemented.
#' @param prop numeric. Proportion of point used for testing autocorrelation between
#' groups (values > 0 and <=1). The smaller this number is, the faster the function will work.
#' Default 0.5
#'
#' @return
#' A list with:
#' \itemize{
#'   \item part: A tibble object with information used in 'data' arguments and a additional column .part with partition group.
#'   \item best_part_info: A tibble with information of the bets partition. It contains the number of the best partition (n_grid), cell size (cell_size), standard deviation of presences (sd_p), standard deviation of absences (sd_a), Moran's I spatial autocorrelation (spa_auto), and environmental similarity based on euclidean distance (env_sim).
#'   \item grid: A SpatRaster object with blocks
#'   }
#'
#' @details write here criteria used for performing the search of the best partition (metrics and quartil selection).
#'
#' @export
#'
#' @importFrom ape Moran.I
#' @importFrom dplyr tibble pull group_by slice_sample select
#' @importFrom flexclust dist2
#' @importFrom stats complete.cases sd
#' @importFrom terra extract res ext vect crs extend values ncell cellFromXY geom
#' @importFrom utils combn
#'
#' @seealso \code{\link{part_classical}}, \code{\link{part_senv}}, \code{\link{get_block}}, and \code{\link{plot_res}}.
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#'
#' # Load datasets
#' data(spp)
#' f <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(f)
#'
#' # Lest practice with a single species
#' single_spp <- spp %>% dplyr::filter(species == "sp3")
#' part <- part_sblock(
#'   env_layer = somevar,
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   min_res_mult = 10,
#'   max_res_mult = 500,
#'   num_grids = 30,
#'   n_part = 2,
#'   prop = 0.5
#' )
#' part
#'
#' part$part
#' part$best_part_info
#' part$grid
#'
#' # Lets explore Grid object
#'
#' plot(part$grid)
#' points(part$part[c("x", "y")],
#'   col = c("blue", "red")[part$part$.part],
#'   cex = 0.5,
#'   pch = 19
#' )
#'
#' terra::res(part$grid)
#' terra::res(somevar)
#'
#' # Note that is a layer with block partition, but it has a
#' # different resolution than the original environmental variables.
#' # In the case you wish have a layer with the same properties
#' # (i.e. resolution, extent, NAs) than your original environmental
#' # variables you can use the \code{\link{get_block}} function.
#'
#' grid_env <- get_block(env_layer = somevar, best_grid = part$grid)
#'
#' plot(grid_env) # this is a block layer with the same layer
#' # properties as environmental variables.
#' points(part$part[c("x", "y")],
#'   col = c("blue", "red")[part$part$.part],
#'   cex = 0.5,
#'   pch = 19
#' )
#' # This layer could be very useful in case you need sample
#' # pseudo_absence or background point
#' # See examples in \code{\link{backgroudp}} and \code{\link{pseudoabs}}
#'
#'
#' # Lest try with a higher number of partition
#' part <- part_sblock(
#'   env_layer = somevar,
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   min_res_mult = 10,
#'   max_res_mult = 500,
#'   num_grids = 30,
#'   n_part = 4,
#'   prop = 0.5
#' )
#'
#' # Lets explore Grid object
#' plot(part$grid, col = gray.colors(4))
#' points(part$part[c("x", "y")],
#'   col = rainbow(n = 4)[part$part$.part],
#'   cex = 0.5,
#'   pch = 19
#' )
#'
#'
#' # Now lets learn use these functions with several species
#' spp2 <- split(spp, spp$species)
#' class(spp2)
#' length(spp2)
#' names(spp2)
#'
#' part_list <- lapply(spp2, function(x) {
#'   result <- part_sblock(
#'     env_layer = somevar,
#'     data = x,
#'     x = "x",
#'     y = "y",
#'     pr_ab = "pr_ab",
#'     min_res_mult = 10,
#'     max_res_mult = 500,
#'     num_grids = 30,
#'     n_part = 2,
#'     prop = 0.5
#'   )
#'   result
#' })
#'
#' # Lets create a single database for all species
#' occ_part <- dplyr::bind_rows(lapply(
#'   part_list,
#'   function(x) x[[1]]
#' ), .id = "species")
#' occ_part
#'
#' # Lets get a the best grid info for all species
#' grid_info <- dplyr::bind_rows(lapply(
#'   part_list,
#'   function(x) x[[2]]
#' ), .id = "species")
#'
#' # Lets get a the best grid layer for all species
#' grid_layer <- lapply(part_list, function(x) x$grid)
#' grid_layer2 <-
#'   lapply(grid_layer, function(x) {
#'     get_block(env_layer = somevar[[1]], best_grid = x)
#'   })
#' grid_layer2 <- terra::rast(grid_layer2)
#' grid_layer2
#' plot(grid_layer2)
#'
#'
#' # Block partition for presences-only database
#' single_spp <- spp %>%
#'   dplyr::filter(species == "sp2", pr_ab == 1)
#' single_spp
#' single_spp$pr_ab %>% unique() # only presences
#'
#' part <- part_sblock(
#'   env_layer = somevar,
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   min_res_mult = 10,
#'   max_res_mult = 500,
#'   num_grids = 30,
#'   n_part = 2,
#'   prop = 0.5
#' )
#'
#' part$part
#' part$best_part_info
#' part$grid
#'
#' plot(part$grid)
#' points(
#'   part$part[c("x", "y")],
#'   col = c("blue", "red")[part$part$.part],
#'   cex = 0.5,
#'   pch = 19
#' )
#' }
#'
part_sblock <- function(env_layer,
                        data,
                        x,
                        y,
                        pr_ab,
                        n_part = 3,
                        min_res_mult = 3,
                        max_res_mult = 200,
                        num_grids = 30,
                        prop = 0.5) {

  # Select columns
  data <- dplyr::tibble(data)
  data <- data[, c(pr_ab, x, y)]
  colnames(data) <- c("pr_ab", "x", "y")

  if (any(!unique(data[, "pr_ab"][[1]]) %in% c(0, 1))) {
    stop(
      "values in pr_ab column did not match with 0 and 1:
unique list values in pr_ab column are: ",
      paste(unique(data[, "pr_ab"]), collapse = " ")
    )
  }

  # Extract data
  data <- dplyr::tibble(data, terra::extract(env_layer, data[, 2:3])[-1])
  filt <- stats::complete.cases(data)
  if (sum(!filt) > 0) {
    data <- data[filt, ]
    message(sum(!filt), " rows were excluded from database because NAs were found")
  }
  rm(filt)

  # Vector with presences and absences
  pa <- data %>%
    dplyr::pull(pr_ab)

  # Vector with grid cell-size used
  cell_size <- seq(terra::res(env_layer[[1]])[1] * min_res_mult,
    terra::res(env_layer[[1]])[1] * max_res_mult,
    length.out = num_grids
  )

  message(
    "The following grid cell sizes will be tested:\n",
    paste(round(cell_size, 2), collapse = " | "),
    "\n"
  )

  # Mask
  message("Creating basic raster mask...\n")

  mask <- env_layer[[1]]
  names(mask) <- "group"
  mask[!is.na(mask)] <- 1


  # Extent
  e <- terra::ext(mask)

  # Start Cluster
  message("Searching for the optimal grid size...\n")

  # Extract coordinates----
  mask2 <- mask
  mask2[] <- 0

  # Eliminate any records with NA
  filt <- stats::complete.cases(data[, c("x", "y")])
  if (any(!filt)) {
    message(sum(!filt), " rows were excluded from database because NAs were found in coordinates")
    data <- data[filt, ]
  }

  rm(filt)
  presences2 <- data

  # Transform the presences points in a DataFrameSpatialPoints
  presences2 <- terra::vect(presences2, geom = c("x", "y"), crs = terra::crs(mask))

  #### Data partitioning using a grid approach ####

  # Create a list of grids based on different raster resolution
  grid <- list() # List of grids

  # raster resolution
  DIM <-
    matrix(0, length(cell_size), 2) # the number of rows and columns of each grid
  colnames(DIM) <- c("R", "C")

  for (i in 1:length(cell_size)) {
    mask3 <- mask2
    terra::res(mask3) <- cell_size[i]
    mask3 <- terra::extend(mask3, y = c(1, 1))
    DIM[i, ] <- dim(mask3)[1:2]
    terra::values(mask3) <- 1 # Add values to cells /
    NAS <-
      c(terra::extract(mask3, presences2)[-1]) # Extract values to test if exist NAs
    if (any(is.na(NAS))) {
      while (any(is.na(NAS))) {
        terra::ext(mask3) <- terra::ext(mask3) + cell_size[i]
        terra::res(mask3) <- cell_size[i] # Give to cells a size
        DIM[i, ] <- dim(mask3)[1:2]
        terra::values(mask3) <- 1
        NAS <- terra::extract(mask3, presences2)[-1]
      }
    }
    grid[[i]] <- mask3
  }
  rm(list = c("mask3", "mask2", "mask"))

  # In this section is assigned the group of each cell
  for (i in 1:length(grid)) {
    if (n_part %% 2 == 0) {
      group <- c(
        rep(1:n_part, DIM[i, 2])[1:DIM[i, 2]],
        rep(c((n_part / 2 + 1):n_part, 1:(n_part / 2)), DIM[i, 2])[1:DIM[i, 2]]
      )
    }
    if (n_part %% 2 == 1) {
      group <- c(
        rep(1:n_part, DIM[i, 2])[1:DIM[i, 2]],
        rep(c((n_part / 3 + 1):n_part, 1:(n_part / 2)), DIM[i, 2])[1:DIM[i, 2]]
      )
    }
    terra::values(grid[[i]]) <- rep(group, length.out = terra::ncell(grid[[i]]))
  }

  # Matrix within each columns represent the partitions of points
  # for each grid resolution
  part <- data.frame(matrix(0, nrow(presences2), length(grid)))
  for (i in 1:length(grid)) {
    part[, i] <- terra::extract(grid[[i]], presences2)[, 2]
  }
  part <- dplyr::tibble(part)

  ### Remove problematic grids based on presences
  # Grids that assigned partitions less than the number of groups will be removed
  pp <- sapply(part[pa == 1, ], function(x) {
    length(unique(x))
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

  cell_size <- cell_size[pp]
  grid <- grid[pp]
  part <- part[, pp]
  names(part) <- names(which(pp == TRUE))

  ### Remove problematic grids based on presences
  # Grids that assigned partitions less than the number of groups will be removed
  if (any(unique(pa) == 0)) {
    pa <- presences2$pr_ab # Vector with presences and absences
    pp <- sapply(part[pa == 0, ], function(x) {
      length(unique(x))
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

    cell_size <- cell_size[pp]
    grid <- grid[pp]
    part <- part[, pp]
    names(part) <- names(which(pp == TRUE))
  }


  if (ncol(part) == 0) {
    message("It was not possible to find a good partition. Try to change values in 'n_part', or in 'min_res_mult', 'max_res_mult', or 'num_grids'")
    return(NA)
  }

  # Ncell
  ncell <- data.frame(matrix(
    0, nrow(presences2),
    length(grid)
  ))
  for (i in 1:length(grid)) {
    ncell[, i] <- terra::cellFromXY(grid[[i]], terra::geom(presences2)[, c("x", "y")])
  }

  # Performance of cells ----
  # SD of number of records per cell size-----

  sd_p <- rep(NA, length(grid))
  if (any(unique(pa) == 0)) {
    sd_a <- rep(NA, length(grid))
  }

  for (i in 1:ncol(part)) {
    if (any(unique(pa) == 0)) {
      sd_a[i] <- stats::sd(table(part[pa == 0, i]))
    }
    sd_p[i] <- stats::sd(table(part[pa == 1, i]))
  }

  # Environmental similarity between train and test based on euclidean  -----
  Env.P <- terra::extract(env_layer, presences2)[-1]
  env_sim <- rep(NA, length(grid))
  for (i in 1:ncol(part)) {
    cmb <- unique(part[, i][[1]]) %>% combn(2)
    Env.P1 <- cbind(part[i], Env.P)
    Env.P1 <- Env.P1[complete.cases(Env.P1), ]
    Env.P1 <- split(Env.P1[, -1], Env.P1[, 1])
    euq_c <- list()
    for (r in 1:ncol(cmb)) {
      euq_c[[r]] <- flexclust::dist2(Env.P1[[cmb[1, r]]], Env.P1[[cmb[2, r]]]) %>% mean()
    }

    env_sim[i] <- euq_c %>%
      unlist() %>%
      mean()
    rm(list = c("Env.P1"))
  }


  # I moran-----
  spa_auto <- rep(NA, length(grid))

  dist <- flexclust::dist2(
    presences2[, c("x", "y")] %>% as.data.frame(),
    presences2[, c("x", "y")] %>% as.data.frame()
  )
  dist <- 1 / dist
  diag(dist) <- 0
  dist[which(dist == Inf)] <- 0

  for (p in 1:ncol(part)) {
    cmb <- unique(part[, p][[1]]) %>% utils::combn(2)
    imoran_grid_c <- rep(NA, ncol(cmb))
    dff <- dplyr::tibble(nrow = 1:nrow(part), data["pr_ab"], group = part[p][[1]])

    for (c in 1:ncol(cmb)) {
      filt <- dff %>%
        dplyr::group_by(group, pr_ab) %>%
        dplyr::slice_sample(prop = prop) %>%
        dplyr::pull(nrow) %>%
        sort()

      odd <- which((part[p][[1]] == cmb[1, c])[filt])
      even <- which((part[p][[1]] == cmb[2, c])[filt])

      dist2 <- dist[filt, filt]
      dist2[odd, odd] <- 0
      dist2[even, even] <- 0

      mins <- apply(dist2, 2, function(x) {
        max(x, na.rm = TRUE)
      })
      for (i in 1:length(mins)) {
        dist2[, i] <- ifelse(dist2[, i] == mins[i], mins[i], 0)
      }

      if (nrow(data) < 3) {
        imoran_grid_c[c] <- NA
      } else {
        im <- sapply(
          data[filt, names(env_layer)],
          function(x) {
            ape::Moran.I(x,
              dist2,
              na.rm = TRUE,
              scaled = TRUE
            )$observed
          }
        )
        imoran_grid_c[c] <- mean(abs(im))
      }
    }
    spa_auto[p] <- mean(imoran_grid_c)
  }

  Opt <-
    if (any(unique(pa) == 0)) {
      data.frame(n_grid = 1:length(cell_size), cell_size = cell_size, round(
        data.frame(spa_auto, env_sim, sd_p, sd_a),
        3
      ))
    } else {
      data.frame(n_grid = 1:length(cell_size), cell_size = cell_size, round(
        data.frame(spa_auto, env_sim, sd_p),
        3
      ))
    }

  # Cleaning those variances based in data divided in a number of partition less than
  # the number of groups

  # SELLECTION OF THE BEST CELL SIZE----
  Opt2 <- Opt
  rownames(Opt2) <- colnames(part)
  Dup <-
    if (any(unique(pa) == 0)) {
      !duplicated(Opt2[c("spa_auto", "env_sim", "sd_p", "sd_a")])
    } else {
      !duplicated(Opt2[c("spa_auto", "env_sim", "sd_p")])
    }

  Opt2 <- Opt2[Dup, ]

  while (nrow(Opt2) > 1) {
    # I MORAN
    if (nrow(Opt2) == 1) {
      break
    }
    Opt2 <-
      Opt2[which(Opt2$spa_auto <= summary(Opt2$spa_auto)[2]), ]
    if (nrow(Opt2) == 1) {
      break
    }
    # Euclidean
    Opt2 <-
      Opt2[which(Opt2$env_sim >= summary(Opt2$env_sim)[5]), ]
    if (nrow(Opt2) == 1) {
      break
    }
    # SD presence
    Opt2 <-
      Opt2[which(Opt2$sd_p <= summary(Opt2$sd_p)[2]), ]
    if (nrow(Opt2) == 2) {
      break
    }
    # SD absences
    if (any(unique(pa) == 0)) {
      Opt2 <-
        Opt2[which(Opt2$sd_a <= summary(Opt2$sd_a)[2]), ]
      if (nrow(Opt2) == 2) {
        break
      }
    }

    if (unique(Opt2$spa_auto) &&
      unique(Opt2$env_sim) && unique(Opt2$sd_p)) {
      Opt2 <- Opt2[nrow(Opt2), ]
    }
  }

  if (nrow(Opt2) > 1) {
    Opt2 <- Opt2[nrow(Opt2), ]
  }


  # Final data.frame result----
  result <- data.frame(data, .part = c(part[, rownames(Opt2)])[[1]])
  result <- result %>% dplyr::select(-names(env_layer))
  colnames(result) <- c("pr_ab", "x", "y", ".part")
  result <- result[c("x", "y", "pr_ab", ".part")]

  grid <- grid[[Opt2$n_grid]]
  names(grid) <- "block"

  # Final data.frame result2----
  out <- list(
    part = dplyr::tibble(result),
    best_part_info = dplyr::tibble(Opt2),
    grid = grid # Optimum size for presences
  )
  return(out)
}
