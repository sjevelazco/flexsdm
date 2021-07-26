#' Spatial band cross-validation
#'
#' @description This function explores different numbers of spatial bands and returns the
#' best one suited for a given presence or presence-absences database. The selection of the best
#' number of bands is performed automatically considering spatial autocorrelation, environmental
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
#' @param type character. Specify bands across different degrees of longitude 'lon' or latitude 'lat'. Default is 'lon'.
#' @param n_part  integer. Number of partition. Default 2, values other than
#' 2 has not yet been implemented.
#' @param min_bands integer. Minimum number of spatial bands to be tested, default 2.
#' @param max_bands integer. Maximum number of spatial bands to be tested, default 20.
#' @param prop numeric. Proportion of point used for testing autocorrelation between
#' groups (values > 0 and <=1). The smaller this number is, the faster the function will work.
#' Default 0.5
#'
#' @return
#' A list with:
#' \itemize{
#'   \item part: A tibble object with information used in 'data' arguments and a additional column .part with partition group.
#'   \item best_part_info: A tibble with information of the bets partition. It contains the number of the best partition (n_grid), number of bands (n_bands), standard deviation of presences (sd_p), standard deviation of absences (sd_a), Moran's I spatial autocorrelation (spa_auto), and environmental similarity based on euclidean distance (env_sim).
#'   \item grid: A SpatRaster object with bands
#'   }
#' @export
#'
#' @importFrom ape Moran.I
#' @importFrom dplyr tibble pull group_by slice_sample select
#' @importFrom flexclust dist2
#' @importFrom stats complete.cases sd
#' @importFrom terra extract res ext vect crs extend values ncell cellFromXY geom
#' @importFrom utils combn
#'
#' @examples
#'
#' \dontrun{
#' require(terra)
#' require(dplyr)
#'
#' # Load datasets
#' data(spp)
#' f <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(f)
#'
#' # Let's practice with a single species
#' single_spp <- spp %>% dplyr::filter(species == "sp3")
#' part <- part_sband(
#'   env_layer = somevar,
#'   data = single_spp,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   type = "lon",
#'   min_bands = 2,
#'   max_bands = 20,
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
#' }
#'
#'
part_sband <- function(env_layer,
                        data,
                        x,
                        y,
                        pr_ab,
                        type = 'lon',
                        n_part = 2,
                        min_bands = 2,
                        max_bands = 20,
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

  # Vector with number of bands to be tested
  n_bands <- seq(min_bands, max_bands)

  message(
    "The following number of bands will be tested:\n",
    paste(round(n_bands, 2), collapse = " | "),
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
  message("Searching for the optimal number of bands...\n")

  # Extract coordinates----
  mask2 <- mask
  mask2[] <- 0
  presences2 <- data

  # Transform the presences points to SpatVect
  presences2 <- terra::vect(presences2, geom = c("x", "y"), crs = terra::crs(mask))

  #### Data partitioning using a grid approach ####

  # Create a list of bands based on user input
  grid <- list() # List of bands

  # Longitude or latitude?

  DIM <-
    matrix(0, length(n_bands), 2) # Rows and columns
   colnames(DIM) <- c("R", "C")

  if (type == 'lon') {
    DIM[, 1] <- n_bands
    DIM[, 2] <- rep(1)

    for (i in 1:length(n_bands)) {
      mask3 <- mask2
      terra::ncol(mask3) <- n_bands[i]
      terra::nrow(mask3) <- 1
      DIM[i, ] <- dim(mask3)[1:2]
      terra::values(mask3) <- 1 # Add values to cells /
      NAS <-
        c(terra::extract(mask3, presences2)[-1]) # Extract values to test if exist NAs
      if (any(is.na(NAS))) {
        while (any(is.na(NAS))) {
          DIM[i, ] <- dim(mask3)[1:2]
          terra::values(mask3) <- 1
          NAS <- terra::extract(mask3, presences2)[-1]
        }
      }
      grid[[i]] <- mask3
    }
    rm(list = c("mask3", "mask2", "mask"))

  }

  if (type == 'lat'){
    DIM[, 1] <- rep(1)
    DIM[, 2] <- n_bands

    for (i in 1:length(n_bands)) {
      mask3 <- mask2
      terra::ncol(mask3) <- 1
      terra::nrow(mask3) <- n_bands[i]
      DIM[i, ] <- dim(mask3)[1:2]
      terra::values(mask3) <- 1 # Add values to cells /
      NAS <-
        c(terra::extract(mask3, presences2)[-1]) # Extract values to test if exist NAs
      if (any(is.na(NAS))) {
        while (any(is.na(NAS))) {
          DIM[i, ] <- dim(mask3)[1:2]
          terra::values(mask3) <- 1
          NAS <- terra::extract(mask3, presences2)[-1]
        }
      }
      grid[[i]] <- mask3
    }
    rm(list = c("mask3", "mask2", "mask"))

  }


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
  # for each option for number of bands
  part <- data.frame(matrix(0, nrow(presences2), length(grid)))
  for (i in 1:length(grid)) {
    part[, i] <- terra::extract(grid[[i]], presences2)[, 2]
  }
  part <- dplyr::tibble(part)

  ### Remove problematic bands based on presences
  # Band options that are assigned partitions less than the number of groups will be removed
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

  n_bands <- n_bands[pp]
  grid <- grid[pp]
  part <- part[, pp]
  names(part) <- names(which(pp == TRUE))

  if (ncol(part) == 0) {
    message("It was not possible to find a good partition. Try to change values in 'n_part', or in 'min_band', 'max_band', or 'num_bands'")
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
        imoran_bands_c[c] <- NA
      } else {
        im <- sapply(
          data[filt, names(env_layer)],
          function(x) {
            suppressMessages(
              ape::Moran.I(x,
                           dist2,
                           na.rm = TRUE,
                           scaled = TRUE
              )$observed
            )
          }
        )
        imoran_grid_c[c] <- mean(abs(im))
      }
    }
    spa_auto[p] <- mean(imoran_grid_c)
  }

  Opt <-
    if (any(unique(pa) == 0)) {
      data.frame(n_grid = 1:length(n_bands), n_bands = n_bands, round(
        data.frame(spa_auto, env_sim, sd_p, sd_a),
        3
      ))
    } else {
      data.frame(n_grid = 1:length(n_bands), n_bands = n_bands, round(
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
