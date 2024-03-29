#' Environmental and spatial cross-validation
#'
#' @description This function explores different numbers of environmental partitions (clusters)
#' based on the K-means clustering algorithm and returns the number of partitions best suited for
#' a given presence, presence-absences, or presence-pseudo-absences database. Selection of the
#' best number of partitions is performed automatically considering spatial autocorrelation,
#' environmental similarity, and the number of presence and/or absence records in each partition.
#'
#' @param env_layer SpatRaster. Raster with environmental
#' variable. This will be used to evaluate spatial autocorrelation and
#' environmental similarity between training and testing partitions. Because this function
#' calculate dissimilarity based on Euclidean distances, it can only be used with continuous
#' variables
#' @param data data.frame. Data.frame or tibble object with presence
#' (or presence-absence, or presences-pseudo-absence) records, and coordinates
#' @param x character. Column name with spatial x coordinates
#' @param y character. Column name with spatial y coordinates
#' @param pr_ab character. Column with presences, presence-absence,
#' or pseudo-absence. Presences must be represented by 1 and absences by 0
#' @param min_n_groups integer. Minimum number of groups to be tested. Default 2.
#' @param max_n_groups integer. Maximum number of groups to be tested. Default 10.
#' @param min_occ numeric. Minimum number of presences or absences in a partition fold.
#' The min_occ value should be base on the amount of predictors in order to avoid over-fitting
#' or error when fitting models for a given fold. Default 10.
#' @param prop numeric. Proportion of point used for testing autocorrelation between
#' groups (values > 0 and <=1). The smaller this number is, the faster the function will work.
#' Default 0.5
#'
#' @return
#' A list with:
#' \itemize{
#'   \item part: A tibble object with information used in 'data' arguments and a additional column
#'   .part with partition group.
#'   \item best_part_info: A tibble with information about the best partition. It contains the
#'   number of partition (n_groups), standard deviation of presences (sd_p), standard deviation
#'   of absences (sd_a), Moran's I spatial autocorrelation (spa_auto) and environmental similarity
#'   based on Euclidean distance (env_sim)
#'   }
#'
#' @details The part_sblock allows test with different numbers of partitions defined in the
#' envirnomental clusters delimited the K-mean cluster algorithm. This function explores a range
#' of environmental clusters and automatically selects best number of cluster for a given given
#' presence, presence-absences, or presence-pseudo-absences dataset. Such selection of number of
#' clusters is based on an optimization procedure that explores partition size in three dimensions
#' determined by spatial autocorrelation (measured by Moran's I), environmental similarity
#' (Euclidean distance), and difference in the amount of data among clusters
#' (Standard Deviation - SD; Velazco et al., 2019). This procedure will cyclically select
#' those partitions with autocorrelation values less than the lowest quartile of Morans I, then
#' those with environmental similarity values greater than the third quartile of the Euclidean
#' distances than those with a difference in the amount of data less than the lowest quartile of SD.
#' This selection is repeated until only one partition is retained (Velazco et al., 2019). The
#' main benefit of this partition selection are i) this is not subjective, ii) balances the
#' environmental similarity and special autocorrelation between partitions, and iii) controls
#' the partition selection with few data that may be problematic for model fitting
#' ("min_occ" argument)..
#'
#' Partitions geographically structured tend to evaluate model transferability more directly than
#' conventional ones (e.g., those performed by \code{\link{part_random}}) (Roberts et al., 2017;
#' Santini et al., 2021), being relevant for models that want to be used for projections in other
#' regions outside the calibration area or for other periods.
#'
#' This function can interact with \code{\link{get_block}}, \code{\link{sample_background}},
#' and \code{\link{sample_pseudoabs}} for sampling background points or pseudo-absences within
#' spatial partition broups
#'
#' @references
#' \itemize{
#' \item Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J., Guillera-Arroita, G.,
#' Hauenstein, S., Lahoz-Monfort, J. J., Schroder, B., Thuiller, W., Warton, D. I., Wintle, B. A.,
#' Hartig, F., & Dormann, C. F. (2017). Cross-validation strategies for data with temporal, spatial,
#'  hierarchical, or phylogenetic structure. Ecography, 40,
#'  913-929. https://doi.org/10.1111/ecog.02881
#' \item Santini, L., Benitez-Lopez, A., Maiorano, L., Cengic, M., & Huijbregts, M. A. J. (2021).
#'  Assessing the reliability of species distribution projections in climate change research.
#'  Diversity and Distributions, ddi.13252. https://doi.org/10.1111/ddi.13252
#' \item Velazco, S. J. E., Villalobos, F., Galvao, F., & De Marco Junior, P. (2019). A dark
#' scenario for Cerrado plant species: Effects of future climate, land use and protected areas
#' ineffectiveness. Diversity and Distributions, 25(4), 660-673. https://doi.org/10.1111/ddi.12886
#' }
#'
#' @export
#' @importFrom dplyr tibble pull bind_cols group_by count mutate filter select slice_sample
#' @importFrom stats complete.cases kmeans sd
#' @importFrom terra extract
#' @importFrom utils combn
#'
#' @seealso \code{\link{part_random}}, \code{\link{part_sblock}}, and \code{\link{part_sband}}
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(ggplot2)
#'
#' f <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(f)
#'
#' # Select a species
#' spp1 <- spp %>% dplyr::filter(species == "sp1")
#'
#' part1 <- part_senv(
#'   env_layer = somevar,
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   min_n_groups = 2,
#'   max_n_groups = 10,
#'   min_occ = 10,
#'   prop = 0.2
#' )
#'
#' part1
#'
#' ggplot(part1$part, aes(x, y, col = factor(.part))) +
#'   geom_point(aes(shape = factor(pr_ab)))
#'
#' ggplot(part1$part, aes(x, y, col = factor(.part))) +
#'   geom_point(aes(shape = factor(pr_ab))) +
#'   facet_wrap(. ~ .part)
#'
#' ggplot(part1$part, aes(x, y, col = factor(.part))) +
#'   geom_point(aes(shape = factor(pr_ab))) +
#'   facet_wrap(. ~ pr_ab)
#' }
part_senv <- function(env_layer,
                      data,
                      x,
                      y,
                      pr_ab,
                      min_n_groups = 2,
                      max_n_groups = 10,
                      min_occ = 10,
                      prop = 0.5) {
  group <- NULL
  # Select columns
  data <- data.frame(data)
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
    dplyr::pull(pr_ab) %>%
    unique()

  # Vector with grid cell-size used
  cell_size <- seq(min_n_groups, max_n_groups)

  message(
    "The following grid cell sizes will be tested:\n",
    paste(round(cell_size, 2), collapse = " | "),
    "\n"
  )

  message("Searching best partition...\n")

  # k-mean algorithm
  part <- list()
  for (i in 1:length(cell_size)) {
    part[[i]] <- stats::kmeans(data[, -1], centers = cell_size[i])$cluster
  }

  names(part) <- paste0(".g", cell_size)

  # Bind groups
  part <- dplyr::bind_cols(part)


  ### Remove problematic partition
  n_records <- apply(part, 2, function(x) {
    dplyr::tibble(x, data[1]) %>%
      dplyr::group_by(x, pr_ab) %>%
      dplyr::count() %>%
      dplyr::mutate(filt = (n <= min_occ)) # 2
  })

  filt <- sapply(n_records, function(x) any(x %>% dplyr::pull("filt")))

  if (sum(!filt) == 0) {
    message("It was not possible to find a good partition. Try to change values in 'min_n_groups' and 'max_n_groups'")
    return(NA)
  }

  n_records <- n_records[!filt]
  part <- part[!filt]

  # Perform SD
  sd_p <- sapply(n_records, function(x) {
    x %>%
      dplyr::filter(pr_ab == 1) %>%
      dplyr::pull(n) %>%
      stats::sd()
  })

  sd_a <- sapply(n_records, function(x) {
    x %>%
      dplyr::filter(pr_ab == 0) %>%
      dplyr::pull(n) %>%
      sd()
  })

  rm(n_records)

  # # Environmental similarity between train and test based on euclidean  -----
  Env.P <- data %>% dplyr::select(-pr_ab, -x, -y)
  env_sim <- rep(NA, ncol(part))
  for (i in 1:ncol(part)) {
    cmb <- unique(part[, i][[1]]) %>% utils::combn(2)

    Env.P1 <- cbind(part[i], Env.P)
    Env.P1 <- Env.P1[stats::complete.cases(Env.P1), ]
    Env.P1 <- stats::complete.cases(Env.P1[, -1], Env.P1[, 1])
    euq_c <- list()
    for (r in 1:ncol(cmb)) {
      euq_c[[r]] <- euc_dist(Env.P1[[cmb[1, r]]], Env.P1[[cmb[2, r]]]) %>% mean()
    }

    env_sim[i] <- euq_c %>%
      unlist() %>%
      mean()
  }


  # # I moran-----
  spa_auto <- rep(NA, ncol(part))

  dist <- euc_dist(
    data[, c("x", "y")] %>% as.data.frame(),
    data[, c("x", "y")] %>% as.data.frame()
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
            morani(x,
              dist2,
              na.rm = TRUE,
              scaled = TRUE
            )
          }
        )
        imoran_grid_c[c] <- mean(abs(im))
      }
    }
    spa_auto[p] <- mean(imoran_grid_c)
  }

  # # SELLECTION OF THE BEST CELL SIZE----
  Opt <-
    if (any(unique(pa) == 0)) {
      data.frame(
        n_parition = 1:ncol(part),
        n_groups = gsub(".g", "", colnames(part)),
        round(
          data.frame(sd_p, sd_a, spa_auto, env_sim),
          3
        )
      )
    } else {
      data.frame(n_groups = gsub(".g", "", colnames(part)), round(
        data.frame(
          sd_p, spa_auto, env_sim
        ),
        3
      ))
    }

  Opt2 <- Opt

  while (nrow(Opt2) > 1) {
    # SD presence
    Opt2 <-
      Opt2[which(Opt2$sd_p <= summary(Opt2$sd_p)[2]), ]
    if (nrow(Opt2) == 1) {
      break
    }
    # SD absences
    if (any(unique(pa) == 0)) {
      Opt2 <-
        Opt2[which(Opt2$sd_a <= summary(Opt2$sd_a)[2]), ]
      if (nrow(Opt2) == 1) {
        break
      }
    }
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

    if ((length(unique(Opt2$spa_auto)) == 1) &&
      (length(unique(Opt2$env_sim)) == 1) &&
      (length(unique(Opt2$sd_p)) == 1)) {
      Opt2 <- Opt2[nrow(Opt2), ]
    }
  }


  # Final data.frame result----
  result <- data.frame(data, .part = c(part[, rownames(Opt2)])[[1]])
  result <- result %>% dplyr::select(-names(env_layer))
  colnames(result) <- c("pr_ab", "x", "y", ".part")
  result <- result[c("x", "y", "pr_ab", ".part")]

  # Final data.frame result2----
  result <- list(
    part = dplyr::tibble(result),
    best_part_info = dplyr::tibble(Opt2)
  )
  return(result)
}
