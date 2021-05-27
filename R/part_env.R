#' Environmental cross validation
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
#' @param min_n_groups
#' @param min_res_mult integer. Minimum number of groups to be tested. Default 2.
#' @param max_n_groups integer. Maximum number of groups to be tested. Default 10.
#' @param prop numeric. Proportion of point used for testing autocorrelation between
#' groups (values > 0 and <=1). The smaller this number is, the faster the function will work.
#' Default 0.5
#'
#' @return
#' A list with:
#' \itemize{
#'   \item part: A tibble object with information used in 'data' arguments and a additional column .part with partition group.
#'   \item best_part_info: A tibble with information of the bets partition. It contains the number of partition (n_groups), standard deviation of presences (sd_p), standard deviation of absences (sd_a), Moran's I spatial autocorrelation (spa_auto) and environmental similarity based on euclidean distance (env_sim)
#'   }
#'
#' @export
#' @importFrom ape Moran.I
#' @importFrom dplyr %>% pull bind_cols tibble group_by count filter select slice_sample
#' @importFrom flexclust dist2
#' @importFrom stats kmeans
#' @importFrom terra extract
#'
#' @seealso \code{\link{part}}, and \code{\link{part_spat}}.
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
#' part1 <- part_env(
#'   env_layer = somevar,
#'   data = spp1,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   min_n_groups = 2,
#'   max_n_groups = 10,
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
part_env <- function(env_layer,
                     data,
                     x,
                     y,
                     pr_ab,
                     min_n_groups = 2,
                     max_n_groups = 10,
                     prop = 0.5) {
  # Select columns
  data <- data.frame(data)
  data <- data[, c(pr_ab, x, y)]
  colnames(data) <- c("pr_ab", "x", "y")
  pa <- data %>%
    dplyr::pull(pr_ab) %>%
    unique() # Test if there are absences

  # Test some about 0 and 1 (TRY TO ADAPT THIS FUNCTION FOR WORKING ONLY FOR PRESENCE)
  if (any(!unique(data[, "pr_ab"]) %in% c(0, 1))) {
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

  # Vector with grid cell-size used
  cell_size <- seq(min_n_groups, max_n_groups)

  message(
    "The following grid cell sizes will be tested:\n",
    paste(round(cell_size, 2), collapse = " | "),
    "\n"
  )

  message("Searching best partition...\n")

  # k-mean algorithm
  km_l <- list()
  for (i in 1:length(cell_size)) {
    km_l[[i]] <- stats::kmeans(data[, -1], centers = cell_size[i])$cluster
  }

  names(km_l) <- paste0(".g", cell_size)

  # Bind groups
  km_l <- dplyr::bind_cols(km_l)


  ### Remove problematic partition
  n_records <- apply(km_l, 2, function(x) {
    dplyr::tibble(x, data[1]) %>%
      dplyr::group_by(x, pr_ab) %>%
      dplyr::count() %>%
      dplyr::mutate(filt = n <= 2) # 2
  })

  filt <- sapply(n_records, function(x) any(x %>% dplyr::pull("filt")))

  if (sum(!filt) == 0) {
    stop("it was not possible to find a good partition")
  }

  n_records <- n_records[!filt]
  km_l <- km_l[!filt]

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

  env_sim <- rep(NA, ncol(km_l))
  for (i in 1:ncol(km_l)) {
    cmb <- unique(km_l[, i][[1]]) %>% combn(2)

    Env.P1 <- cbind(km_l[i], Env.P)
    Env.P1 <- split(Env.P1[, -1], Env.P1[, 1])
    euq_c <- list()
    for (r in 1:ncol(cmb)) {
      euq_c[[r]] <- flexclust::dist2(Env.P1[[cmb[1, r]]], Env.P1[[cmb[2, r]]]) %>% mean()
    }

    env_sim[i] <- euq_c %>%
      unlist() %>%
      mean()
  }


  # # I moran-----
  spa_auto <- rep(NA, ncol(km_l))

  dist <- flexclust::dist2(
    data[, c("x", "y")] %>% as.data.frame(),
    data[, c("x", "y")] %>% as.data.frame()
  )
  dist <- 1 / dist
  diag(dist) <- 0
  dist[which(dist == Inf)] <- 0

  for (p in 1:ncol(km_l)) {
    cmb <- unique(km_l[, p][[1]]) %>% utils::combn(2)
    imoran_grid_c <- rep(NA, ncol(cmb))
    dff <- dplyr::tibble(nrow = 1:nrow(km_l), data["pr_ab"], group = km_l[p][[1]])

    for (c in 1:ncol(cmb)) {
      filt <- dff %>%
        dplyr::group_by(group, pr_ab) %>%
        dplyr::slice_sample(prop = prop) %>%
        dplyr::pull(nrow) %>%
        sort()

      odd <- which((km_l[p][[1]] == cmb[1, c])[filt])
      even <- which((km_l[p][[1]] == cmb[2, c])[filt])

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
        imoran_grid_c[c] <- mean(im)
      }
    }
    spa_auto[p] <- mean(imoran_grid_c)
  }

  spa_auto <- abs(spa_auto)

  # # SELLECTION OF THE BEST CELL SIZE----
  Opt <-
    if (any(unique(pa) == 0)) {
      data.frame(
        n_parition = 1:ncol(km_l),
        n_groups = gsub(".g", "", colnames(km_l)),
        round(
          data.frame(sd_p, sd_a, spa_auto, env_sim),
          3
        )
      )
    } else {
      data.frame(n_groups = gsub(".g", "", colnames(km_l)), round(
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
    if (any(pa == 0)) {
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

    if (unique(Opt2$spa_auto) &&
      unique(Opt2$env_sim) && unique(Opt2$sd_p)) {
      Opt2 <- Opt2[nrow(Opt2), ]
    }
  }


  # Final data.frame result----
  result <- data.frame(data, .part = c(km_l[, rownames(Opt2)]))
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
