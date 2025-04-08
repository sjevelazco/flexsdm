#' Calculate environmental distance between presences and projection data
#'
#' @param training_data data.frame or tibble with environmental conditions of
#' presence used for constructing models
#' @param projection_data SpatRaster, data.frame or tibble with environmental condition used for projecting a model (e.g.,
#' a larger, encompassing region, a spatially separate region, or a different time period).
#' If data.frame or tibble is used function will return a tibble object.
#' Otherwise, as SpatRaster object.
#' @param metric character. Metric used for measuring distance. Default "domain".
#'
#' @return
#' A SpatRaster or tibble object with the nearest environmental distance between presences and projection data.
#' So far only the Domain algorithm (based on the Gower distance; Carpenter et al., 1993) has been implemented.
#'
#' @export
#'
#' @references
#' \itemize{
#' \item Carpenter, G., Gillison, A.N., Winter, J., 1993. DOMAIN: a flexible modelling procedure for mapping potential distributions of plants and animals. Biodiversity & Conservation 2, 667–680
#' }
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' require(terra)
#' data(spp)
#' f <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(f)
#'
#' # Let's use only two variables to turn more evident the pater in the environmental space
#' somevar <- somevar[[1:2]]
#' names(somevar) <- c("aet", "cwd")
#'
#'
#' spp$species %>% unique()
#' sp <- spp %>%
#'   dplyr::filter(species == "sp3", pr_ab == 1) %>%
#'   dplyr::select(x, y, pr_ab)
#'
#' # Get environmental condition of presences
#' sp_pa_2 <- sdm_extract(
#'   data = sp,
#'   x = "x",
#'   y = "y",
#'   env_layer = somevar
#' )
#' sp_pa_2
#'
#' # Measure degree of extrapolation based on Mahalanobis and
#' # for a projection area based on a SpatRaster object
#' envdist <-
#'   map_env_dist(
#'     training_data = sp_pa_2,
#'     projection_data = somevar,
#'     metric = "domain"
#'   )
#' plot(envdist, main = "Domain")
#' p_extra(
#'   training_data = sp_pa_2,
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   extra_suit_data = envdist,
#'   projection_data = somevar,
#'   geo_space = FALSE,
#'   prop_points = 0.8,
#'   alpha_p = 0.9,
#'   color_p = "red",
#'   color_gradient = c("#000033", "#1400FF", "#C729D6", "#FF9C63", "#FFFF60")
#' )
#'
#' }
map_env_dist <- function(
    training_data,
    projection_data,
    metric = "domain") {
  cell <- y <- x <- cluster <- . <-  NULL

  # Gower distance

  if (!metric %in% c("domain")) {
    stop("metric argument must be used with 'domain'")
  }
  if (length(metric) > 1) {
    stop("metric argument must be used with 'domain'")
  }
  if (any("data.frame" == class(training_data))) {
    training_data_pr_ab <- training_data[names(projection_data)] %>% na.omit()
    training_data <- training_data_pr_ab[names(projection_data)] %>%
      na.omit()
  }
  v0 <- unique(c(names(training_data), names(projection_data)))
  v0 <- sort(v0)
  if (!all(c(
    all(names(training_data) %in% names(projection_data)),
    all(names(projection_data) %in% names(training_data))
  ))) {
    stop(
      "colnames of dataframes of env_records, training_data, and projection_data\n        do not match each other",
      "\nraster layers names:", "\n", paste(sort(unique(unlist(v0))),
        collapse = "\n"
      )
    )
  }
  if (any("data.frame" %in% class(training_data))) {
    training_data <- training_data[v0]
  } else {
    training_data <- training_data[[v0]]
  }
  if (any("SpatRaster" == class(projection_data))) {
    projection_data <- projection_data[[v0]]
    extraraster <- projection_data[[1]]
    extraraster[!is.na(extraraster)] <- 0
  } else {
    projection_data <- projection_data[v0]
  }
  env_calib2 <- terra::as.data.frame(training_data,
    xy = FALSE,
    na.rm = TRUE
  ) %>% dplyr::as_tibble()
  env_proj2 <- terra::as.data.frame(projection_data,
    xy = FALSE,
    na.rm = TRUE
  )
  if (any("SpatRaster" == class(projection_data))) {
    ncell <- rownames(env_proj2) %>% as.numeric()
  }

  env_proj2 <- env_proj2 %>% dplyr::as_tibble()
  set <- c(seq(1, nrow(env_proj2), 100), nrow(env_proj2) + 1)

  extra <- list()
  for (x in seq_len((length(set) - 1))) {
    rowset <- set[x]:(set[x + 1] - 1)
    if (metric == "domain") {
      extra[[x]] <- min_gower_rcpp(data1 = env_calib2[v0], data2 = env_proj2[rowset, v0])
    }
  }
  # Set zero any value below 0
  extra <- unlist(extra)
  # extra <- 1 - extra
  # extra[extra < 0] <- 0


  if (any("SpatRaster" == class(projection_data))) {
    env_proj2 <- data.frame(distance = extra)

    extraraster[ncell] <- env_proj2[, "distance"]
    names(extraraster) <- "domain"
    rm(extra)
    return(extraraster)
  } else {
    return(extra)
  }
}
