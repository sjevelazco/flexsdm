#' Measure model extrapolation
#'
#' @description Measure extrapolation comparing environmental data used for modeling calibration and area for
#' model projection. This function use the approach proposed by xxx et al., in prep
#' (EXPERIMENTAL)
#'
#' @param training_data SpatRaster or tibble with environmental conditions of the calibration area or the
#' presence and absence (or background points or pseudo-absences) used for constructing models
#' @param projection_data SpatRaster with environmental condition used for projecting a model (e.g., a larger, encompassing region, a spatially separate region, or a different time period)
#' @param n_cores numeric. Number of cores use for parallelization. Default 1
#' @param aggreg_factor positive integer. Aggregation factor expressed as number of cells in each
#'  direction to reduce raster resolution. Use value higher than 1 would be useful when
#'  measuring extrapolation using a raster with a high number of cells. The resolution of output will be
#'  the same as raster object used in 'projection_data' argument. Default 1, i.e., by default, no changes
#'  will be made to the resolution of the environmental variables.
#'
#'
#' @return
#' A SpatRaster object with extrapolation values measured in percentage of extrapolation (relative Euclidean distance)
#'
#' @seealso \code{\link{extra_truncate}}, \code{\link{p_extra}}, \code{\link{p_pdp}}
#'
#' @export
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr summarise_all
#' @importFrom foreach foreach "%dopar%"
#' @importFrom parallel makeCluster stopCluster
#' @importFrom stats sd
#' @importFrom terra mask aggregate as.data.frame resample
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' require(terra)
#'
#' data(spp)
#' f <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(f)
#' names(somevar) <- c("aet", "cwd", "tmx", "tmn")
#'
#'
#' spp$species %>% unique()
#' sp <- spp %>%
#'   dplyr::filter(species == "sp3", pr_ab == 1) %>%
#'   dplyr::select(x, y, pr_ab)
#'
#' # Calibration area based on some criterion such as dispersal ability
#' ca <- calib_area(sp, x = "x", y = "y", method = c("bmcp", width = 50000), crs = crs(somevar))
#'
#' plot(somevar[[1]])
#' points(sp)
#' plot(ca, add = T)
#'
#'
#' # Sampling pseudo-absences
#' set.seed(10)
#' psa <- sample_pseudoabs(
#'   data = sp,
#'   x = "x",
#'   y = "y",
#'   n = nrow(sp) * 2, # selecting number of pseudo-absence points twice number of presences
#'   method = "random",
#'   rlayer = somevar,
#'   calibarea = ca
#' )
#'
#' # Merge presences and abasences databases to get a complete calibration data
#' sp_pa <- dplyr::bind_rows(sp, psa)
#' sp_pa
#'
#' # Get environmental condition of calibration area
#' sp_pa_2 <- sdm_extract(data = sp_pa, x = "x", y = "y", env_layer = somevar)
#' sp_pa_2
#'
#' # Measure extrapolation based on calibration data (presence and pseudo-absences)
#' extr <-
#'   extra_eval(
#'     training_data = sp_pa_2,
#'     projection_data = somevar,
#'     n_cores = 1,
#'     aggreg_factor = 1
#'   )
#' plot(extr, main = "Extrapolation pattern")
#'
#'
#'
#' # Let's fit, predict and truncate a model with extra_truncate
#' sp_pa_2 <- part_random(
#'   data = sp_pa_2,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 5)
#' )
#'
#' a_model <- fit_glm(
#'   data = sp_pa_2,
#'   response = "pr_ab",
#'   predictors = c("aet", "cwd", "tmx", "tmn"),
#'   partition = ".part",
#'   thr = c("max_sorensen")
#' )
#'
#' predsuit <- sdm_predict(models = a_model, pred = somevar, thr = "max_sorensen")
#' predsuit # list with a raster with two layer
#' plot(predsuit[[1]])
#'
#' # Truncate a model based on a given value of extrapolation based on SHAPE metric
#' par(mfrow = c(1, 2))
#' plot(extr, main = "Extrapolation")
#' plot(predsuit[[1]][[1]], main = "Suitability")
#' par(mfrow = c(1, 1))
#'
#' predsuit_2 <- extra_truncate(
#'   suit = predsuit[[1]],
#'   extra = extr,
#'   threshold = c(50, 100, 200)
#' )
#' predsuit_2 # a list of continuous and binayr models with differnt truncated at different
#' # extrapolation thresholds
#'
#' plot(predsuit_2$`50`)
#' plot(predsuit_2$`100`)
#' plot(predsuit_2$`200`)
#'
#' # See also functions p_extra (to explore extrapolation and suitability paterns in the
#' # geographical and environmental space) and p_pdp to construct partial depence plots
#' }
extra_eval <- function(training_data, projection_data, n_cores = 1, aggreg_factor = 1) {
  Value <- val <- . <- x <- NULL

  if (any("data.frame" == class(training_data))) {
    training_data <- training_data[names(projection_data)]
  }

  # Get variable names
  v0 <- unique(c(names(training_data), names(projection_data)))
  v0 <- sort(v0)

  # Test rasters variable names
  if (!all(c(
    all(names(training_data) %in% names(projection_data)),
    all(names(projection_data) %in% names(training_data))
  ))) {
    stop(
      "colnames of dataframes of env_records, training_data, and projection_data
        do not match each other",
      "\nraster layers names:",
      "\n",
      paste(sort(unique(unlist(
        v0
      ))), collapse = "\n")
    )
  }


  # Sort in the same way layer in both raster
  if (any("data.frame" %in% class(training_data))) {
    training_data <- training_data[v0]
  } else {
    training_data <- training_data[[v0]]
  }
  projection_data <- projection_data[[v0]]

  # Layer base
  extraraster <- terra::mask(!is.na(projection_data[[1]]), projection_data[[1]])

  if (aggreg_factor == 1) {
    aggreg_factor <- NULL
  }
  if (!is.null(aggreg_factor)) {
    disag <- extraraster
    if (any("SpatRaster" == class(training_data))) {
      training_data <- terra::aggregate(training_data, fact = aggreg_factor, na.rm = TRUE)
    }
    projection_data <- terra::aggregate(projection_data, fact = aggreg_factor, na.rm = TRUE)
    extraraster <- terra::aggregate(extraraster, fact = aggreg_factor, na.rm = TRUE)
  }


  # Transform raster to df
  env_calib2 <-
    terra::as.data.frame(training_data, xy = FALSE, na.rm = TRUE)
  env_proj2 <-
    terra::as.data.frame(projection_data, xy = FALSE, na.rm = TRUE)

  # save coordinates and cell number
  ncell <- rownames(env_proj2) %>% as.numeric()

  # Standardization
  # standardization based on calibration area
  s_center <- colMeans(env_calib2)
  s_scale <- apply(env_calib2, 2, stats::sd)

  for (i in 1:ncol(env_calib2)) {
    env_calib2[i] <- (env_calib2[i] - s_center[i]) / s_scale[i]
  }
  for (i in 1:ncol(env_proj2)) {
    env_proj2[i] <- (env_proj2[i] - s_center[i]) / s_scale[i]
  }

  # Measure extrapolation - Euclidean distance
  set <-
    c(
      seq(
        1, nrow(env_proj2),
        200
      ),
      nrow(env_proj2) + 1
    )

  if (n_cores == 1) {
    extra <- lapply(seq_len((length(set) - 1)), function(x) {
      rowset <- set[x]:(set[x + 1] - 1)
      auclidean <-
        euc_dist(
          env_proj2[rowset, v0], # env_proj2 prediction environmental conditions outside
          env_calib2[v0]
        ) # training_data environmental conditions used as references
      auclidean <- sapply(data.frame(t(auclidean)), min)
      return(auclidean)
    })
  } else {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    extra <- foreach::foreach(x = seq_len((length(
      set
    ) - 1)), .export = "euc_dist", .combine = "c") %dopar% {
      rowset <- set[x]:(set[x + 1] - 1)
      auclidean <-
        euc_dist(
          env_proj2[rowset, v0], # env_proj2 prediction environmental conditions outside
          env_calib2[v0]
        ) # training_data environmental conditions used as references
      auclidean <- sapply(data.frame(t(auclidean)), min)
      auclidean
    }

    parallel::stopCluster(cl)
  }


  extra <- unlist(extra)
  env_proj2 <-
    data.frame(distance = extra, env_proj2)
  rm(extra)


  # Euclidean distance between points used for calibration and its centroid
  base_stand_distance <- env_calib2 %>%
    dplyr::summarise_all(., mean) %>%
    euc_dist(env_calib2, .) %>%
    mean()

  # Standardization of projection points
  env_proj2 <-
    data.frame(
      scaled_distace = env_proj2$distance / base_stand_distance *
        100,
      env_proj2
    )

  # Result
  extraraster[ncell] <- env_proj2[, "scaled_distace"]

  if (!is.null(aggreg_factor)) {
    extraraster <- terra::resample(x = extraraster, y = disag)
    extraraster <- terra::mask(extraraster, disag)
  }

  names(extraraster) <- "extrapolation"
  return(extraraster)
}
