#' Measure model extrapolation
#'
#' @description Measure extrapolation comparing data used for modeling calibration and area for
#' model projection. This function use the approach proposed by Velazco et al., in prep
#'
#' @param env_calib SpatRaster with environmental conditions of the calibration area or the
#' presence and absence points localities used for constructing models
#' @param env_proj SpatRaster with environmental condition used for projecting a model (e.g., a bigger region, other region, or time period)
#' @param n_cores numeric. Number of cores use for parallelization. Default 1
#' @param aggreg_factor positive integer. Aggregation factor expressed as number of cells in each
#'  direction to reduce raster resolution. Use value higher than 1 would be interesting when
#'  measuring extrapolation with raster with a high number of cells. The resolution of output will be
#'  the same as raster object used in 'env_proj' argument. Default 1, i.e., by default, no changes
#'  will be made to the resolution of the environmental variables
#'
#'
#' @return
#' A SpatRaster object with extrapolation values measured in percentage
#'
#' @seealso \code{\link{extra_exclude}}
#'
#' @export
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr summarise_all
#' @importFrom flexclust dist2
#' @importFrom foreach foreach
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
#'
#' spp$species %>% unique()
#' sp <- spp %>%
#'   dplyr::filter(species == "sp3", pr_ab == 1) %>%
#'   dplyr::select(x, y)
#'
#' # Accessible area
#' ca <- calib_area(sp, x = "x", y = "y", method = c("buffer", width = 30000))
#'
#' plot(somevar$CFP_1)
#' points(sp)
#' plot(ca, add = T)
#'
#' # Get environmental condition of calibration area
#' somevar_ca <- somevar %>%
#'   crop(., ca) %>%
#'   mask(., ca)
#' plot(somevar_ca)
#'
#' xp <-
#'   extra_eval(
#'     env_calib = somevar_ca,
#'     env_proj = somevar,
#'     n_cores = 1,
#'     aggreg_factor = 3
#'   )
#' plot(xp)
#' }
extra_eval <- function(env_calib, env_proj, n_cores = 1, aggreg_factor = 1) {
  . <- x <- NULL
  # Get variable names
  v0 <- unique(c(names(env_calib), names(env_proj)))
  v0 <- sort(v0)

  # Test rasters variable names
  if (!all(c(
    all(names(env_calib) %in% names(env_proj)),
    all(names(env_proj) %in% names(env_calib))
  ))) {
    stop(
      "colnames of dataframes of env_records, env_calib, and env_proj
        do not match each other",
      "\nraster layers names:",
      "\n",
      paste(sort(unique(unlist(
        v0
      ))), collapse = "\n")
    )
  }


  # Sort in the same way layer in both raster
  env_calib <- env_calib[[v0]]
  env_proj <- env_proj[[v0]]

  # Layer base
  extraraster <- terra::mask(!is.na(env_proj[[1]]), env_proj[[1]])

  if (aggreg_factor == 1) {
    aggreg_factor <- NULL
  }
  if (!is.null(aggreg_factor)) {
    disag <- extraraster
    env_calib <- terra::aggregate(env_calib, fact = aggreg_factor, na.rm = TRUE)
    env_proj <- terra::aggregate(env_proj, fact = aggreg_factor, na.rm = TRUE)
    extraraster <- terra::aggregate(extraraster, fact = aggreg_factor, na.rm = TRUE)
  }


  # Transform raster to df
  env_calib2 <-
    terra::as.data.frame(env_calib, xy = FALSE, na.rm = TRUE)
  env_proj2 <-
    terra::as.data.frame(env_proj, xy = FALSE, na.rm = TRUE)

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
        flexclust::dist2(
          env_proj2[rowset, v0], # env_proj2 prediction environmental conditions outside
          env_calib2[v0]
        ) # env_calib environmental conditions used as references
      auclidean <- sapply(data.frame(t(auclidean)), min)
      return(auclidean)
    })
  } else {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    extra <- foreach::foreach(x = seq_len((length(
      set
    ) - 1)), .combine = "c") %dopar% {
      rowset <- set[x]:(set[x + 1] - 1)
      auclidean <-
        flexclust::dist2(
          env_proj2[rowset, v0], # env_proj2 prediction environmental conditions outside
          env_calib2[v0]
        ) # env_calib environmental conditions used as references
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
    flexclust::dist2(env_calib2, .) %>%
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
