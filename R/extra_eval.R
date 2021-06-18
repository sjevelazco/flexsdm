#' Measure model extrapolation
#'
#' @param env_calib SpatRaster with environmental conditions of the calibration area or the
#' presence and absence points localities used for constructing models
#' @param env_proj SpatRaster with environmental condition used for projecting a model
#'
#' @return A SpatRaster object with extrapolation values measured in percentage
#'
#' @seealso \code{\link{extra_correct}}
#'
#' @importFrom dplyr summarise_all
#' @importFrom flexclust dist2
#' @importFrom terra mask as.data.frame
#'
#' @export
#' @examples
extra_eval <- function(env_calib, env_proj) {
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
        round(nrow(env_proj2) / 1000)
      ),
      nrow(env_proj2) + 1
    )

  # TODO write code for processin this part in parallel
  extra <- lapply(seq_len((length(set) - 1)), function(x) {
    rowset <- set[x]:(set[x + 1] - 1)
    auclidean <-
      flexclust::dist2(
        env_proj2[rowset, v0], # env_proj2 prediction environmental conditions outside
        env_calib2[v0]
      ) # env_calib environmental conditions used as references
    auclidean <- apply(auclidean, 1, min)
    return(auclidean)
  })

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

  # Result -----
  extraraster[ncell] <- env_proj2[, "scaled_distace"]

  return(extraraster)
}
