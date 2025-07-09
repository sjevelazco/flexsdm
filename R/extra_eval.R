#' Measure model extrapolation based on Shape extrapolation metric
#'
#' @description Measure extrapolation comparing environmental data used for modeling calibration
#' and area for model projection. This function use the Shape metric
#' proposed by \href{https://doi.org/10.1111/ecog.06992}{Velazco et al., 2023}
#'
#'
#' @param training_data data.frame or tibble with environmental conditions of
#' presence and absence (or background points or pseudo-absences) used for constructing models
#' @param pr_ab character. Column name with presence and absence (or background points or
#' pseudo-absences) data (i.e., 1 and 0)
#' @param metric character. Metric used to measure degree of extrapolation. Default = mahalanobis.
#' \itemize{
#'   \item mahalanobis: Degree of extrapolation is calculated based on Mahalanobis distance.
#'   \item euclidean: Degree of extrapolation is calculated based on Euclidean distance.
#'   }
#' @param univar_comb logical. If true, the function will add a layer or column to distinguish
#' between univariate (i.e., projection data outside the range of training conditions) and
#' combinatorial extrapolation (i.e., projection data within the range of training conditions)
#' using values 1 and 2, respectively. Default FALSE
#' @param projection_data SpatRaster, data.frame or tibble with environmental condition used for projecting a model (e.g.,
#' a larger, encompassing region, a spatially separate region, or a different time period).
#' If data.frame or tibble is used function will return a tibble object.
#' Otherwise, as SpatRaster object.
#' @param n_cores numeric. Number of cores use for parallelization. Default 1
#' @param aggreg_factor positive integer. Aggregation factor expressed as number of cells in each
#'  direction to reduce raster resolution. Use value higher than 1 would be useful when
#'  measuring extrapolation using a raster with a high number of cells. The resolution of output will be
#'  the same as raster object used in 'projection_data' argument. Default 1, i.e., by default, no changes
#'  will be made to the resolution of the environmental variables.
#'
#'
#' @details This function measure model extrapolation base on the Shape metric
#' (\href{https://doi.org/10.1111/ecog.06992}{Velazco et al., 2023}).
#' Shape is a model-agnostic approach that calculates the extrapolation
#' degree for a given projection data point by its multivariate distance to the nearest training
#' data point. Such distances are relativized by a factor that reflects the dispersion of the
#' training data in environmental space. Distinct from other approaches (e.g.,
#' MESS-Multivariate Environmental Similarity Surfaces, EO-Environmental Overlap,
#' MOP-Mobility-Oriented Parity, EXDET-Extrapolation Detection, or AOA-Area of Applicability),
#' Shape incorporates an adjustable threshold to control the binary discrimination between
#' acceptable and unacceptable extrapolation degrees (see \code{\link{extra_truncate}}).
#'
#' See this \href{https://sjevelazco.github.io/flexsdm/articles/v06_Extrapolation_example.html}{vignette at flexsdm website}
#' for further details about Shape metric, model truncation, and tools to explore model extrapolation.
#'
#' @references
#' \itemize{
#' \item Velazco, S.J.E., Brooke, M.R., De Marco Jr., P., Regan, H.M. and Franklin, J. 2023.
#' How far can I extrapolate my species distribution model? Exploring Shape, a novel method.
#' Ecography: e06992. https://doi.org/10.1111/ecog.06992
#' }
#'
#' @return
#' A SpatRaster or tibble object with extrapolation values measured by Shape metric. Also it
#' is possible estimate univariate and combinatorial extrapolation metric (see `univar_comb` argument).
#'
#' @seealso \code{\link{extra_truncate}}, \code{\link{p_extra}}, \code{\link{p_pdp}}, \code{\link{p_bpdp}}
#'
#' @export
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr summarise_all
#' @importFrom foreach foreach "%dopar%"
#' @importFrom parallel makeCluster stopCluster
#' @importFrom stats sd mahalanobis cov
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
#' ca <- calib_area(sp,
#'   x = "x", y = "y",
#'   method = c("bmcp", width = 50000),
#'   crs = crs(somevar)
#' )
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
#'   n = nrow(sp) * 2,
#'   method = "random",
#'   rlayer = somevar,
#'   calibarea = ca
#' )
#'
#' # Merge presences and absences databases to get a complete calibration data
#' sp_pa <- dplyr::bind_rows(sp, psa)
#' sp_pa
#'
#' # Get environmental condition of calibration area
#' sp_pa_2 <- sdm_extract(
#'   data = sp_pa,
#'   x = "x",
#'   y = "y",
#'   env_layer = somevar
#' )
#' sp_pa_2
#'
#' # Measure degree of extrapolation based on Mahalanobis and
#' # for a projection area based on a SpatRaster object
#' extr <-
#'   extra_eval(
#'     training_data = sp_pa_2,
#'     projection_data = somevar,
#'     pr_ab = "pr_ab",
#'     n_cores = 1,
#'     aggreg_factor = 1,
#'     metric = "mahalanobis"
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
#' predsuit <- sdm_predict(
#'   models = a_model,
#'   pred = somevar,
#'   thr = "max_sorensen"
#' )
#' predsuit # list with a raster with two layer
#' plot(predsuit[[1]])
#'
#' # Truncate a model based on a given value of extrapolation
#' # using 'extra_truncate' function
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
#' predsuit_2 # a list of continuous and binary models with
#' # different truncated at different extrapolation thresholds
#'
#' plot(predsuit_2$`50`)
#' plot(predsuit_2$`100`)
#' plot(predsuit_2$`200`)
#'
#'
#' ## %######################################################%##
#' ####        Measure degree of extrapolation for         ####
#' ####        projection area based on data.frame         ####
#' ## %######################################################%##
#'
#' extr_df <-
#'   extra_eval(
#'     training_data = sp_pa_2,
#'     projection_data = as.data.frame(somevar, xy = TRUE),
#'     pr_ab = "pr_ab",
#'     n_cores = 1,
#'     aggreg_factor = 1,
#'     metric = "mahalanobis"
#'   )
#' extr_df
#' # see 'p_extra()' to explore extrapolation or suitability pattern in the
#' # environmental and/or geographical space
#'
#' ## %######################################################%##
#' ####             Explore Shape metric with              ####
#' ####     univariate and combinatorial extrapolation     ####
#' ## %######################################################%##
#' extr <-
#'   extra_eval(
#'     training_data = sp_pa_2,
#'     projection_data = somevar,
#'     pr_ab = "pr_ab",
#'     n_cores = 1,
#'     aggreg_factor = 1,
#'     metric = "mahalanobis",
#'     univar_comb = TRUE
#'   )
#'
#' extr
#' plot(extr) # In the second layer, values equal to 1 and 2
#' # depict univariate and combinatorial extrapolation, respectively
#' }
extra_eval <-
  function(training_data,
           pr_ab,
           projection_data,
           metric = "mahalanobis",
           univar_comb = FALSE,
           n_cores = 1,
           aggreg_factor = 1) {
    Value <- val <- . <- x <- extrapolation <- NULL

    if (!metric %in% c("euclidean", "mahalanobis")) {
      stop("metric argument must be used with 'euclidean' or 'mahalanobis'")
    }
    if (length(metric) > 1) {
      stop("metric argument must be used with 'euclidean or 'mahalanobis'")
    }


    if (any("data.frame" == class(training_data))) {
      training_data_pr_ab <- training_data[c(names(projection_data), pr_ab)] %>%
        na.omit()
      training_data <- training_data_pr_ab[names(projection_data)] %>%
        na.omit()
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

    if (any("SpatRaster" == class(projection_data))) {
      projection_data <- projection_data[[v0]]
      # Layer base
      extraraster <- projection_data[[1]]
      extraraster[!is.na(extraraster)] <- 0

      # Change resolution of raster
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
    } else {
      projection_data <- projection_data[v0]
    }



    # Transform raster to df
    env_calib2 <-
      terra::as.data.frame(training_data, xy = FALSE, na.rm = TRUE)
    if("SpatRaster" == class(projection_data)){
      env_proj2 <-
        terra::as.data.frame(projection_data, xy = FALSE, na.rm = TRUE)
    } else {
      env_proj2 <- projection_data
      # names(env_proj2) <- names(projection_data)
    }


    # save coordinates and cell number
    if (any("SpatRaster" == class(projection_data))) {
      ncell <- rownames(env_proj2) %>% as.numeric()
    }

    # Standardization
    # standardization based on training data
    s_center <- colMeans(env_calib2)
    s_scale <- apply(env_calib2, 2, stats::sd)

    for (i in 1:ncol(env_calib2)) {
      env_calib2[i] <- (env_calib2[i] - s_center[i]) / s_scale[i]
      # if (metric == "mahalanobis_pres") {
      #   training_data_pr_ab[i] <-
      #     (training_data_pr_ab[i] - s_center[i]) / s_scale[i]
      # }
    }

    # Calculate covariance matrix based on presences for mahalanobis_pres
    # if (metric == "mahalanobis_pres") {
    #   s_cov <-  stats::cov(training_data_pr_ab[training_data_pr_ab[pr_ab] == 1, names(env_calib2)])
    # }

    for (i in 1:ncol(env_proj2)) {
      env_proj2[i] <- (env_proj2[i] - s_center[i]) / s_scale[i]
    }

    # Measure extrapolation - Euclidean distance
    set <- c(seq(1, nrow(env_proj2), 200), nrow(env_proj2) + 1)

    if (n_cores == 1) {
      extra <- lapply(seq_len((length(set) - 1)), function(x) {
        rowset <- set[x]:(set[x + 1] - 1)
        if (metric == "euclidean") {
          envdist <-
            euc_dist(
              env_proj2[rowset, v0], # env_proj2 environmental conditions used to predict
              env_calib2[v0]
            ) # training_data environmental conditions used as references
          envdist <- sapply(data.frame(t(envdist)), min)
        }
        if (metric == "mahalanobis") {
          envdist <-
            mah_dist(
              x = env_proj2[rowset, v0], # env_proj2 environmental conditions used to predict
              y = env_calib2[v0], # training_data environmental conditions used as references
              cov = stats::cov(env_calib2) # covariance matrix based on presences and absences
            )
          envdist <- sapply(data.frame(t(envdist)), min)
        }
        # if (metric == "mahalanobis_pres") {
        #   envdist <-
        #     mah_dist(
        #       x = env_proj2[rowset, v0], # env_proj2 environmental conditions used to predict
        #       y = env_calib2[v0], # training_data environmental conditions used as references
        #       cov = s_cov # covariance matrix based on presences
        #     )
        #   envdist <- sapply(data.frame(t(envdist)), min)
        # }
        return(envdist)
      })
    } else {
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)

      extra <- foreach::foreach(x = seq_len((length(
        set
      ) - 1)), .export = c("euc_dist"), .combine = "c") %dopar% {
        rowset <- set[x]:(set[x + 1] - 1)
        if (metric == "euclidean") {
          envdist <-
            euc_dist(
              env_proj2[rowset, v0], # env_proj2 environmental conditions used to predict
              env_calib2[v0]
            ) # training_data environmental conditions used as references
          envdist <- sapply(data.frame(t(envdist)), min)
        }
        if (metric == "mahalanobis") {
          envdist <-
            mah_dist(
              x = env_proj2[rowset, v0], # env_proj2 environmental conditions used to predict
              y = env_calib2[v0], # training_data environmental conditions used as references
              cov =  stats::cov(env_calib2) # covariance matrix based on presences and absences
            )
          envdist <- sapply(data.frame(t(envdist)), min)
        }
        # if (metric == "mahalanobis_pres") {
        #   envdist <-
        #     mah_dist(
        #       x = env_proj2[rowset, v0], # env_proj2 environmental conditions used to predict
        #       y = env_calib2[v0], # training_data environmental conditions used as references
        #       cov = s_cov # covariance matrix based on presences
        #     )
        #   envdist <- sapply(data.frame(t(envdist)), min)
        # }
        envdist
      }
      parallel::stopCluster(cl)
    }


    extra <- unlist(extra)
    if (any("SpatRaster" == class(projection_data))) {
      env_proj2 <-
        data.frame(distance = extra)
    } else {
      env_proj2 <-
        cbind(data.frame(distance = extra), env_proj2)
    }
    rm(extra)


    # Euclidean distance between points used for calibration and its centroid
    if (metric == "euclidean") {
      base_stand_distance <- env_calib2 %>%
        dplyr::summarise_all(., mean) %>%
        euc_dist(env_calib2, .) %>%
        mean()
    }
    # if (metric == "mahalanobis_pres") {
    #   base_stand_distance <- env_calib2 %>%
    #     dplyr::summarise_all(., mean) %>%
    #     mah_dist(x = env_calib2, y = ., cov = s_cov) %>%
    #     mean()
    # }
    if (metric == "mahalanobis") {
      base_stand_distance <- env_calib2 %>%
        dplyr::summarise_all(., mean) %>%
        mah_dist(x = env_calib2, y = ., cov = stats::cov(env_calib2)) %>%
        mean()
    }

    # Standardization of projection points
    env_proj2 <-
      data.frame(extrapolation = env_proj2$distance / base_stand_distance * 100) %>%
      cbind(., env_proj2) %>%
      dplyr::select(-distance)

    # Result
    if (any("SpatRaster" == class(projection_data))) {
      extraraster[ncell] <- env_proj2[, "extrapolation"]
      if (!is.null(aggreg_factor)) {
        extraraster <- terra::resample(x = extraraster, y = disag)
        extraraster <- terra::mask(extraraster, disag)
      }
      names(extraraster) <- "extrapolation"

      # Univariate and combinatorial extrapolation
      if (univar_comb) {
        rng <- apply(training_data, 2, range, na.rm = TRUE)
        univar_ext <- projection_data
        for (i in 1:terra::nlyr(projection_data)) {
          univar_ext[[i]] <- (projection_data[v0[i]] <= rng[, v0[i]][1] |
            projection_data[v0[i]] >= rng[, v0[i]][2])
        }
        univar_comb_r <- sum(univar_ext)
        univar_comb_r[univar_comb_r > 0] <- 2
        univar_comb_r[univar_comb_r == 0] <- 1
        names(univar_comb_r) <- "uni_comb"
        extraraster <- c(extraraster, univar_comb_r)
      }
    } else {
      for (i in 2:length(s_center)) { # process from the 2nd column to skip extrapolation column
        env_proj2[i] <- env_proj2[i] * s_scale[i] + s_center[i]
      }
      # Univariate and combinatorial extrapolation
      if (univar_comb) {
        rng <- apply(training_data, 2, range, na.rm = TRUE)
        univar_ext <- projection_data
        for (i in 1:ncol(projection_data)) {
          univar_ext[, v0[i]] <- (projection_data[, v0[i]] <= rng[, v0[i]][1] |
            projection_data[, v0[i]] >= rng[, v0[i]][2])
        }
        univar_comb_r <- apply(univar_ext, 1, sum)
        univar_comb_r[univar_comb_r > 0] <- 2
        univar_comb_r[univar_comb_r == 0] <- 1
        env_proj2 <- env_proj2 %>%
          dplyr::mutate(univar_comb_r) %>%
          dplyr::relocate(extrapolation, univar_comb = univar_comb_r)
      }
      return(dplyr::as_tibble(env_proj2))
    }
    return(extraraster)
  }
