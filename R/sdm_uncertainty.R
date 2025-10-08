#' Calculate species distribution model uncertainty
#'
#' @description This function calculates the uncertainty of a species distribution model by performing a bootstrap procedure.
#' It refits the model multiple times on resampled data and then calculates the standard deviation of the predictions across
#' all iterations.
#'
#' @param models A flexsdm model object from `fit_*` and `tune_*` functions.
#' @param training_data A data.frame or tibble with presence-absence/background data and predictors.
#' @param background A data.frame or tibble with background data, used for Maxent and Gausian Process models. Default NULL.
#' @param response character. Column name of the response variable.
#' @param projection_data A SpatRaster object with the environmental layers for projection.
#' @param iteration numeric. The number of bootstrap iterations. A large number of iterations will increase computation time, especially for large areas or high-resolution data. Default 50.
#' @param n_cores numeric. The number of cores to use for parallel processing. Default 5.
#' @param clamp logical. If TRUE, predictors and features are restricted to the range seen during model training. Only for Maxent models. Default TRUE.
#' @param pred_type character. The type of response required. Available options are: "link", "exponential", "cloglog", and "logistic". Only for Maxent models. Default "cloglog".
#'
#' @return A SpatRaster object with a single layer representing the model uncertainty, calculated as the standard deviation of the bootstrap predictions.
#' @export
#' @seealso \code{\link{sdm_predict}}
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#' 
#' data(spp)
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#' 
#' sp_data <- spp %>% filter(species == "sp3")
#' sp_data <- part_random(
#'   data = sp_data,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 3)
#' )
#' 
#' sp_data <- sdm_extract(
#'     data = sp_data,
#'     x = "x",
#'     y = "y",
#'     env_layer = somevar
#'   )
#' 
#' m <- fit_svm(
#'   data = sp_data,
#'   response = "pr_ab",
#'   predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"), 
#'     partition = ".part"
#' )
#' 
#' unc <- sdm_uncertainty(
#'   models = m,
#'   training_data = sp_data,
#'   response = "pr_ab",
#'   projection_data = somevar,
#'   iteration = 10,
#'   n_cores = 2
#' )
#' 
#' plot(unc)
#' }
#' @export
sdm_uncertainty <- function(
    models,
    training_data,
    background = NULL,
    response,
    projection_data,
    iteration = 50,
    n_cores = 5,
    clamp = TRUE,
    pred_type = "cloglog") {
  # Extract model names object
  if ("domain" %in% names(models[[1]])) {
    clss <- "domain"
  } else {
    clss <- class(models[[1]])[1] %>%
      tolower() %>%
      gsub(".formula", "", .)
  }


  # Prepare training dataset
  training_data <- training_data %>% dplyr::select(-starts_with(".part"))

  # Prepare raster
  proj_data_df <- as.data.frame(projection_data)
  proj_data_df <-  proj_data_df[complete.cases(proj_data_df),]
  r <- projection_data[[!terra::is.factor(projection_data)]][[1]]
  r <- terra::ifel(!is.na(r), 0, NA)

  # Prepare predictors
  pr_c <- models$predictors %>%
    dplyr::select(dplyr::starts_with("c")) %>%
    unlist()
  pr_f <- models$predictors %>%
    dplyr::select(dplyr::starts_with("f")) %>%
    unlist()
  names(pr_c) <- NULL
  names(pr_f) <- NULL


  #### Bootstrap approach ####
  if (clss == "domain") {
    r_list <- NULL
    for (ii in 1:iteration) {
      set.seed(ii)
      db <- part_random(
        data = training_data,
        pr_ab = response,
        method = c(method = "boot", replicates = "1", proportion = "0.8")
      ) %>% dplyr::filter(grepl("train", .part1))

      m <- switch(clss,
        "domain" = fit_dom(data = db, response = response, predictors = pr_c, predictors_f = pr_f, partition = NULL)$model
      )
      r_list[[ii]] <- switch(clss,
        "domain" = 1-map_env_dist(training_data = m$domain, projection_data = proj_data_df, metric = "domain")
      )
    }
  } else {
    my_cluster <- makeCluster(n_cores)
    registerDoParallel(my_cluster)

    r_list <- foreach::foreach(
      ii = 1:iteration,
      .packages = c("dplyr", "mgcv", "kernlab"), # Packages needed on each worker
      .export = c("part_random", "predict.graf", "predict_maxnet", "extract.maxnet.classes", "fit_gam", "fit_gau", "fit_glm", "fit_gbm", "fit_dom", "fit_max", "fit_net", "fit_raf", "fit_svm", "map_env_dist"),
      .errorhandling = "pass" # 
    ) %dopar% {
      set.seed(ii)
      db <- part_random(
        data = training_data,
        pr_ab = response,
        method = c(method = "boot", replicates = "1", proportion = "0.8")
      ) %>% dplyr::filter(grepl("train", .part1))

      m <- switch(clss,
        "gam" = fit_gam(data = db, response = response, predictors = pr_c, predictors_f = pr_f, fit_formula = models$model$formula, partition = NULL)$model,
        "glm" = fit_glm(data = db, response = response, predictors = pr_c, predictors_f = pr_f, fit_formula = models$model$formula, partition = NULL)$model,
        "gbm" = fit_gbm(data = db, response = response, predictors = pr_c, predictors_f = pr_f, partition = NULL, n_trees = models$performance$n.trees, n_minobsinnode = models$performance$n.minobsinnode, shrinkage = models$performance$shrinkage)$model,
        "graf" = fit_gau(data = db, response = response, predictors = pr_c, predictors_f = pr_f, background = background, partition = NULL)$model,
        "nnet" = fit_net(data = db, response = response, predictors = pr_c, predictors_f = pr_f, fit_formula = NULL, partition = NULL, size = models$model$n[2], decay = models$model$decay)$model,
        "randomforest" = fit_raf(data = db, response = response, predictors = pr_c, predictors_f = pr_f, fit_formula = NULL, partition = NULL, mtry = models$model$mtry, ntree = models$model$ntree)$model,
        "ksvm" = fit_svm(data = db, response = response, predictors = pr_c, predictors_f = pr_f, fit_formula = NULL, partition = NULL, sigma = models$model@kernelf@kpar$sigma, C = models$model@param$C)$model,
        "maxnet" = {
          terms <- models$model$beta %>% rownames()
          nums <- lapply(strsplit(terms, ":"), function(x) x[2]) %>% unlist()
          suppressWarnings(is_num <- !is.na(as.numeric(nums)))
          terms <- c(terms[!is_num], sapply(strsplit(terms[is_num], ":"), function(x) x[1]) %>% unique()) %>%
            paste(collapse = " + ")
          terms <- paste("~", terms) %>% as.formula()
          terms <- extract.maxnet.classes(terms)
          fit_max(data = db, response = response, predictors = pr_c, predictors_f = pr_f, fit_formula = NULL, partition = NULL, background = background, clamp = clamp, pred_type = pred_type, classes = terms, regmult = 1)$model
        }
      )

      result <- switch(clss,
        "gam" = mgcv::predict.gam(m, proj_data_df, type = "response", se.fit = FALSE),
        "glm" = stats::predict(m, proj_data_df, type = "response"),
        "gbm" = stats::predict(m, proj_data_df, type = "response"),
        "graf" = predict.graf(object = m, newdata = proj_data_df[, names(m$obsx)], type = "response", CI = NULL)[, 1],
        "nnet" = stats::predict(m, proj_data_df, type = "raw")[, 1],
        "randomforest" = stats::predict(m, proj_data_df, type = "prob")[, 2],
        "ksvm" = kernlab::predict(m, newdata = proj_data_df, type = "prob")[, 2],
        "maxnet" = predict_maxnet(m, newdata = proj_data_df, clamp = clamp, type = pred_type)[, 1]
      )

      result
    }
    stopCluster(my_cluster)
  }

  names(r_list) <- 1:length(r_list)
  r_list <- data.frame(r_list) %>% apply(1, sd, na.rm = TRUE)
  r[rownames(proj_data_df) %>% as.numeric()] <- r_list
  rm(r_list)

  names(r) <- "uncertainty"
  return(r)
}
