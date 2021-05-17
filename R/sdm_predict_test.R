#' Spatial model prediction and ensemble EXPERIMENTAL
#'
#' @description This function allows the prediction of one or more models built with the fit_ or tune_ function set. It can return continuous or continuous and binary predictions for one or more thresholds
#'
#' @param models list. A list a single or several models fitted with some of fit_ or tune_ functions
#' @param pred stack or brick. Raster layer with predictor variables. Names of raster layers must exactly match those used in model fitting.
#' @param thr character. Binarize projection. Default NULL, i.e. function returns only continuous projection. If used with 'selected_thr', function returns continuous and binarized models used in the 'thr' argument of some of fit_ or tune_ functions. It used with "all_thr" , function returns continuous and binarized for all threshold types.
#' @param calib_area SpatialPolygon or SpatialPolygonDataFrame. Spatial polygon used for restrinc prediction into a given region. Default = NULL
#' @param clamp logical. It is set with TRUE, predictors and features are restricted to the range seen during model training. Only valid for Maxent model (see tune_mx and fit_mx)
#' @param pred_type character. Type of response required available "link", "exponential", "cloglog" and "logistic". Default "cloglog". Only valid for Maxent model (see tune_mx and fit_mx)
#' @param ensemble character. Method used to ensemble different models It is possible to use more than one method. A vector must be provided for this argument. For meansup, meanw or pcasup method it is necessary provide an evaluation metric to ensemble arguments (i.e., AUC, Kappa, TSS, Jaccard, Sorensen or Fpb) see below. (default NULL):
#'   \itemize{
#'   \item mean: Simple average of the different models. Usage ensemble=c(method='mean').
#'   \item meanw: Weighted average of models based on their performance. An evaluation metric and threshold type must be provided. Usage ensemble = c(method='meanw', metric='TSS', thr = 'max_sens_spec').
#'   \item meansup: Average of the best models (e.g., TSS over the average). An evaluation metric must be provided. Usage ensemble=c(method='meansup', metric='TSS').
#'   \item meanthr: pca performed only with those cells with suitability values above the selected threshold. Usage ensemble=c(method='pcathr').
#'   \item median: Median of the different models. Usage ensemble = c(method = 'median')
#'   }
#'
#'  In the case of use more than one ensemble method it is necessary concatenate the names of ensemble methods within the argument, e.g., ensemble=c(method=c('mean', 'median')) or ensemble = c(method=c('mean', 'meanw', 'meansup', 'meanthr', 'median'), metric='TSS', thr = 'max_sens_spec')
#'
#' @return
#'
#' A list of Raster or RasterStack with continuous or continuous and binary predictions
#'
#' @export
#'
#' @importFrom dplyr mutate across left_join bind_rows filter pull
#' @importFrom gam predict.Gam
#' @importFrom GRaF predict.graf
#' @importFrom kernlab predict
#' @importFrom stats predict
#' @importFrom terra crop as.data.frame values rast mean weighted.mean median nlyr ifel
#' @examples
sdm_predict_test <- function(models, pred, thr, calib_area = NULL, clamp = TRUE, pred_type = "cloglog", ensemble = NULL) {

  #### Prepare datasets ####
  # Crop projection area
  if (!is.null(calib_area)) {
    pred2 <- terra::crop(pred, calib_area)
  }


  # Prepare model list
  m <- lapply(models, function(x) x[[1]])
  names(m) <- paste0("m_", 1:length(m))

  # Extract model names object
  clss <- sapply(m, function(x) class(x)[1]) %>%
    tolower() %>%
    gsub(".formula", "", .)

  # Create list for storing raster for current condition
  model_c <- as.list(names(m))
  names(model_c) <- names(m)

  ## Transform raster to data.frame
  pred_df <- terra::as.data.frame(pred, cells = FALSE, na.rm = TRUE)

  #                                                          #
  ####          Prediction for different models           ####
  #                                                          #

  #### graf models ####
  wm <- which(clss == "graf")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      suppressWarnings(model_c[[i]] <- GRaF::predict.graf(m[[i]], pred, type = "response", CI = NULL))
    }
  }


  #### gam models ####
  wm <- which(clss == "gam")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      r <- pred[[1]]
      terra::values(r) <- NA

      # Test factor levels
      f <- which(sapply(m[[i]]$data, class) == "factor")
      if (f > 0) {
        for (ii in 1:length(f)) {
          vf <- m[[i]]$data[, f[ii]] %>% unique()
          vf2 <- pred_df[, names(f[ii])] %>% unique()
          vfilter <- list()
          if (sum(!vf2 %in% vf) > 0) {
            vfilter[[ii]] <- !pred_df[, names(f[ii])] %in% vf
          }
        }
        if (length(vfilter) > 0) {
          if (length(vfilter) > 1) {
            vfilter <- vapply(do.call("rbind", vfilter), any, logical(1))
          } else {
            vfilter <- vfilter[[1]]
          }
        } else {
          vfilter <- 0
        }
      }

      if (sum(vfilter) > 0) {
        v <- rep(0, nrow(pred_df))
        v[!vfilter] <-
          gam::predict.Gam(m[[i]], pred_df[!vfilter, ], type = "response")
        r[as.numeric(rownames(pred_df))] <- v
        rm(v)
      } else {
        r[as.numeric(rownames(pred_df))] <-
          gam::predict.Gam(m[[i]], pred_df, type = "response")
      }

      model_c[[i]] <- r
    }
  }

  #### glm models ####
  wm <- which(clss == "glm")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      r <- pred[[1]]
      terra::values(r) <- NA

      # Test factor levels
      f <- which(sapply(m[[i]]$data, class) == "factor")
      if (f > 0) {
        for (ii in 1:length(f)) {
          vf <- m[[i]]$data[, f[ii]] %>% unique()
          vf2 <- pred_df[, names(f[ii])] %>% unique()
          vfilter <- list()
          if (sum(!vf2 %in% vf) > 0) {
            vfilter[[ii]] <- !pred_df[, names(f[ii])] %in% vf
          }
        }
        if (length(vfilter) > 0) {
          if (length(vfilter) > 1) {
            vfilter <- vapply(do.call("rbind", vfilter), any, logical(1))
          } else {
            vfilter <- vfilter[[1]]
          }
        } else {
          vfilter <- 0
        }
      }

      if (sum(vfilter) > 0) {
        v <- rep(0, nrow(pred_df))
        v[!vfilter] <-
          stats::predict(m[[i]], pred_df[!vfilter, ], type = "response")
        r[as.numeric(rownames(pred_df))] <- v
        rm(v)
      } else {
        r[as.numeric(rownames(pred_df))] <-
          stats::predict(m[[i]], pred_df, type = "response")
      }

      model_c[[i]] <- r
    }
  }

  #### gbm models ####
  wm <- which(clss == "gbm")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      r <- pred[[1]]
      terra::values(r) <- NA
      r[as.numeric(rownames(pred_df))] <-
        suppressMessages(stats::predict(m[[i]], pred_df, type = "response"))

      model_c[[i]] <- r
      rm(i)
    }
  }


  #### maxnet models ####
  wm <- which(clss == "maxnet")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      r <- pred[[1]]
      terra::values(r) <- NA
      r[as.numeric(rownames(pred_df))] <-
        predict_maxnet(object = m[[i]], newdata = pred_df, clamp = clamp, type = pred_type)
      model_c[[i]] <- r
      rm(i)
    }
  }

  #### nnet class ####
  wm <- which(clss == "nnet")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      r <- pred[[1]]
      terra::values(r) <- NA

      # Test factor levels
      f <- (m[[i]]$xlevels)

      if (length(f) > 0) {
        for (ii in 1:length(f)) {
          vf <- f[[ii]] %>%
            unique() %>%
            as.numeric()
          vf2 <- pred_df[, names(f[ii])] %>% unique()
          vfilter <- list()
          if (sum(!vf2 %in% vf) > 0) {
            vfilter[[ii]] <- !pred_df[, names(f[ii])] %in% vf
          }
        }
        if (length(vfilter) > 0) {
          if (length(vfilter) > 1) {
            vfilter <- vapply(do.call("rbind", vfilter), any, logical(1))
          } else {
            vfilter <- vfilter[[1]]
          }
        } else {
          vfilter <- 0
        }
      }

      if (sum(vfilter) > 0) {
        v <- rep(0, nrow(pred_df))
        v[!vfilter] <-
          stats::predict(m[[i]], pred_df[!vfilter, ], type = "raw")
        r[as.numeric(rownames(pred_df))] <- v
        rm(v)
      } else {
        r[as.numeric(rownames(pred_df))] <-
          stats::predict(m[[i]], pred_df, type = "raw")
      }

      model_c[[i]] <- r
    }
  }


  #### randomforest class ####
  wm <- which(clss == "randomforest")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      r <- pred[[1]]
      terra::values(r) <- NA

      # Test factor levels
      f <-
        m[[i]]$forest$xlevels[sapply(m[[i]]$forest$xlevels, function(x) {
          class(x) == "character"
        })]
      if (length(f) > 0) {
        for (ii in 1:length(f)) {
          vf <- f[[ii]] %>%
            unique() %>%
            as.numeric()
          vf2 <- pred_df[, names(f[ii])] %>% unique()
          vfilter <- list()
          if (sum(!vf2 %in% vf) > 0) {
            vfilter[[ii]] <- !pred_df[, names(f[ii])] %in% vf
          }
        }
        if (length(vfilter) > 0) {
          if (length(vfilter) > 1) {
            vfilter <- vapply(do.call("rbind", vfilter), any, logical(1))
          } else {
            vfilter <- vfilter[[1]]
          }
        } else {
          vfilter <- 0
        }
      }


      if (sum(vfilter) > 0) {
        v <- rep(0, nrow(pred_df))
        v[!vfilter] <-
          stats::predict(m[[i]], pred_df[!vfilter, ] %>%
            dplyr::mutate(dplyr::across(
              .cols = names(f),
              .fns = ~ droplevels(.)
            )),
          type = "prob"
          )[, 2]
        r[as.numeric(rownames(pred_df))] <- v
        rm(v)
      } else {
        r[as.numeric(rownames(pred_df))] <-
          stats::predict(m[[i]], pred_df, type = "prob")[, 2]
      }

      model_c[[i]] <- r
    }
  }

  #### ksvmj class ####
  wm <- which(clss == "ksvm")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      r <- pred[[1]]
      terra::values(r) <- NA

      # Test factor levels
      f_n <- which(sapply(pred_df, class) == "factor") %>% names()
      f_n2 <- m[[i]]@xmatrix[[1]] %>% colnames()
      f <- lapply(f_n, function(x) {
        gsub(x, "", grep(x, f_n2, value = TRUE))
      })
      names(f) <- f_n

      if (length(f) > 0) {
        for (ii in 1:length(f)) {
          vf <- f[[ii]] %>%
            unique() %>%
            as.numeric()
          vf2 <- pred_df[, names(f[ii])] %>% unique()
          vfilter <- list()
          if (sum(!vf2 %in% vf) > 0) {
            vfilter[[ii]] <- !pred_df[, names(f[ii])] %in% vf
          }
        }
        if (length(vfilter) > 0) {
          if (length(vfilter) > 1) {
            vfilter <- vapply(do.call("rbind", vfilter), any, logical(1))
          } else {
            vfilter <- vfilter[[1]]
          }
        } else {
          vfilter <- 0
        }
      }

      if (sum(vfilter) > 0) {
        v <- rep(0, nrow(pred_df))
        v[!vfilter] <-
          kernlab::predict(m[[i]], pred_df[!vfilter, ] %>%
            dplyr::mutate(dplyr::across(
              .cols = names(f),
              .fns = ~ droplevels(.)
            )), type = "prob")[, 2]
        r[as.numeric(rownames(pred_df))] <- v
        rm(v)
      } else {
        r[as.numeric(rownames(pred_df))] <-
          kernlab::predict(m[[i]], pred_df, type = "prob")[, 2]
      }
      model_c[[i]] <- r
    }
  }

  df <- data.frame(
    alg = c("gam", "graf", "glm", "gbm", "maxnet", "nnet", "randomforest", "ksvm"),
    names = c("gam", "gau", "glm", "gbm", "mx", "nne", "rf", "svm")
  )

  names(model_c) <- dplyr::left_join(data.frame(alg = clss), df, by = "alg")[, 2]
  model_c <- mapply(function(x, n) {
    names(x) <- n
    x
  }, model_c, names(model_c))


  #### Thresholds ####
  if (!is.null(thr)) {
    if (thr == "selected_thr") {
      thr_df <- lapply(models, function(x) {
        x[["selected_thresholds"]]
      })
    }
    if (thr == "all_thr") {
      thr_df <- lapply(models, function(x) {
        x[["all_thresholds"]]
      })
    }

    model_b <- list()
    for (i in 1:length(model_c)) {
      model_b[[i]] <-
        lapply(thr_df[[i]]$values, function(x) {
          model_c[[i]] >= x
        }) %>% terra::rast()
      names(model_b[[i]]) <- names(thr_df[[i]]$values)
    }
    names(model_b) <- names(model_c)
    mf <- function(x, x2) {
      terra::rast(list(x, x2))
    }
    result <- mapply(mf, model_c, model_b, SIMPLIFY = FALSE)
    return(result)
  } else {
    result <- model_c
    return(result)
  }

  #                                                          #
  ####                      Ensemble                      ####
  #                                                          #
  if (!is.null(ensemble)) {
    ense_type <- ensemble[grep("method", names(ensemble))]

    # Model performance
    perf <- sapply(models, function(x) x["performance"])
    names(perf) <- names(model_c)

    thrm <- sapply(models, function(x) x["selected_thresholds"])
    names(thrm) <- names(model_c)

    if (any(ense_type %in% c("meanw", "meansup", "pcasup"))) {
      perf_2 <-
        lapply(perf, function(x) {
          x[, c("threshold", paste0(ensemble["metric"], "_mean"))]
        }) %>%
        dplyr::bind_rows(.id = "algo")

      thrm_2 <-
        lapply(thrm, function(x) {
          dplyr::filter(x, threshold == ensemble["thr" == names(ensemble)])
        }) %>%
        dplyr::bind_rows(.id = "algo")

      if (!any("thr" == names(ensemble)) || !any("metric" == names(ensemble))) {
        stop("it is necessary provide a performance metric and threshold type for this ensemble method,\n e.g. ensemble = c(method='meanw', metric='TSS', thr = 'max_sens_spec')")
      }

      # WRITE HERE A TEST TO DETECT IF EXIST A METRIC FOR ONE OR MORE USED THRESHOLDS
    }

    # Prepare SpatRaster and data.frame
    p_algo <- as(terra::rast(lapply(model_c, function(x) x[[1]])), "SpatRaster")

    # Mean
    if (any(ense_type == "mean")) {
      en_mean <- terra::mean(p_algo)
      names(en_mean) <- "en_mean"
    }

    # Weighted mean
    if (any(ense_type == "meanw")) {
      w <- perf_2 %>%
        dplyr::filter(threshold == ensemble["thr"]) %>%
        dplyr::pull(3)
      en_meanw <- terra::weighted.mean(p_algo, w = w)
      names(en_meanw) <- "en_meanw"
    }

    # Mean sup
    if (any(ense_type == "meansup")) {
      w <- perf_2 %>%
        dplyr::filter(threshold == ensemble["thr"]) %>%
        dplyr::pull(3)
      en_meansup <- terra::mean(p_algo[[w > mean(w)]])
      names(en_meansup) <- "en_meansup"
    }

    # Median
    if (any(ense_type == "median")) {
      en_median <- terra::median(p_algo)
      names(en_median) <- "en_median"
    }

    # Mean threshold
    if (any(ense_type == "meanthr")) {
      en_meanthr <- p_algo
      for (i in 1:terra::nlyr(p_algo)) {
        en_meanthr[[i]] <- terra::ifel(p_algo[[i]] > thrm_2$values[[i]], p_algo[[i]], 0)
      }
      en_meanthr <- terra::mean(en_meanthr)
      names(en_meanthr) <- "en_meanthr"
    }

    en_object <- c("en_mean", "en_meansup", "en_meanthr", "en_median", "en_meanw") %>%
      sapply(., exists) %>%
      which() %>%
      names()
    result <- terra::rast(sapply(en_object, get))

    # Binarize ensemble (FOR THIS PART OF THE CODE IT IS NECESSARY PROVIDE THRESHOLD VALUES FROM fit_ and tune_ function families)
  } else {
    return(result)
  }
}
