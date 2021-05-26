#' Spatial predictions of individual and ensemble model
#'
#' @description This function allows the prediction of one or more models built with the fit_ or tune_ function set. It can return continuous or continuous and binary predictions for one or more thresholds
#'
#' @param models list of individual models or a fit_ensemble output. A list a single or several models fitted with some of fit_ or tune_ functions or object returned by the \code{\link{fit_ensemble}} function.
#' @param pred SpatRaster. Raster layer with predictor variables. Names of layers must exactly match those used in model fitting.
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is possible to use more than one threshold type. It is mandatory to use the same threshold/s used to fit the models. The next threshold area available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which Jaccard is the highest.
#'   \item max_sorensen: The threshold at which Sorensen is highest.
#'   \item max_fpb: The threshold at which FPB is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value used to fit the models.
#'   \item all: Will be used all the threshold filled in the 'models' argument. It is not possible to combine with previews one
#'   }
#' Usage thr = c('lpt', 'max_sens_spec', 'max_jaccard'), thr=c('lpt', 'max_sens_spec', 'sensitivity'), or thr='all'. If no threshold is specified (i.e., thr = NULL) function will return continuous prediction only. Default NULL
#' @param con_thr logical. If true will be returned predictions with suitability values above threshold/s. Default = FALSE
#' @param predict_area SpatVector, SpatialPolygon, or SpatialPolygonDataFrame. Spatial polygon used for restring prediction into a given region. Default = NULL
#' @param clamp logical. It is set with TRUE, predictors and features are restricted to the range seen during model training. Only valid for Maxent model (see tune_mx and fit_mx)
#' @param pred_type character. Type of response required available "link", "exponential", "cloglog" and "logistic". Default "cloglog". Only valid for Maxent model (see tune_mx and fit_mx)
#'
#' @return
#'
#' A list of SpatRaster with continuous and/or binary predictions
#'
#' @export
#'
#' @seealso \code{\link{fit_ensemble}}
#'
#' @importFrom dplyr mutate across left_join pull bind_rows filter all_of select
#' @importFrom gam predict.Gam
#' @importFrom GRaF predict.graf
#' @importFrom kernlab predict
#' @importFrom stats predict median
#' @importFrom terra vect crop mask as.data.frame values rast app lapp
#' @examples
sdm_predict <-
  function(models,
           pred,
           thr = NULL,
           con_thr = FALSE,
           predict_area = NULL,
           clamp = TRUE,
           pred_type = "cloglog") {
    if (all(names(ens) %in% c("models", "thr_metric", "predictors", "performance"))) {
      message("Predicting ensembles")
      ensembles <- models
      models <- NULL
    } else {
      message("Predicting individual models")
      ensembles <- NULL
    }

    #### Prepare datasets ####
    # Crop and mask projection area
    if (!is.null(predict_area)) {
      if (class(ca_1) %in% c("SpatialPolygons", "SpatialPolygonsDataFrame")) {
        predict_area <- terra::vect(predict_area)
      }
      pred <-
        pred %>%
        terra::crop(., predict_area) %>%
        terra::mask(., predict_area)
    }

    # Transform raster to data.frame
    pred_df <-
      terra::as.data.frame(pred, cells = FALSE, na.rm = TRUE)


    #### Model predictions
    if (!is.null(models)) {
      # Prepare model list
      m <- lapply(models, function(x) {
        x[[1]]
      })
      names(m) <- paste0("m_", 1:length(m))

      # Extract model names object
      clss <- sapply(m, function(x) {
        class(x)[1]
      }) %>%
        tolower() %>%
        gsub(".formula", "", .)
    }

    if (!is.null(ensembles)) {
      # Prepare model list
      m <- lapply(ensembles$models, function(x) x[["model"]])
      names(m) <- paste0("m_", 1:length(m))

      # Extract model names object
      clss <- sapply(m, function(x) {
        class(x)[1]
      }) %>%
        tolower() %>%
        gsub(".formula", "", .)
    }


    ## %######################################################%##
    #                                                          #
    ####          Prediction for different models           ####
    #                                                          #
    ## %######################################################%##

    # Create list for storing raster for current condition
    model_c <- as.list(names(m))
    names(model_c) <- names(m)

    #### graf models ####
    wm <- which(clss == "graf")
    if (length(wm) > 0) {
      wm <- names(wm)
      for (i in wm) {
        suppressWarnings(model_c[[i]] <-
          GRaF::predict.graf(m[[i]], pred, type = "response", CI = NULL))
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
          predict_maxnet(
            object = m[[i]],
            newdata = pred_df,
            clamp = clamp,
            type = pred_type
          )
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


    rm(pred_df)

    df <- data.frame(
      alg = c(
        "gam",
        "graf",
        "glm",
        "gbm",
        "maxnet",
        "nnet",
        "randomforest",
        "ksvm"
      ),
      names = c("gam", "gau", "glm", "gbm", "max", "net", "raf", "svm")
    )

    names(model_c) <-
      dplyr::left_join(data.frame(alg = clss), df, by = "alg")[, 2]
    model_c <- mapply(function(x, n) {
      names(x) <- n
      x
    }, model_c, names(model_c))



    # Predict ensemble
    if (!is.null(ensembles)) {
      model_c <- terra::rast(model_c) # stack individual models
      ens_perf <- ensembles[["performance"]] # get performance of ensembles
      ens_method <- ens_perf %>%
        dplyr::pull(model) %>%
        unique() # get ensemble methods
      single_model_perf <- lapply(ensembles[[1]], function(x) { # Get performance of individual model used for performing ensemble
        x[["performance"]]
      }) %>% dplyr::bind_rows()

      # Threshold and metric values for performing some ensembles
      if (any(ens_method %in% c("meanw", "meansup", "meanthr"))) {
        weight_data <- single_model_perf %>%
          dplyr::filter(threshold == dplyr::all_of(ensembles$thr_metric[1])) %>%
          dplyr::select(model, thr_value, dplyr::all_of(ensembles$thr_metric[2]))
      }

      ensemble_c <- as.list(ens_method)
      names(ensemble_c) <- ensemble_c

      if (any("mean" == ens_method)) {
        ensemble_c[["mean"]] <- terra::app(model_c, fun = mean, cores = 1)
      }

      if (any("meanw" == ens_method)) {
        ensemble_c[["meanw"]] <- terra::app(model_c * weight_data[[3]], fun = mean, cores = 1)
      }

      if (any("meansup" == ens_method)) {
        ensemble_c[["meansup"]] <-
          terra::app(model_c[[which(weight_data[[3]] >= mean(weight_data[[3]]))]],
            fun = mean, cores = 1
          )
      }

      if (any("meanthr" == ens_method)) {
        mf <- function(x, y) {
          terra::lapp(x, function(x) ifelse(x >= y, x, 0))
        }
        model_c2 <- mapply(mf, model_c, weight_data[[2]], SIMPLIFY = FALSE) %>% terra::rast()
        ensemble_c[["meanthr"]] <- terra::app(model_c2, fun = mean, cores = 1)
        rm(model_c2)
      }

      if (any("median" == ens_method)) {
        ensemble_c[["median"]] <- terra::app(model_c, fun = stats::median, cores = 1)
      }

      ensemble_c <- mapply(function(x, y) {
        names(x) <- y
        x
      }, ensemble_c, ens_method, SIMPLIFY = FALSE)

      model_c <- ensemble_c
      models <- split(ensembles[["performance"]], ensembles[["performance"]]$model)
      models <- lapply(models, function(x) {
        x <- list(x)
        names(x) <- "performance"
        x
      })
      rm(ensembles)
      rm(ensemble_c)
    }



    #### Get binary predictions ####
    if (is.null(thr)) {
      return(model_c)
    } else {
      if (any("all" == thr)) {
        thr_df <- lapply(models, function(x) {
          x[["performance"]]
        })
      } else {
        thr_df <- lapply(models, function(x) {
          x[["performance"]] %>%
            dplyr::filter(threshold %in% thr)
        })
      }

      model_b <- list()

      for (i in 1:length(model_c)) {
        model_b[[i]] <-
          lapply(thr_df[[i]]$thr_value, function(x) {
            model_c[[i]] >= x
          }) %>% terra::rast()
        names(model_b[[i]]) <- names(thr_df[[i]]$thr_value)
      }

      names(model_b) <- names(model_c)

      # Return suitability values above thresholds
      if (con_thr) {
        for (i in 1:length(model_b)) {
          nms <- names(model_b[[i]])
          model_b[[i]] <- model_b[[i]] * model_c[[i]]
          names(model_b[[i]]) <- nms
        }
      }

      mf <- function(x, x2) {
        terra::rast(list(x, x2))
      }

      result <- mapply(mf, model_c, model_b, SIMPLIFY = FALSE)
      return(result)
    }
  }
