#' Calculate permutation-based variable importance scores for SDMs
#'
#' @description This function calculates permutation-based variable importance scores for species distribution models (SDMs) based on the
#'
#' @param models list of one or more models fitted with fit_ or tune_ functions. In case use models fitted with fit_ensemble or esm_ family function only one model could be used. Usage models = mglm or models = list(mglm, mraf, mgbm)
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of
#' predictor variables.
#' Usage predictors = c("aet", "cwd", "tmin")
#' @param n_sim integer. The number of Monte Carlo replications to perform. Default is 50. The results from each replication are averaged together (the standard deviation will also be returned).
#' @param n_cores numeric. Number of cores use for parallelization. Default 1
#' @param thr character. Threshold criterion used to get binary suitability values (i.e. 0,1).
#' Used for threshold-dependent performance metrics.
#' It is possible to use more than one threshold type.
#' A vector must be provided for this argument. The following threshold criteria are available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the Sensitivity and Specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the Sensitivity and Specificity
#'   is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which the Jaccard index is the highest.
#'   \item max_sorensen: The threshold at which the Sorensen index is the highest.
#'   \item max_fpb: The threshold at which FPB (F-measure on presence-background data) is the highest.
#'   \item sensitivity: Threshold based on a specified Sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers
#'   to Sensitivity value. If a sensitivity value is not specified, the
#'    default value is 0.9
#'   }
#'   #' If more than one threshold type is used, concatenate threshold types,
#' e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec',
#' 'sensitivity', sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity').
#' Function will use all thresholds if no threshold type is specified
#'
#' @param clamp logical. If TRUE, predictors and features are restricted to the range seen during
#'  model training.
#' @param pred_type character. Type of response required available "link", "exponential", "cloglog"
#' and "logistic". Default "cloglog"
#'
#' @return a tibble with the columns:
#' \itemize{
#' \item model: model name
#' \item threshold: threshold names
#' \item predictors: predictor names
#' \item from TPR to IMAE: performance metrics
#' }
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom maxnet maxnet maxnet.formula
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr select arrange mutate group_by summarise_all bind_rows left_join
#' @importFrom foreach foreach
#' @importFrom gbm predict.gbm
#' @importFrom kernlab predict
#' @importFrom mgcv predict.gam
#' @importFrom parallel makeCluster stopCluster
#' @importFrom stats predict.glm predict
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(tidyr)
#' require(dplyr)
#'
#' data(abies)
#' abies
#'
#' data(backg)
#' backg
#'
#' # In this example we will partition the data using the k-fold method
#'
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 5)
#' )
#'
#' backg2 <- part_random(
#'   data = backg,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 5)
#' )
#'
#' max_t1 <- fit_max(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   background = backg2,
#'   thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
#'   clamp = TRUE,
#'   classes = "default",
#'   pred_type = "cloglog",
#'   regmult = 1
#' )
#'
#' net_t1 <- fit_net(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
#' )
#'
#' svm_f1 <- fit_svm(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
#' )
#'
#' vip_t <- sdm_varimp(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth", "landform"),
#'   models = list(max_t1, net_t1, svm_f1),
#'   clamp = TRUE,
#'   pred_type = "cloglog",
#'   thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
#'   n_sim = 50,
#'   n_cores = 5
#' )
#'
#' vip_t
#'
#' # Plot the variable importance for AUC TSS and SORENSEN for the
#' # threshold that maximizes Sorensen metric and Maxent
#' vip_t %>%
#'   pivot_longer(
#'     cols = TPR:IMAE,
#'     names_to = "metric",
#'     values_to = "value"
#'   ) %>%
#'   dplyr::filter(threshold == "max_sorensen") %>%
#'   dplyr::filter(metric %in% c("AUC", "TSS", "SORENSEN")) %>%
#'   dplyr::filter(model == "max") %>%
#'   ggplot(aes(x = reorder(predictors, value), y = value, fill = predictors)) +
#'   geom_col(
#'     col = "black",
#'     show.legend = FALSE
#'   ) +
#'   facet_wrap(~metric, scales = "free_x") +
#'   labs(x = "Predictors", y = "Variable Importance") +
#'   theme_classic() +
#'   coord_flip()
#' }
#'
sdm_varimp <- function(models,
                       data,
                       response,
                       predictors,
                       n_sim = 50,
                       n_cores = 1,
                       thr = NULL,
                       clamp = TRUE,
                       pred_type = "cloglog") {
  . <- model <- threshold <- TPR <- IMAE <- TSS <- BOYCE <- NULL

  # Function to calculate performance
  var_perf <- function(x) {
    p <- sdm_eval(
      p = x$pred[x$pr_ab == 1],
      a = x$pred[x$pr_ab == 0],
      thr = thr
    )


    p %>%
      as.data.frame() %>%
      dplyr::select(threshold, TPR:IMAE) %>%
      dplyr::arrange(threshold) %>%
      dplyr::mutate(TSS = TSS + 1, BOYCE = BOYCE + 1)
  }
  # Function to calculate variable importance
  var_prep <- function(x, x2) {
    x <- as.data.frame(x)
    x <- suppressMessages(
      x %>%
        dplyr::group_by(threshold) %>%
        dplyr::select(-threshold) %>%
        dplyr::summarise_all(mean) %>%
        dplyr::arrange(threshold)
    )
    perf_n <- names(x[-1])

    for (k in perf_n) {
      if (k == "OR") {
        x2[k] <- abs(x2[k] - x[k])
      } else {
        x2[k] <- x2[k] - x[k]
      }
    }
    return(x2)
  }

  # if (is.null(names(models))) {
  #   ensembles <- NULL
  #   esm <- NULL
  # } else if (all(names(models) %in% c("esm_model", "predictors", "performance"))) {
  #   esm <- models
  #   models <- NULL
  #   ensembles <- NULL
  # }

  if (is.null(names(models))) {
    message("Calculating variable importance for a list of individual models")
    ensembles <- NULL
    esm <- NULL
  } else if (all(names(models) %in% c("models", "thr_metric", "predictors", "performance"))) {
    message("Calculating variable importance for ensembles")
    ensembles <- models
    models <- NULL
    esm <- NULL
  } else if (all(names(models) %in% c("esm_model", "predictors", "performance"))) {
    message("Calculating variable importance for ensemble of small models")
    esm <- models
    models <- NULL
    ensembles <- NULL
  } else {
    message("Calculating variable importance for individual models")
    models <- list(models)
    ensembles <- NULL
    esm <- NULL
  }

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

  #### Ensemble predictions
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

  #### ESM predictions
  if (!is.null(esm)) {
    # Prepare model list
    m <- esm$esm_model
    names(m) <- paste0("m_", 1:length(m))

    # Extract model names object
    clss <- sapply(m, function(x) {
      class(x)[1]
    }) %>%
      tolower() %>%
      gsub(".formula", "", .)
  }

  var_i <- list()

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  wm <- which(clss == "gam")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      # prepare real data
      pred_test <- data.frame(
        pr_ab = data.frame(data)[, response],
        pred = suppressWarnings(
          mgcv::predict.gam(
            m[[i]],
            newdata = data,
            type = "response",
            se.fit = FALSE
          )
        )
      )
      var_r <- var_perf(pred_test)

      # Null model
      vi <- list()
      for (jj in 1:length(predictors)) {
        data_2 <- data

        var_rand <- foreach::foreach(
          j = 1:n_sim,
          .combine = rbind,
          .export = c("sdm_eval", "boyce"),
          .packages = c("dplyr")
        ) %dopar% {
          data_2[, predictors[jj]] <- data[sample(nrow(data)), predictors[jj]]

          x <- data.frame(
            pr_ab = data.frame(data_2)[, response],
            pred = suppressWarnings(
              mgcv::predict.gam(
                m[[i]],
                newdata = data_2,
                type = "response",
                se.fit = FALSE
              )
            )
          )

          var_perf(x)
        }

        # Calculate variable importance
        vi[[jj]] <- var_prep(var_rand, var_r)
      }

      names(vi) <- predictors
      var_i[[i]] <- dplyr::bind_rows(vi, .id = "predictors") %>%
        dplyr::mutate(model = i)
    }
  }

  #### gau models ####
  wm <- which(clss == "graf")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      # prepare real data
      pred_test <- data.frame(
        pr_ab = data.frame(data)[, response],
        pred = suppressWarnings(
          predict.graf(
            m[[i]],
            newdata = data.frame(data[, names(m[[i]]$obsx)]),
            type = "response",
            CI = NULL
          )[, 1]
        )
      )
      var_r <- var_perf(pred_test)

      # Null model
      vi <- list()
      for (jj in 1:length(predictors)) {
        data_2 <- data

        var_rand <- foreach::foreach(
          j = 1:n_sim,
          .combine = rbind,
          .export = c("sdm_eval", "boyce", "predict.graf"),
          .packages = c("dplyr")
        ) %dopar% {
          data_2[, predictors[jj]] <- data[sample(nrow(data)), predictors[jj]]

          x <- data.frame(
            pr_ab = data.frame(data_2)[, response],
            pred = suppressWarnings(
              predict.graf(
                m[[i]],
                newdata = data.frame(data_2[, names(m[[i]]$obsx)]),
                type = "response",
                CI = NULL
              )[, 1]
            )
          )

          var_perf(x)
        }

        # Calculate variable importance
        vi[[jj]] <- var_prep(var_rand, var_r)
      }
      names(vi) <- predictors
      var_i[[i]] <- dplyr::bind_rows(vi, .id = "predictors") %>%
        dplyr::mutate(model = i)
    }
  }

  #### glm models ####
  wm <- which(clss == "glm")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      # prepare real data
      pred_test <- data.frame(
        pr_ab = data.frame(data)[, response],
        pred = suppressWarnings(
          stats::predict.glm(
            m[[i]],
            newdata = data,
            type = "response",
            se.fit = FALSE
          )
        )
      )

      var_r <- var_perf(pred_test)

      # Null model
      vi <- list()
      for (jj in 1:length(predictors)) {
        data_2 <- data

        var_rand <- foreach::foreach(
          j = 1:n_sim,
          .combine = rbind,
          .export = c("sdm_eval", "boyce"),
          .packages = c("dplyr")
        ) %dopar% {
          data_2[, predictors[jj]] <- data[sample(nrow(data)), predictors[jj]]

          x <- data.frame(
            pr_ab = data.frame(data_2)[, response],
            pred = suppressWarnings(
              stats::predict.glm(
                m[[i]],
                newdata = data_2,
                type = "response",
                se.fit = FALSE
              )
            )
          )

          var_perf(x)
        }

        # Calculate variable importance
        vi[[jj]] <- var_prep(var_rand, var_r)
      }

      names(vi) <- predictors
      var_i[[i]] <- dplyr::bind_rows(vi, .id = "predictors") %>%
        dplyr::mutate(model = i)
    }
  }

  #### gbm models ####
  wm <- which(clss == "gbm")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      # prepare real data
      pred_test <- data.frame(
        pr_ab = data.frame(data)[, response],
        pred = suppressMessages(
          gbm::predict.gbm(
            m[[i]],
            newdata = data,
            type = "response",
            se.fit = FALSE
          )
        )
      )

      var_r <- var_perf(pred_test)

      # Null model
      vi <- list()
      for (jj in 1:length(predictors)) {
        data_2 <- data
        var_rand <- foreach::foreach(
          j = 1:n_sim,
          .combine = rbind,
          .export = c("sdm_eval", "boyce"),
          .packages = c("dplyr")
        ) %dopar% {
          data_2[, predictors[jj]] <- data[sample(nrow(data)), predictors[jj]]

          x <- data.frame(
            pr_ab = data.frame(data_2)[, response],
            pred = suppressMessages(
              gbm::predict.gbm(
                m[[i]],
                newdata = data_2,
                type = "response",
                se.fit = FALSE
              )
            )
          )

          var_perf(x)
        }

        # Calculate variable importance
        vi[[jj]] <- var_prep(var_rand, var_r)
      }

      names(vi) <- predictors
      var_i[[i]] <- dplyr::bind_rows(vi, .id = "predictors") %>%
        dplyr::mutate(model = i)
    }
  }

  #### raf models ####
  wm <- which(clss == "randomforest")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      # prepare real data
      pred_test <- data.frame(
        pr_ab = data.frame(data)[, response],
        pred = stats::predict(m[[i]], newdata = data, type = "prob")[, 2]
      )

      var_r <- var_perf(pred_test)

      # Null model
      vi <- list()
      for (jj in 1:length(predictors)) {
        data_2 <- data
        var_rand <- foreach::foreach(
          j = 1:n_sim,
          .combine = rbind,
          .export = c("sdm_eval", "boyce"),
          .packages = c("dplyr")
        ) %dopar% {
          data_2[, predictors[jj]] <- data[sample(nrow(data)), predictors[jj]]

          x <- data.frame(
            pr_ab = data.frame(data_2)[, response],
            pred = stats::predict(m[[i]], newdata = data_2, type = "prob")[, 2]
          )

          var_perf(x)
        }

        # Calculate variable importance
        vi[[jj]] <- var_prep(var_rand, var_r)
      }

      names(vi) <- predictors
      var_i[[i]] <- dplyr::bind_rows(vi, .id = "predictors") %>%
        dplyr::mutate(model = i)
    }
  }

  #### max models ####
  wm <- which(clss == "maxnet")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      # prepare real data
      pred_test <- data.frame(
        pr_ab = data.frame(data)[, response],
        pred = predict_maxnet(
          object = m[[i]],
          newdata = data.frame(data),
          clamp = clamp,
          type = pred_type
        )
      )

      var_r <- var_perf(pred_test)

      # Null model
      vi <- list()
      for (jj in 1:length(predictors)) {
        data_2 <- data
        var_rand <- foreach::foreach(
          j = 1:n_sim,
          .combine = rbind,
          .export = c("sdm_eval", "boyce", "predict_maxnet"),
          .packages = c("dplyr")
        ) %dopar% {
          data_2[, predictors[jj]] <- data[sample(nrow(data)), predictors[jj]]

          x <- data.frame(
            pr_ab = data.frame(data_2)[, response],
            predict_maxnet(
              object = m[[i]],
              newdata = data.frame(data_2),
              clamp = clamp,
              type = pred_type
            )
          )

          var_perf(x)
        }

        # Calculate variable importance
        vi[[jj]] <- var_prep(var_rand, var_r)
      }

      names(vi) <- predictors
      var_i[[i]] <- dplyr::bind_rows(vi, .id = "predictors") %>%
        dplyr::mutate(model = i)
    }
  }

  #### net models ####
  wm <- which(clss == "nnet")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      # prepare real data
      pred_test <- data.frame(
        pr_ab = data.frame(data)[, response],
        pred = suppressMessages(stats::predict(
          m[[i]],
          newdata = data, type = "raw"
        ))
      )

      var_r <- var_perf(pred_test)

      # Null model
      vi <- list()
      for (jj in 1:length(predictors)) {
        data_2 <- data
        var_rand <- foreach::foreach(
          j = 1:n_sim,
          .combine = rbind,
          .export = c("sdm_eval", "boyce"),
          .packages = c("dplyr")
        ) %dopar% {
          data_2[, predictors[jj]] <- data[sample(nrow(data)), predictors[jj]]

          x <- data.frame(
            pr_ab = data.frame(data_2)[, response],
            pred = suppressMessages(stats::predict(
              m[[i]],
              newdata = data_2, type = "raw"
            ))
          )

          var_perf(x)
        }

        # Calculate variable importance
        vi[[jj]] <- var_prep(var_rand, var_r)
      }

      names(vi) <- predictors
      var_i[[i]] <- dplyr::bind_rows(vi, .id = "predictors") %>%
        dplyr::mutate(model = i)
    }
  }

  #### svm models ####
  wm <- which(clss == "ksvm")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      # prepare real data
      pred_test <- data.frame(
        pr_ab = data.frame(data)[, response],
        pred = suppressMessages(kernlab::predict(
          m[[i]],
          newdata = data, type = "prob",
        )[, 2])
      )


      var_r <- var_perf(pred_test)

      # Null model
      vi <- list()
      for (jj in 1:length(predictors)) {
        data_2 <- data
        var_rand <- foreach::foreach(
          j = 1:n_sim,
          .combine = rbind,
          .export = c("sdm_eval", "boyce"),
          .packages = c("dplyr")
        ) %dopar% {
          data_2[, predictors[jj]] <- data[sample(nrow(data)), predictors[jj]]

          x <- data.frame(
            pr_ab = data.frame(data_2)[, response],
            pred = suppressMessages(
              kernlab::predict(m[[i]], newdata = data_2, type = "prob", )[, 2]
            )
          )


          var_perf(x)
        }

        vi[[jj]] <- var_prep(var_rand, var_r)
      }

      names(vi) <- predictors
      var_i[[i]] <- dplyr::bind_rows(vi, .id = "predictors") %>%
        dplyr::mutate(model = i)
    }
  }

  parallel::stopCluster(cl)

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
  names(var_i) <- dplyr::left_join(data.frame(alg = clss), df, by = "alg")[, 2]

  if (!is.null(esm)) {
    names(var_i) <- paste("esm", names(var_i), 1:length(var_i), sep = "_")
  }

  result <- dplyr::bind_rows(var_i, .id = "model") %>%
    relocate(model, threshold, predictors) %>%
    as_tibble()

  # In result, for metrics except IAME, if the value is less than 0, transform that value to 0
  nsm <- result %>%
    dplyr::select(TPR:IMAE) %>%
    names()
  result[, nsm][result[, nsm] < 0] <- 0

  return(result)
}
