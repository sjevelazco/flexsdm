#' Ensemble model fitting and validation
#'
#' @param models list. A list of models fitted with fit_ or tune_ function family. Models used for ensemble must have the same presences-absences records, partition methods, and threshold types.
#' @param ens_method character. Method used to create ensemble of different models. A vector must be provided for this argument. For meansup, meanw or pcasup method, it is necessary to provide an evaluation metric and threshold in 'metric' and 'thr_model' arguments respectively. By default all of the following ensemble methods will be performed:
#'   \itemize{
#'   \item mean: Simple average of the different models.
#'   \item meanw: Weighted average of models based on their performance. An evaluation metric and threshold type must be provided.
#'   \item meansup: Average of the best models (those with the evaluation metric above the average). An evaluation metric must be provided.
#'   \item meanthr: Averaging performed only with those cells with suitability values above the selected threshold.
#'   \item median: Median of the different models.}
#'   Usage ensemble = "meanthr". If several ensemble methods are to be implemented it is necessary to concatenate them, e.g., ensemble = c("meanw", "meanthr", "median")
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is useful for threshold-dependent performance metrics. It is possible to use more than one threshold criterion. A vector must be provided for this argument. The following threshold criteria are available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which Jaccard is the highest.
#'   \item max_sorensen: The threshold at which Sorensen is highest.
#'   \item max_fpb: The threshold at which FPB (F-measure on presence-background data) is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to sensitivity value. If a sensitivity values is not specified, default is 0.9.
#'   }
#' In the case of using more than one threshold type it is necessary concatenate threshold types, e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function will use all thresholds if no threshold is specified.
#'
#' @param thr_model character. This threshold is needed for conduct meanw, meandsup, and meanthr ensemble methods. It is mandatory to use only one threshold, and this must be the same threshold used to fit all the models used in the "models" argument. Usage thr_model = 'equal_sens_spec'
#' @param metric character. Performance metric used for selecting the best combination of hyper-parameter values. One of the following metrics can be used: SORENSEN, JACCARD, FPB, TSS, KAPPA, AUC, IMAE, and BOYCE. Default TSS. Usage metric = BOYCE
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item models: A list of models used for performing ensemble.
#' \item thr_metric: Threshold and metric specified in the function.
#' \item predictors: A tibble of quantitative (column names with c) and qualitative (column names with f) variables used in each models.
#' \item performance: A tibble with performance metrics (see \code{\link{sdm_eval}}).
#' \item performance_part: Performance metric for each replica and partition (see \code{\link{sdm_eval}}).
#' Those metrics that are threshold-dependent are calculated based on the threshold specified in the argument.
#' }
#'
#'
#' @importFrom dplyr %>% filter pull bind_rows mutate inner_join select starts_with bind_cols across group_by summarise left_join relocate
#' @importFrom stats median sd
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' require(terra)
#'
#' # Environmental variables
#' somevar <-
#'   system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' # Species occurrences
#' data("spp")
#' set.seed(1)
#' some_sp <- spp %>%
#'   dplyr::filter(species == "sp2") %>%
#'   sdm_extract(
#'     data = .,
#'     x = "x",
#'     y = "y",
#'     env_layer = somevar,
#'     variables = names(somevar),
#'     filter_na = TRUE
#'   ) %>%
#'   part_random(
#'     data = .,
#'     pr_ab = "pr_ab",
#'     method = c(method = "kfold", folds = 3)
#'   )
#'
#'
#' # gam
#' mglm <- fit_glm(
#'   data = some_sp,
#'   response = "pr_ab",
#'   predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
#'   partition = ".part",
#'   poly = 2
#' )
#' mraf <- fit_raf(
#'   data = some_sp,
#'   response = "pr_ab",
#'   predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
#'   partition = ".part",
#' )
#' mgbm <- fit_gbm(
#'   data = some_sp,
#'   response = "pr_ab",
#'   predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
#'   partition = ".part"
#' )
#'
#' # Fit and validate ensemble model
#' mensemble <- fit_ensemble(
#'   models = list(mglm, mraf, mgbm),
#'   ens_method = "meansup",
#'   thr = NULL,
#'   thr_model = "max_sens_spec",
#'   metric = "TSS"
#' )
#' 
#' mensemble
#' 
#' }
fit_ensemble <-
  function(models,
           ens_method = c("mean", "meanw", "meansup", "meanthr", "median"),
           thr = NULL,
           thr_model = NULL,
           metric = NULL) {
    . <- thr_value <- pr_ab <- rnames <- replicates <- TPR <- IMAE <- n_absences <- model <- thr_value <- n_presences <- n_absences <- NULL

    if (any(c("meanw", "meansup", "meanthr") %in% ens_method)) {
      if (is.null(thr_model) | is.null(metric)) {
        stop("for 'meanw', 'meansup', or 'meanthr' ensemble methods it is necessary to provide a threshold type in 'thr_model' and 'metric' argument")
      }
    }

    #### Models names
    nms <- paste0("m_", 1:length(models))


    if (any(c("meanw", "meansup", "meanthr") %in% ens_method)) {
      #### Performance metric
      metric <- paste0(metric, "_mean")

      #### Model performances
      perf <- sapply(models, function(x) {
        x[["performance"]] %>%
          dplyr::filter(threshold == {{ thr_model }}) %>%
          dplyr::pull({{ metric }})
      })

      #### Model thresholds
      thr_v <- sapply(models, function(x) {
        x[["performance"]] %>%
          dplyr::filter(threshold == {{ thr_model }}) %>%
          dplyr::pull(thr_value)
      })

      names(perf) <- names(thr_v) <- nms
    }

    # Variables used in each models
    variables <- lapply(models, function(x) {
      x$predictors
    }) %>%
      dplyr::bind_rows()

    #### Extract and merge suitability databases
    data_ens <- sapply(models, function(x) {
      x["data_ens"]
    })

    data_ens <- mapply(function(x, cn) {
      colnames(x)[colnames(x) %in% "pred"] <- cn
      x
    }, data_ens, nms, SIMPLIFY = FALSE)

    data_ens <- lapply(data_ens, function(x) {
      x %>% dplyr::mutate(pr_ab = pr_ab %>%
        as.character() %>%
        as.double())
    })

    data_ens2 <-
      dplyr::full_join(data_ens[[1]],
        data_ens[[2]],
        by = c("rnames", "replicates", "part", "pr_ab")
      )
    if (length(data_ens) > 2) {
      for (i in 3:length(data_ens)) {
        data_ens2 <-
          dplyr::full_join(data_ens2,
            data_ens[[i]],
            by = c("rnames", "replicates", "part", "pr_ab")
          )
      }
    }
    rm(data_ens)

    #### Extract predicted suitability of each model
    values <- data_ens2 %>%
      dplyr::select(dplyr::starts_with("m_"))

    #### Remove suitability values from data_ens2
    data_ens2 <- data_ens2 %>% dplyr::select(-dplyr::starts_with("m_"))

    #### Perform ensembles
    v <- as.list(rep("x", 5))
    names(v) <- c("mean", "meanw", "meansup", "meanthr", "median")

    if (any("mean" == ens_method)) {
      v[["mean"]] <- apply(values, 1, function(x) {
        mean(x, na.rm = TRUE)
      })
    }
    if (any("meanw" == ens_method)) {
      v[["meanw"]] <- mapply(function(x, v) {
        (x * v)
      }, values, perf, SIMPLIFY = TRUE) %>%
        apply(., 1, function(x) {
          sum(x, na.rm = TRUE)
        }) / sum(perf)
    }
    if (any("meansup" == ens_method)) {
      v[["meansup"]] <- apply(values[, perf >= mean(perf)], 1, function(x) {
        mean(x, na.rm = TRUE)
      })
    }
    if (any("meanthr" == ens_method)) {
      v[["meanthr"]] <- mapply(function(x, v) {
        ifelse(x >= v, x, 0)
      }, values, thr_v, SIMPLIFY = TRUE) %>%
        apply(., 1, function(x) {
          mean(x, na.rm = TRUE)
        })
    }
    if (any("median" == ens_method)) {
      v[["median"]] <- apply(values, 1, function(x) {
        stats::median(x, na.rm = TRUE)
      })
    }

    v <- v[!sapply(v, is.character)]
    v <- dplyr::bind_rows(v)
    data_ens2 <- dplyr::bind_cols(data_ens2, v)
    rm(list = c("values", "v"))

    #### Calculate ensemble performance
    p_names <- data_ens2 %>%
      dplyr::pull(replicates) %>%
      unique()
    np <- length(p_names)

    ##### average ensemble prediction for calculating model threshold
    data_ens3 <- data_ens2 %>%
      group_by(rnames, pr_ab) %>%
      summarise(dplyr::across(
        {{ ens_method }},
        list(mean = function(x) {
          mean(x)
        })
      ), .groups = "drop") %>%
      select(-rnames)
    colnames(data_ens3) <- gsub("_mean", "", colnames(data_ens3))

    #### Objects to store outputs
    eval_partial_list <- list()
    eval_partial_list_2 <- as.list(ens_method)
    names(eval_partial_list_2) <- ens_method
    threshold <- ensemble <- as.list(rep(NA, length(ens_method)))
    names(threshold) <- names(ensemble) <- ens_method

    pb <- utils::txtProgressBar(min = 0, max = length(ens_method), style = 3)
    for (g in ens_method) {
      for (h in 1:np) {
        # message("\n", "Replica number: ", h, "/", np)

        pred_test <- data_ens2 %>%
          dplyr::filter(replicates == p_names[h])
        pred_test <- split(pred_test, f = pred_test$part)
        np2 <- length(pred_test)
        eval_partial <- as.list(rep(NA, np2))

        for (i in 1:np2) {
          # message("Partition number: ", i, "/", np2)
          # Validation of model
          tryCatch(
            {
              suppressWarnings(eval <-
                sdm_eval(
                  p = pred_test[[i]] %>% dplyr::filter(pr_ab == 1) %>% dplyr::pull({{ g }}),
                  a = pred_test[[i]] %>% dplyr::filter(pr_ab == 0) %>% dplyr::pull({{ g }}),
                  thr = thr
                ))
              eval_partial[[i]] <- eval
              rm(eval)
            },
            error = function(cond) {}
          )
        }
        names(eval_partial) <- 1:np2

        # Create final database with parameter performance
        eval_partial <- eval_partial[sapply(eval_partial, function(x) !is.null(dim(x)))] %>%
          dplyr::bind_rows(., .id = "partition")
        eval_partial_list[[h]] <- eval_partial
      }

      eval_partial <- eval_partial_list %>%
        dplyr::bind_rows(., .id = "replica")

      eval_partial_list_2[[g]] <- eval_partial

      eval_final <- eval_partial %>%
        dplyr::group_by(threshold) %>%
        dplyr::summarise(dplyr::across(
          TPR:IMAE,
          list(mean = mean, sd = stats::sd)
        ), .groups = "drop")

      ensemble[[g]] <- eval_final

      threshold[[g]] <- sdm_eval(
        p = data_ens3 %>% dplyr::filter(pr_ab == 1) %>% dplyr::pull({{ g }}),
        a = data_ens3 %>% dplyr::filter(pr_ab == 0) %>% dplyr::pull({{ g }}),
        thr = thr
      )
      utils::setTxtProgressBar(pb, which(ens_method == g))
    }
    close(pb)

    # Threshold
    threshold <- lapply(threshold, function(x) x %>% dplyr::select(threshold:n_absences)) %>%
      dplyr::bind_rows(., .id = "model")

    # Bind ensemble performance
    ensemble <- dplyr::bind_rows(ensemble, .id = "model")
    eval_partial_list_2 <- dplyr::bind_rows(eval_partial_list_2, .id = "model")
    
    #### Model object
    m <- lapply(models, function(x) x[c("model", "performance")])
    names(m) <- nms

    result <- list(
      models = m,
      thr_metric = c(thr_model, metric),
      predictors = variables,
      performance = dplyr::left_join(ensemble, threshold, by = c("model", "threshold")) %>%
        dplyr::relocate(model, threshold, thr_value, n_presences, n_absences),
      performance_part = eval_partial_list_2
    )

    return(result)
  }
