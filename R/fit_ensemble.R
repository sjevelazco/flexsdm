#' Model assemble and validation
#'
#' @param models list. A list of models fitted with fit_ or tune_ function family. Models used for ensemble must have the same presences-absences records, partition methods, threshold types.
#' @param ensemble character. Method used to ensemble different models. A vector must be provided for this argument. For meansup, meanw or pcasup method, it is necessary to provide an evaluation metric to ensemble arguments (metric argument below). By default will be performed all ensemble methods:
#'   \itemize{
#'   \item mean: Simple average of the different models. Usage ensemble=c(method='mean').
#'   \item meanw: Weighted average of models based on their performance. An evaluation metric and threshold type must be provided..
#'   \item meansup: Average of the best models (e.g., TSS over the average). An evaluation metric must be provided.
#'   \item meanthr: Average performed only with those cells with suitability values above the selected threshold.
#'   \item median: Median of the different models.
#'   }
#'   Usage ensemble = "meanthr". In case is needed perform different ensemble method it is necessary concatenate them, e.g., ensemble = c("meanw", "meanthr", "median")
#' @param thr character. Threshold type used to get binary suitability values (i.e. 0,1). In fit_ensemble, only it possible to use one threshold type. It is mandatory that the threshold used was the same for all models (see thr argument in fit_ and tune_ family function). Threshold available: lpt, equal_sens_spec, max_sens_spec, max_kappa, max_jaccard, max_sorensen, max_fpb, and sensitivity. Usage thr = 'equal_sens_spec'
#' @param metric character. Performance metric used for selecting the best combination of hyper-parameter values. One of the next metrics can be used SORENSEN, JACCARD, FPB, TSS, KAPPA, AUC, IMAE, and BOYCE. Default TSS. Usage metric = BOYCE
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A list of models used for performing ensemble.
#' \item predictors: A tibble of quantitative (column names with c) and qualitative (column names with f) variables used in each models.
#' \item performance: A tibble with performance metric (see \code{\link{sdm_eval}}).
#' Those threshold dependent metrics are calculated based on the threshold specified in thr argument .
#' \item selected_thresholds: A tibble with value of the threshold selected.
#' \item threshold_table: A tibble with all threshold values.
#' }
#'
#'
#' @importFrom dplyr filter all_of pull bind_rows mutate inner_join select starts_with group_by summarise across everything
#' @importFrom stats sd
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @examples
fit_ensemble <-
  function(models,
           ens_method = c("mean", "meanw", "meansup", "meanthr", "median"),
           thr,
           metric = "TSS") {

    #### Models names
    nms <- paste0("m_", 1:length(models))

    #### Performance metric
    metric <- paste0(metric, "_mean")

    #### Model performances
    perf <- sapply(models, function(x) {
      x[["performance"]] %>%
        dplyr::filter(threshold == dplyr::all_of(thr)) %>%
        dplyr::pull(dplyr::all_of(metric))
    })

    #### Model thresholds
    thr_v <- sapply(models, function(x) {
      x[["selected_thresholds"]] %>%
        dplyr::filter(threshold == dplyr::all_of(thr)) %>%
        dplyr::pull(values)
    })

    names(perf) <- names(thr_v) <- nms

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
      dplyr::inner_join(data_ens[[1]],
        data_ens[[2]],
        by = c("rnames", "replicates", "part", "pr_ab")
      )
    if (length(data_ens) > 2) {
      for (i in 3:length(data_ens)) {
        data_ens2 <-
          dplyr::inner_join(data_ens2,
            data_ens[[2]],
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
    v <- as.list(rep('x', 5))
    names(v) <- c('mean', 'meanw', 'meansup', 'meanthr', 'median')

    if(any('mean' == ens_method)) {
      v[['mean']] <- apply(values, 1, function(x) {
        mean(x, na.rm = TRUE)
      })
    }
    if(any('meanw' == ens_method)) {
      v[['meanw']] <- mapply(function(x, v) {
        (x * v)
      }, values, perf, SIMPLIFY = TRUE) %>%
        apply(., 1, function(x) {
          mean(x, na.rm = TRUE)
        })
    }
    if(any('meansup' == ens_method)) {
      v[['meansup']] <- apply(values[, perf >= mean(perf)], 1, function(x) {
        mean(x, na.rm = TRUE)
      })
    }
    if(any('meanthr' == ens_method)) {
      v[['meanthr']] <- mapply(function(x, v) {
        ifelse(x >= v, x, 0)
      }, values, thr_v, SIMPLIFY = TRUE) %>%
        apply(., 1, function(x) {
          mean(x, na.rm = TRUE)
        })
    }
    if(any('median' == ens_method)) {
      v[['median']] <- apply(values, 1, function(x) {
        median(x, na.rm = TRUE)
      })
    }

    v <- v[!sapply(v, is.character)]
    v <- dplyr::bind_rows(v)
    data_ens2 <- dplyr::bind_cols(data_ens2, v)
    rm(list=c('values', 'v'))

    #### Calculate ensemble performance
    p_names <- data_ens2 %>%
      dplyr::pull(replicates) %>%
      unique()
    np <- length(p_names)

    eval_partial_list <- list()
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
        eval_partial <- list()

        for (i in 1:np2) {
          # message("Partition number: ", i, "/", np2)
          # Validation of model

          eval <-
            sdm_eval(
              p = pred_test[[i]] %>% dplyr::filter(pr_ab == 1) %>% dplyr::pull(dplyr::all_of(g)),
              a = pred_test[[i]] %>% dplyr::filter(pr_ab == 0) %>% dplyr::pull(dplyr::all_of(g)),
              thr = thr
            )
          eval_partial[[i]] <- eval$selected_thresholds
          names(eval_partial)[i] <- i
        }

        # Create final database with parameter performance
        eval_partial <- eval_partial %>%
          dplyr::bind_rows(., .id = "partition")
        eval_partial_list[[h]] <- eval_partial
      }

      eval_partial <- eval_partial_list %>%
        dplyr::bind_rows(., .id = "replica")

      eval_final <- eval_partial %>%
        dplyr::group_by(threshold) %>%
        dplyr::select(-c(replica:partition, values:n_absences)) %>%
        dplyr::summarise(dplyr::across(
          dplyr::everything(),
          list(mean = mean, sd = stats::sd)
        ), .groups = "drop")

      ensemble[[g]] <- eval_final

      threshold[[g]] <- sdm_eval(
        p = data_ens2 %>% dplyr::filter(pr_ab == 1) %>% dplyr::pull(dplyr::all_of(g)),
        a = data_ens2 %>% dplyr::filter(pr_ab == 0) %>% dplyr::pull(dplyr::all_of(g)),
        thr = thr
      )
      utils::setTxtProgressBar(pb, which(ens_method == g))
    }
    close(pb)

    # Threshold
    st <- lapply(threshold, function(x) x$selected_thresholds) %>%
      dplyr::bind_rows(., .id = "model")
    threshold <- lapply(threshold, function(x) x$all_thresholds) %>%
      dplyr::bind_rows(., .id = "model")

    ensemble <- dplyr::bind_rows(ensemble, .id = "model")

    #### Model object
    m <- lapply(models, function(x) x[[1]])
    names(m) <- nms

    result <- list(
      model = m,
      predictors = variables,
      performance = ensemble,
      selected_thresholds = st %>% dplyr::select(model, threshold:values),
      all_thresholds = threshold %>% dplyr::select(model, threshold:values)
    )
  }
