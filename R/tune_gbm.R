#' Function for constructing Generalized Boosted Regression models with exploration of hyper-parameters
#'
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous or discrete variables). Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param partition character. Column name with training and validation partition groups.
#' @param grid data.frame. Provide a data frame object with algorithm hyper-parameters values to be tested. It Is recommended to generate this data.frame with grid() function. Hyper-parameters needed for tuning are 'n.trees', 'shrinkage', and 'n.minobsinnode'.
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The next threshold area available:
#' \itemize{
#'   \item LPT: The highest threshold at which there is no omission. Usage thr=c(type='LPT').
#'   \item EQUAL_SENS_SPEC: Threshold at which the sensitivity and specificity are equal.
#'   \item MAX_TSS: Threshold at which the sum of the sensitivity and specificity is the highest.
#'   Usage thr=c(type='MAX_TSS').
#'   \item MAX_KAPPA: The threshold at which kappa is the highest ("max kappa"). Usage thr=c(type='MAX_KAPPA').
#'   \item MAX_JACCARD: The threshold at which Jaccard is the highest. Usage thr=c(type='MAX_JACCARD').
#'   \item MAX_SORENSEN: The threshold at which Sorensen is highest. Usage thr=c(type='MAX_SORENSEN').
#'   \item MAX_FPB: The threshold at which FPB is highest. Usage thr=c(type='MAX_FPB').
#'   \item SENSITIVITY: A threshold value specified by user. Usage thr=c(type='SENSITIVITY', sens='0.6'). 'sens' refers to models will be binarized using this suitability value.
#'   }
#' @param metric character. Performance metric used for selecting the best combination of hyper-parameter values. Can be used one of the next metrics SORENSEN, JACCARD, FPB, TSS, KAPPA, AUC, and BOYCE. TSS is used as default.
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "gbm" class object. This object can be used for predicting.
#' \item tune_performance: Performance metric (see \code{\link{enm_eval}}) for each combination of the hyper-parameters.
#' \item best_hyper: Hyper-parameters values and performance metric (see \code{\link{enm_eval}}) for the best hyper-parameters combination.
#' \item selected_threshold: Value of the threshold selected.
#' \item threshold_table: Value of all threshold.
#' }
#'
#' @export
#'
#' @importFrom dplyr select starts_with bind_rows tibble group_by_at summarise across everything pull
#' @importFrom gbm gbm predict.gbm
#' @importFrom stats formula na.omit
#'
#' @family aggregate functions
#' @seealso \code{\link{tune_mx}}, \code{\link{tune_nnet}}, \code{\link{tune_rf}}, and \code{\link{tune_svm}}.
#'
#' @examples
#' \dontrun{
#' data(abies_db)
#' abies_db
#'
#' # We will partition the data with the k-fold method
#'
#' abies_db2 <- data_part(
#' data = abies_db,
#' p_a = 'pr_ab',
#' bg_data = NULL,
#' bg_a = NULL,
#' method = c(method = "KFOLD", folds = 10)
#' )
#'
#' # pr_ab columns is species presence and absences (i.e. the response variable)
#' # partition is columns with partition groups for performing 5-fold cross validation
#' # from aet to landform are the predictors variables (landform is a qualitative variable)
#'
#' # Hyper-parameter values for tuning
#' tune_grid <-
#'   expand.grid(
#'     n.trees = c(20, 50, 100, 200),
#'     shrinkage = c(0.1, 0.5, 1),
#'     n.minobsinnode = c(1, 3, 5, 7, 9)
#'   )
#'
#' gbm_t <-
#'   tune_gbm(
#'     data = abies_db2,
#'     response = "pr_ab",
#'     predictors = c(
#'       "aet", "cwd", "tmin", "ppt_djf",
#'       "ppt_jja", "pH", "awc", "depth", "percent_clay"
#'     ),
#'     predictors_f = c("landform"),
#'     partition = "partition",
#'     grid = tune_grid,
#'     thr = "MAX_TSS",
#'     metric = "TSS",
#'   )
#'
#' # Outputs
#' gbm_t$model
#' gbm_t$tune_performance
#' gbm_t$best_hyper
#' gbm_t$selected_threshold
#' gbm_t$threshold_table
#'
#' # Graphical exploration of performance of each hyper-parameter setting
#' require(ggplot2)
#' ggplot(gbm_t$tune_performance, aes(n.minobsinnode, TSS_mean)) +
#'   geom_point(aes(col = factor(n.trees)))
#'
#' pg <- position_dodge(width = 0.5)
#' ggplot(gbm_t$tune_performance, aes(factor(n.minobsinnode),
#'   TSS_mean,
#'   col = factor(shrinkage)
#' )) +
#'   geom_errorbar(aes(ymin = TSS_mean - TSS_sd, ymax = TSS_mean + TSS_sd), width = 0.2, position = pg) +
#'   geom_point(position = pg) +
#'   geom_line(
#'     data = gbm_t$tune_performance,
#'     aes(as.numeric(factor(n.minobsinnode)),
#'       TSS_mean,
#'       col = factor(shrinkage)
#'     ), position = pg
#'   ) +
#'   facet_wrap(. ~ n.trees) +
#'   theme(legend.position = "bottom")
#' }
#'
tune_gbm <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           partition,
           grid = NULL,
           thr = NULL,
           metric = "TSS",
           ...) {
    data <- data.frame(data)

    if (is.null(predictors_f)) {
      data <- data %>%
        dplyr::select(response, predictors, dplyr::starts_with(partition))
      data <- data.frame(data)
    } else {
      data <- data %>%
        dplyr::select(response, predictors, predictors_f, dplyr::starts_with(partition))
      data <- data.frame(data)
      for (i in predictors_f) {
        data[, i] <- as.factor(data[, i])
      }
    }

    # Formula
    Fmula <- stats::formula(paste(
      response, "~",
      paste(c(predictors, predictors_f), collapse = " + ")
    ))

    # Prepare grid when grid=default or NULL
    if (is.null(grid)) {
      nv <- length(stats::na.omit(c(predictors, predictors_f)))
      grid <- data.frame(n.trees = 100, shrinkage = 0.1, n.minobsinnode = 10)
    }
    if (class(grid) == "character") {
      nv <- length(stats::na.omit(c(predictors, predictors_f)))
      if (grid == "defalut") {
        grid <- expand.grid(
          n.trees = c(20, 50, 100, 200, 500),
          shrinkage = c(0.01, 0.1, 0.3),
          n.minobsinnode = c(1:10, 20, 30) # try define a best default tuning set
        )
      }
    }

    # Test hyper-parameters names
    hyperp <- names(grid)
    if (!all(c("n.trees", "shrinkage", "n.minobsinnode") %in% hyperp)) {
      stop("Database used in 'grid' argument has to contain this columns for tunning: 'n.trees', 'shrinkage', 'n.minobsinnode'")
    }

    grid$tune <- 1:nrow(grid)

    # Fit models
    np <- ncol(data %>% dplyr::select(dplyr::starts_with(partition)))
    p_names <- names(data %>% dplyr::select(dplyr::starts_with(partition)))
    eval_partial_list <- list()
    for (h in 1:np) {
      message("Replica number: ", h, "/", np)

      out <- pre_tr_te(data, p_names, h)
      train <- out$train
      test <- out$test
      np2 <- out$np2
      rm(out)

      eval_partial <- list()

      for (i in 1:np2) {
        message("Partition number: ", i, "/", np2)
        mod <- as.list(rep(NA, nrow(grid)))
        names(mod) <- 1:nrow(grid)
        for (ii in 1:nrow(grid)) {
          set.seed(1)
          try(mod[[ii]] <-
            suppressMessages(
              gbm::gbm(
                Fmula,
                data = train[[i]],
                distribution = "bernoulli",
                bag.fraction = 1, # Explore more this parameter
                n.trees = grid$n.trees[ii],
                interaction.depth = 1,
                shrinkage = grid$shrinkage[ii],
                n.minobsinnode = grid$n.minobsinnode[ii]
              )
            ))
        }

        filt <- sapply(mod, function(x) class(x) == "gbm")
        mod <- mod[filt]
        grid2 <- grid[filt, ]

        # Predict for presences absences data
        pred_test <-
          lapply(mod, function(x) {
            data.frame(
              pr_ab = test[[i]][, response],
              pred = suppressMessages(
                gbm::predict.gbm(
                  x,
                  newdata = test[[i]],
                  type = "response"
                )
              )
            )
          })

        # Validation of parameter combination
        eval <- list()
        for (ii in 1:length(pred_test)) {
          eval[[ii]] <-
            enm_eval(
              p = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 1],
              a = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 0],
              thr = thr
            )
        }

        eval <- dplyr::bind_rows(lapply(eval, function(x) x$selected_threshold))
        eval <- dplyr::tibble(cbind(grid2, eval))
        eval_partial[[i]] <- eval
      }

      # Create final database with parameter performance
      names(eval_partial) <- 1:np2
      eval_partial <- eval_partial %>%
        dplyr::bind_rows(., .id = "partition")
      eval_partial_list[[h]] <- eval_partial
    }

    eval_partial <- eval_partial_list %>%
      dplyr::bind_rows(., .id = "replica")

    eval_final <- eval_partial %>%
      dplyr::select(-replica, -partition, -c(tune:n_absences)) %>%
      dplyr::group_by_at(hyperp) %>%
      dplyr::summarise(dplyr::across(
        dplyr::everything(),
        list(mean = mean, sd = sd)
      ), .groups = "drop")

    # Find the bets parameter setting
    filt <- eval_final %>% dplyr::pull(paste0(metric, "_mean"))
    filt <- which.max(filt)
    best_tune <- eval_final[filt, ]
    best_hyperp <- eval_final[filt, hyperp]


    # Fit final models with best settings
    set.seed(1)
    mod <-
      suppressMessages(
        gbm::gbm(
          Fmula,
          data = data,
          distribution = "bernoulli",
          n.trees = best_hyperp$n.trees,
          interaction.depth = 1,
          shrinkage = best_hyperp$shrinkage,
          n.minobsinnode = best_hyperp$n.minobsinnode
        )
      )


    pred_test <- data.frame(
      pr_ab = data[, response],
      pred = suppressMessages(gbm::predict.gbm(
        mod,
        newdata = data,
        type = "response"
      ))
    )

    threshold <- enm_eval(
      p = pred_test$pred[pred_test$pr_ab == 1],
      a = pred_test$pred[pred_test$pr_ab == 0],
      thr = thr
    )

    result <- list(
      model = mod,
      tune_performance = eval_final,
      best_hyper_performance = best_tune,
      selected_threshold = threshold[[1]] %>% dplyr::select(threshold:TNR ),
      threshold_table = threshold[[2]]
    )
    return(result)
  }
