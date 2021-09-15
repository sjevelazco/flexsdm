#' Fit and validate Generalized Boosted Regression models with exploration of
#' hyper-parameters that optimize performance
#'
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative predictor variables
#' (i.e. continuous variables). Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables
#' (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor
#' variables (e.g. formula(pr_ab ~ aet + ppt_jja + pH + awc + depth + landform)).
#' Note that the variables used here must be consistent with those used in
#' response, predictors, and predictors_f arguments. Default is NULL.
#' @param partition character. Column name with training and validation partition groups.
#' @param grid data.frame. A data frame object with algorithm hyper-parameter values to be tested.
#' It Is recommended to generate this data.frame with the grid() function. Hyper-parameters needed
#' for tuning are 'n.trees', 'shrinkage', and 'n.minobsinnode'.
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1) needed for
#' threshold-dependent performance metrics. It is possible to use more than one threshold type.
#' Provide a vector for this argument. The following threshold types are available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the
#'   highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which the Jaccard index is the highest.
#'   \item max_sorensen: The threshold at which the Sorensen index is highest.
#'   \item max_fpb: The threshold at which FPB (F-measure on presence-background data) is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to
#'   sensitivity value. If no sensitivity value is specified, the default used is 0.9
#'   }
#' If more than one threshold type is used they must be concatenate, e.g.,
#' thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity',
#' sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function will use all threshold
#' types if no threshold is specified.
#' @param metric character. Performance metric used for selecting the best combination of
#' hyper-parameter values. The following metrics can be used: SORENSEN, JACCARD, FPB, TSS,
#' KAPPA, AUC, and BOYCE. TSS is used as the default.
#' @param n_cores numeric. Number of cores use for parallelization. Default 1
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "gbm" class object. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names)
#' variables use for modeling.
#' \item performance: Hyper-parameter values and performance metric (see \code{\link{sdm_eval}})
#' for the best hyper-parameter combination.
#' \item hyper_performance: Performance metric (see \code{\link{sdm_eval}}) for each combination of
#' the hyper-parameters.
#' \item data_ens: Predicted suitability for each test partition based on the best model. This
#' database is used in \code{\link{fit_ensemble}}
#' }
#'
#' @export
#'
#' @importFrom dplyr %>% select starts_with bind_rows tibble group_by_at summarise across everything pull
#' @importFrom gbm gbm predict.gbm
#' @importFrom stats formula na.omit
#'
#' @seealso \code{\link{tune_max}}, \code{\link{tune_net}}, \code{\link{tune_raf}},
#' and \code{\link{tune_svm}}.
#'
#' @examples
#' \dontrun{
#' data(abies)
#' abies
#'
#' # Partition the data with the k-fold method
#'
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 5)
#' )
#'
#' # pr_ab is the name of the column with species presence and absences (i.e. the response variable)
#' # from aet to landform are the predictors variables (landform is a qualitative variable)
#'
#' # Hyper-parameter values for tuning
#' tune_grid <-
#'   expand.grid(
#'     n.trees = c(20, 50, 100),
#'     shrinkage = c(0.1, 0.5, 1),
#'     n.minobsinnode = c(1, 3, 5, 7, 9)
#'   )
#'
#' gbm_t <-
#'   tune_gbm(
#'     data = abies2,
#'     response = "pr_ab",
#'     predictors = c(
#'       "aet", "cwd", "tmin", "ppt_djf", "ppt_jja",
#'       "ppt_jja", "pH", "awc", "depth"
#'     ),
#'     predictors_f = c("landform"),
#'     partition = ".part",
#'     grid = tune_grid,
#'     thr = "max_sens_spec",
#'     metric = "TSS",
#'    n_cores = 1
#'   )
#'
#' # Outputs
#' gbm_t$model
#' gbm_t$predictors
#' gbm_t$performance
#' gbm_t$data_ens
#' gbm_t$hyper_performance
#'
#' # Graphical exploration of performance of each hyper-parameter setting
#' require(ggplot2)
#' pg <- position_dodge(width = 0.5)
#' ggplot(gbm_t$hyper_performance, aes(factor(n.minobsinnode),
#'   TSS_mean,
#'   col = factor(shrinkage)
#' )) +
#'   geom_errorbar(aes(ymin = TSS_mean - TSS_sd, ymax = TSS_mean + TSS_sd),
#'     width = 0.2, position = pg
#'   ) +
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
           fit_formula = NULL,
           partition,
           grid = NULL,
           thr = NULL,
           metric = "TSS",
           n_cores = 1) {
    . <- model <- TPR <- IMAE <- thr_value <- n_presences <- n_absences <- NULL
    variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))

    data <- data.frame(data)

    # Test response variable
    r_test <- (data %>% dplyr::pull(response) %>% unique() %>% na.omit())
    if ((!all(r_test %in% c(0, 1)))) {
      stop("values of response variable do not match with 0 and 1")
    }

    if (is.null(predictors_f)) {
      data <- data %>%
        dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), dplyr::starts_with(partition))
      data <- data.frame(data)
    } else {
      data <- data %>%
        dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), dplyr::all_of(predictors_f), dplyr::starts_with(partition))
      data <- data.frame(data)
      for (i in predictors_f) {
        data[, i] <- as.factor(data[, i])
      }
    }

    # Remove NAs
    complete_vec <- stats::complete.cases(data[, c(response, unlist(variables))])
    if (sum(!complete_vec) > 0) {
      message(sum(!complete_vec), " rows were excluded from database because NAs were found")
      data <- data %>% dplyr::filter(complete_vec)
    }
    rm(complete_vec)

    # Formula
    if (is.null(fit_formula)) {
      formula1 <- stats::formula(paste(
        response, "~",
        paste(c(predictors, predictors_f), collapse = " + ")
      ))
    } else {
      formula1 <- fit_formula
    }
    message(
      "Formula used for model fitting:\n",
      Reduce(paste, deparse(formula1)) %>% gsub(paste("  ", "   ", collapse = "|"), " ", .),
      "\n"
    )

    # Prepare grid when grid=default or NULL
    if (is.null(grid)) {
      nv <- length(stats::na.omit(c(predictors, predictors_f)))
      grid <- expand.grid(
        n.trees = c(20, 50, 100, 200, 500),
        shrinkage = c(0.01, 0.1, 0.3),
        n.minobsinnode = c(1:10, 20, 30) # TODO try define a best default tuning set
      )
      message("Hyper-parameter values were not provided, default values will be used")
      message("n.trees = c(20, 50, 100, 200, 500")
      message("shrinkage = c(0.01, 0.1, 0.3)")
      message("n.minobsinnode = c(1:10, 20, 30)")
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

    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    for (h in 1:np) {
      message("Replica number: ", h, "/", np)

      out <- pre_tr_te(data, p_names, h)
      train <- out$train
      test <- out$test
      np2 <- out$np2
      rm(out)

      # eval_partial <- as.list(rep(NA, np2))

      eval_partial <- foreach::foreach(i = 1:np2, .export=c('sdm_eval', 'boyce'), .packages = c("dplyr")) %dopar%{
        # message("Partition number: ", i, "/", np2)
        mod <- as.list(rep(NA, nrow(grid)))
        names(mod) <- 1:nrow(grid)
        for (ii in 1:nrow(grid)) {
          set.seed(1)
          try(mod[[ii]] <-
            suppressMessages(
              gbm::gbm(
                formula1,
                data = train[[i]],
                distribution = "bernoulli",
                bag.fraction = 0.9,
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
        tnames <- apply(grid2, 1, function(x) paste(x, collapse = "_"))

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
            sdm_eval(
              p = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 1],
              a = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 0],
              thr = thr
            ) %>% dplyr::tibble(model = "gbm", .)
        }

        names(eval) <- tnames
        eval <- dplyr::bind_rows(eval, .id = "tnames")

        eval <-
          dplyr::tibble(dplyr::left_join(dplyr::mutate(grid2, tnames),
            eval,
            by = "tnames"
          )) %>%
          dplyr::select(-tnames)
        eval
      }
      parallel::stopCluster(cl)

      # Create final database with parameter performance
      names(eval_partial) <- 1:np2
      eval_partial <- eval_partial[sapply(eval_partial, function(x) !is.null(dim(x)))] %>%
        dplyr::bind_rows(., .id = "partition")
      eval_partial_list[[h]] <- eval_partial
    }

    eval_partial <- eval_partial_list %>%
      dplyr::bind_rows(., .id = "replica")

    eval_final <- eval_partial %>%
      dplyr::group_by_at(c(hyperp, "model", "threshold")) %>%
      dplyr::summarise(dplyr::across(
        TPR:IMAE,
        list(mean = mean, sd = sd)
      ), .groups = "drop")

    # Find the bets parameter setting
    filt <- eval_final %>% dplyr::pull(paste0(metric, "_mean"))
    filt <- which.max(filt)
    best_tune <- eval_final[filt, ]
    best_hyperp <- eval_final[filt, hyperp]


    # Get data for ensemble
    pred_test_ens <- fit_gbm(
      data = data,
      response = response,
      predictors = predictors,
      predictors_f = predictors_f,
      partition = partition,
      thr = thr,
      fit_formula = formula1,
      n_trees = best_tune$n.trees,
      n_minobsinnode = best_tune$n.minobsinnode,
      shrinkage = best_tune$shrinkage
    )[["data_ens"]]


    # Fit final models with best settings
    set.seed(1)
    mod <-
      suppressMessages(
        gbm::gbm(
          formula1,
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

    threshold <- sdm_eval(
      p = pred_test$pred[pred_test$pr_ab == 1],
      a = pred_test$pred[pred_test$pr_ab == 0],
      thr = thr
    )

    result <- list(
      model = mod,
      predictors = variables,
      performance = dplyr::left_join(best_tune, threshold[1:4], by = "threshold") %>%
        dplyr::relocate(dplyr::all_of(hyperp), model, threshold, thr_value, n_presences, n_absences),
      hyper_performance = eval_final,
      data_ens = pred_test_ens
    )
    return(result)
  }
