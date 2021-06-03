#' Fit and validate Random Forest models with exploration of hyper-parameters
#'
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous or discrete variables). Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor
#' variables (e.g. formula(pr_ab ~ aet + ppt_jja + pH + awc + depth + landform)).
#' Note that the variables used here must be consistent with those used in
#' response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param grid data.frame. Provide a data frame object with algorithm hyper-parameters values to be tested. It Is recommended to generate this data.frame. Hyper-parameter needed for tuning is 'mtry'. The maximum mtry must be equal to total number of predictors
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The next threshold area available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which Jaccard is the highest.
#'   \item max_sorensen: The threshold at which Sorensen is highest.
#'   \item max_fpb: The threshold at which FPB is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to sensitivity value. If it is not specified a sensitivity values, function will use by default 0.9
#'   }
#' In the case of use more than one threshold type it is necessary concatenate threshold types, e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function will use all thresholds if no threshold is specified
#' @param metric character. Performance metric used for selecting the best combination of hyper-parameter values. Can be used one of the next metrics SORENSEN, JACCARD, FPB, TSS, KAPPA, AUC, and BOYCE. TSS is used as default.
#'
#' @importFrom dplyr %>% select starts_with bind_rows tibble group_by_at summarise across everything pull
#' @importFrom randomForest randomForest
#' @importFrom stats formula na.omit predict
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "randomForest" class object. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c colum names) and qualitative (f colum names) variables use for modeling.
#' \item performance: Hyper-parameters values and performance metric (see \code{\link{sdm_eval}}) for the best hyper-parameters combination.
#' \item hyper_performance: Performance metric (see \code{\link{sdm_eval}}) for each combination of the hyper-parameters.
#' \item data_ens: Predicted suitability for each test partition based on the best model. This database is used in \code{\link{fit_ensemble}}
#' }
#'
#' @export
#'
#' @seealso \code{\link{tune_gbm}}, \code{\link{tune_max}}, \code{\link{tune_net}}, and \code{\link{tune_svm}}.
#'
#' @examples
#' \dontrun{
#' data(abies_db)
#' abies_db
#'
#' # We will partition the data with the k-fold method
#'
#' abies_db2 <- part(
#'   data = abies_db,
#'   pr_ab = "pr_ab",
#'   bg_data = NULL,
#'   bg_a = NULL,
#'   method = c(method = "kfold", folds = 5)
#' )
#'
#' tune_grid <-
#'   expand.grid(mtry = seq(1, 7, 1))
#'
#' rf_t <-
#'   tune_raf(
#'     data = abies_db2,
#'     response = "pr_ab",
#'     predictors = c(
#'       "aet", "cwd", "tmin", "ppt_djf",
#'       "ppt_jja", "pH", "awc", "depth", "percent_clay"
#'     ),
#'     predictors_f = c("landform"),
#'     partition = ".part",
#'     grid = tune_grid,
#'     thr = "max_sens_spec",
#'     metric = "TSS",
#'   )
#'
#' # Outputs
#' rf_t$model
#' rf_t$predictors
#' rf_t$performance
#' rf_t$hyper_performance
#' rf_t$data_ens
#' }
#'
tune_raf <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           grid = NULL,
           thr = NULL,
           metric = "TSS") {
    . <- model <- TPR <- IMAE <- thr_value <- n_presences <- n_absences <- NULL
    variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))

    data <- data.frame(data)

    # Transform response variable as factor
    data[, response] <- as.factor(data[, response])

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
      grid <- data.frame(mtry = round(sqrt(nv)))
    }
    if (class(grid) == "character") {
      nv <- length(stats::na.omit(c(predictors, predictors_f)))
      if (grid == "defalut") {
        grid <- expand.grid(mtry = seq(2, nv, 1))
      }
    }

    # Test hyper-parameters names
    hyperp <- names(grid)
    if (!all(c("mtry") %in% hyperp)) {
      stop("Database used in 'grid' argument has to contain this columns for tunning: 'mtry'")
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
          try(
            mod[[ii]] <-
              randomForest::randomForest(
                formula1,
                data = train[[i]],
                mtry = grid$mtry[ii],
                ntree = 500,
                importance = FALSE
              )
          )
        }


        filt <- sapply(mod, function(x) length(class(x)) > 1)
        mod <- mod[filt]
        grid2 <- grid[filt, ]
        tnames <- apply(grid2, 1, function(x) paste(x, collapse = "_"))

        # Predict for presences absences data
        pred_test <-
          lapply(mod, function(x) {
            data.frame(
              pr_ab = test[[i]][, response],
              pred = stats::predict(
                x,
                newdata = test[[i]],
                type = "prob"
              )[, 2]
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
            ) %>% dplyr::tibble(model = "raf", .)
        }

        names(eval) <- tnames
        eval <- dplyr::bind_rows(eval, .id = "tnames")

        eval <-
          dplyr::tibble(dplyr::left_join(dplyr::mutate(grid2, tnames),
            eval,
            by = "tnames"
          )) %>%
          dplyr::select(-tnames)
        eval_partial[[i]] <- eval
      }

      # Create final database with parameter performance 1
      names(eval_partial) <- 1:np2
      eval_partial <- eval_partial %>%
        dplyr::bind_rows(., .id = "partition")
      eval_partial_list[[h]] <- eval_partial
    }

    # Create final database with parameter performance 2
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
    pred_test_ens <- fit_raf(
      data = data,
      response = response,
      predictors = predictors,
      predictors_f = predictors_f,
      partition = partition,
      thr = thr,
      fit_formula = formula1,
      mtry = best_tune$mtry
    )[["data_ens"]]

    # Fit final models with best settings
    set.seed(1)
    mod <-
      randomForest::randomForest(
        formula1,
        data = data,
        mtry = best_hyperp$mtry,
        ntree = 500,
        importance = TRUE
      )

    pred_test <- data.frame(
      pr_ab = data[, response],
      pred = stats::predict(
        mod,
        newdata = data,
        type = "prob"
      )[, 2]
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