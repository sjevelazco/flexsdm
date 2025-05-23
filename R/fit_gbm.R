#' Fit and validate Generalized Boosted Regression models
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous variables).
#' Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative
#' predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor
#' variables (e.g. formula(pr_ab ~ aet + ppt_jja + pH + awc + depth + landform)).
#' Note that the variables used here must be consistent with those used in
#' response, predictors, and predictors_f arguments. Default is NULL.
#' @param partition character. Column name with training and validation partition groups.
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1) needed for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The following threshold criteria are available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which the Jaccard index is the highest.
#'   \item max_sorensen: The threshold at which the Sorensen index is highest.
#'   \item max_fpb: The threshold at which FPB (F-measure on presence-background data) is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to sensitivity value. If a sensitivity value is not specified, the default used is 0.9
#'   }
#' If more than one threshold type is used they must be concatenated, e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function will use all thresholds if no threshold is specified.
#'
#' @param n_trees Integer specifying the total number of trees to fit.
#' This is equivalent to the number of iterations and the number of basis
#' functions in the additive expansion. Default is 100.
#' @param n_minobsinnode Integer specifying the minimum number of
#' observations in the terminal nodes of the trees. Note that this
#' is the actual number of observations, not the total weight.
#' The default value used is nrow(data)*0.5/4
#' @param shrinkage Numeric. This parameter applied to each tree in the
#' expansion. Also known as the learning rate or step-size reduction;
#' 0.001 to 0.1 usually works, but a smaller learning rate typically
#' requires more trees. Default is 0.1.
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "gbm" class object from gbm package. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: Performance metric (see \code{\link{sdm_eval}}).
#' Threshold dependent metrics are calculated based on the threshold specified in thr argument.
#' \item data_ens: Predicted suitability for each test partition based on the best model. This database is used in \code{\link{fit_ensemble}}
#' }
#'
#' @seealso \code{\link{fit_gam}}, \code{\link{fit_gau}}, \code{\link{fit_glm}},
#' \code{\link{fit_max}}, \code{\link{fit_net}}, \code{\link{fit_raf}}, and \code{\link{fit_svm}}.
#'
#' @export
#'
#' @importFrom dplyr %>% select all_of starts_with bind_rows summarise across everything
#' @importFrom gbm gbm predict.gbm
#' @importFrom stats formula sd
#'
#' @examples
#' \dontrun{
#' data("abies")
#'
#' # Using k-fold partition method
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 10)
#' )
#' abies2
#'
#' gbm_t1 <- fit_gbm(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
#' )
#' gbm_t1$model
#' gbm_t1$predictors
#' gbm_t1$performance
#' gbm_t1$data_ens
#'
#' # Using bootstrap partition method
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "boot", replicates = 10, proportion = 0.7)
#' )
#' abies2
#'
#' gbm_t2 <- fit_gbm(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("ppt_jja", "pH", "awc"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = "max_sens_spec"
#' )
#' gbm_t2
#' }
fit_gbm <- function(data,
                    response,
                    predictors,
                    predictors_f = NULL,
                    fit_formula = NULL,
                    partition,
                    thr = NULL,
                    n_trees = 100,
                    n_minobsinnode = as.integer(nrow(data) * 0.5 / 4),
                    shrinkage = 0.1) {
  . <- model <- TPR <- IMAE <- rnames <- thr_value <- n_presences <- n_absences <- NULL
  variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))

  data <- data.frame(data)
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


  # Fit models
  np <- ncol(data %>% dplyr::select(dplyr::starts_with(partition)))
  p_names <- names(data %>% dplyr::select(dplyr::starts_with(partition)))
  eval_partial_list <- list()
  pred_test_ens <- data %>%
    dplyr::select(dplyr::starts_with(partition)) %>%
    apply(., 2, unique) %>%
    data.frame() %>%
    as.list() %>%
    lapply(., function(x) {
      x <- stats::na.exclude(x)
      x[!(x %in% c("train-test", "test"))] %>% as.list()
    })

  for (h in 1:np) {
    message("Replica number: ", h, "/", np)

    out <- pre_tr_te(data, p_names, h)
    train <- out$train
    test <- out$test
    np2 <- out$np2
    rm(out)

    eval_partial <- as.list(rep(NA, np2))
    pred_test <- list()
    mod <- list()

    for (i in 1:np2) {
      tryCatch({
        message("Partition number: ", i, "/", np2)
        set.seed(1)
        mod[[i]] <-
          suppressMessages(
            gbm::gbm(
              formula1,
              data = train[[i]],
              distribution = "bernoulli",
              n.trees = n_trees,
              n.minobsinnode = n_minobsinnode,
              shrinkage = shrinkage,
              bag.fraction = 0.9
            )
          )

        pred_test <- data.frame(
          pr_ab = test[[i]][, response],
          pred = suppressMessages(
            gbm::predict.gbm(
              mod[[i]],
              newdata = test[[i]],
              type = "response",
              se.fit = FALSE
            )
          )
        )

        pred_test_ens[[h]][[i]] <- pred_test %>%
          dplyr::mutate(rnames = rownames(test[[i]]))

        # Validation of model
        eval <-
          sdm_eval(
            p = pred_test$pred[pred_test$pr_ab == 1],
            a = pred_test$pred[pred_test$pr_ab == 0],
            thr = thr
          )
        eval_partial[[i]] <- dplyr::tibble(model = "gbm", eval)
      })
    }

    # Create final database with parameter performance
    names(eval_partial) <- 1:np2
    eval_partial <-
      eval_partial[sapply(eval_partial, function(x) !is.null(dim(x)))] %>%
      dplyr::bind_rows(., .id = "partition")
    eval_partial_list[[h]] <- eval_partial
  }

  eval_partial <- eval_partial_list %>%
    dplyr::bind_rows(., .id = "replica")

  eval_final <- eval_partial %>%
    dplyr::group_by(model, threshold) %>%
    dplyr::summarise(dplyr::across(
      TPR:IMAE,
      list(mean = mean, sd = stats::sd)
    ), .groups = "drop")

  # Bind data for ensemble
  pred_test_ens <-
    lapply(pred_test_ens, function(x) {
      bind_rows(x, .id = "part")
    }) %>%
    bind_rows(., .id = "replicates") %>%
    dplyr::tibble() %>%
    dplyr::relocate(rnames)

  # Fit final models with best settings
  set.seed(1)
  suppressMessages(mod <-
    gbm::gbm(formula1,
      data = data,
      distribution = "bernoulli",
      n.trees = n_trees,
      n.minobsinnode = n_minobsinnode,
      shrinkage = shrinkage
    ))

  pred_test <- data.frame(
    pr_ab = data.frame(data)[, response],
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
    performance = dplyr::left_join(eval_final, threshold[1:4], by = "threshold") %>% dplyr::relocate(model, threshold, thr_value, n_presences, n_absences),
    data_ens = pred_test_ens
  )
  return(result)
}
