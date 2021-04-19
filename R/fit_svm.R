#' Fit and validate Support Vector Machine models
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous or discrete variables).
#' Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative
#' predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param partition character. Column name with training and validation partition groups.
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1).
#' It is useful for threshold-dependent performance metrics.
#' It is possible to use more than one threshold type. It is necessary to provide a
#' vector for this argument. The next threshold area available:
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
#' @param fit_formula formula. A formula object with response and predictor
#' variables (e.g. forumla(pr_ab ~ aet + ppt_jja + pH + awc + depth + landform)).
#' Note that the variables used here must be consistent with those used in
#' response, predictors, and predictors_f arguments
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "ksvm" class object. This object can be used for predicting.
#' \item performance: Performance metric (see \code{\link{enm_eval}}).
#' Those threshold dependent metric are calculated based on the threshold specified in thr argument .
#' \item selected_threshold: Value of the threshold selected.
#' \item threshold_table: Value of all threshold.
#' }
#'
#' @export
#'
#' @importFrom dismo predict
#' @importFrom dplyr select all_of starts_with bind_rows group_by summarise across everything
#' @importFrom kernlab ksvm
#' @importFrom stats formula sd
#'
#' @examples
#' \dontrun{
#' data('abies_db')
#'
#' # Using KFOLD partition method
#' abies_db2 <- data_part(
#'   data = abies_db,
#'   p_a = 'pr_ab',
#'   method = c(method = "KFOLD", folds = 10)
#' )
#' abies_db2
#'
#' svm_t1 <- fit_svm(data = abies_db2,
#'                   response = "pr_ab",
#'                   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'                   predictors_f = c("landform"),
#'                   partition = ".part",
#'                   thr = c("MAX_TSS", "EQUAL_SENS_SPEC", "MAX_SORENSEN"),
#'                   fit_formula = NULL)
#'
#' svm_t1$model
#' svm_t1$performance
#' svm_t1$selected_threshold
#' svm_t1$threshold_table
#'
#' # Using BOOTS partition method and only with presence-absence
#' # and get performance for several method
#' abies_db2 <- data_part(
#'   data = abies_db,
#'   p_a = 'pr_ab',
#'   method = c(method = "BOOT", replicates=10,  proportion=0.7)
#' )
#' abies_db2
#'
#' svm_t2 <- fit_svm(data = abies_db2,
#'                   response = "pr_ab",
#'                   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'                   predictors_f = c("landform"),
#'                   partition = ".part",
#'                   thr = c("MAX_TSS", "EQUAL_SENS_SPEC", "MAX_SORENSEN"),
#'                   fit_formula = NULL)
#' svm_t2
#' }
#'
fit_svm <- function(data,
                   response,
                   predictors,
                   predictors_f = NULL,
                   partition,
                   thr = NULL,
                   fit_formula = NULL,
                   ...) {
  data <- data.frame(data)

  # Transform response variable as factor
  data[, response] <- as.factor(data[, response])

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

  # Formula
  if (is.null(fit_formula)) {
    formula1 <- stats::formula(paste(
      response, "~",
      paste(c(predictors, predictors_f), collapse = " + ")
    ))
    message(
      "Formula used for model fitting:\n",
      Reduce(paste, deparse(formula1)) %>% gsub(paste("  ", "   ", collapse = "|"), " ", .)
    )
  } else {
    formula1 <- fit_formula
    message(
      "Formula used for model fitting:\n",
      Reduce(paste, deparse(formula1)) %>% gsub(paste("  ", "   ", collapse = "|"), " ", .)
    )
  }


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
    pred_test <- list()
    mod <- list()

    for (i in 1:np2) {
      message("Partition number: ", i, "/", np2)
      set.seed(1)
      try(
        mod[[i]] <-
          kernlab::ksvm(
            formula1,
            data = train[[i]],
            type = "C-svc",
            kernel = "rbfdot",
            # kpar = list(sigma = grid$sigma[ii]),
            # C = grid$C[ii],
            prob.model = TRUE
          )
      )


      pred_test <- try(data.frame(
        pr_ab = test[[i]][, response],
        pred = suppressMessages(
          dismo::predict(
            mod[[i]],
            newdata = test[[i]],
            type = "prob",
          )[,2]
        )
      ))

      # Validation of model
      eval <-
        enm_eval(
          p = pred_test$pred[pred_test$pr_ab == 1],
          a = pred_test$pred[pred_test$pr_ab == 0],
          thr = thr
        )
      eval_partial[[i]] <- eval$selected_threshold
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
    dplyr::group_by(threshold) %>%
    dplyr::select(-c(replica:partition, values:n_absences)) %>%
    dplyr::summarise(dplyr::across(
      dplyr::everything(),
      list(mean = mean, sd = stats::sd)
    ), .groups = "drop")

  # Fit final models with best settings
  set.seed(1)
  suppressMessages(mod <-
                     kernlab::ksvm(
                       formula1,
                       data = data,
                       type = "C-svc",
                       kernel = "rbfdot",
                       # kpar = list(sigma = grid$sigma[ii]),
                       # C = grid$C[ii],
                       prob.model = TRUE
                     )
  )

  pred_test <- data.frame(
    pr_ab = data[, response],
    pred = dismo::predict(
      mod,
      newdata = data,
      type = "prob"
    )[,2]
  )

  threshold <- enm_eval(
    p = pred_test$pred[pred_test$pr_ab == 1],
    a = pred_test$pred[pred_test$pr_ab == 0],
    thr = thr
  )

  result <- list(
    model = mod,
    performance = eval_final,
    selected_threshold = threshold[[1]] %>% dplyr::select(threshold:values),
    threshold_table = threshold[[2]] %>% dplyr::select(threshold:values)
  )
  return(result)
}