#' Fit and validate Support Vector Machine models
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous or discrete variables).
#' Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative
#' predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor
#' variables (e.g. formula(pr_ab ~ aet + ppt_jja + pH + awc + depth + landform)).
#' Note that the variables used here must be consistent with those used in
#' response, predictors, and predictors_f arguments
#' @param partition character. Column name with training and validation partition groups.
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
#'
#' @param sigma numeric. Inverse kernel width for the Radial Basis kernel function "rbfdot". Default "automatic".
#' @param C numeric. Cost of constraints violation this is the 'C'-constant of the regularization term in the Lagrange formulation. Default 1
#' @param ...
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "ksvm" class object. This object can be used for predicting.
#' \item predictors: A character with quantitative (elements names with c) and qualitative (elements names with f) variables use for modeling.
#' \item performance: Performance metric (see \code{\link{sdm_eval}}).
#' Those threshold dependent metric are calculated based on the threshold specified in thr argument .
#' \item selected_thresholds: Value of the threshold selected.
#' \item all_thresholds: Value of all threshold.
#' }
#'
#' @details This function constructs 'C-svc' classification type and uses
#' Radial Basis kernel "Gaussian" function (rbfdot). See details details in \link[kernlab]{ksvm}.
#' @export
#'
#' @importFrom dplyr select all_of starts_with bind_rows group_by summarise across everything
#' @importFrom kernlab ksvm predict
#' @importFrom stats formula sd
#'
#' @examples
#' \dontrun{
#' data("abies_db")
#'
#' # Using k-fold partition method
#' abies_db2 <- data_part(
#'   data = abies_db,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 10)
#' )
#' abies_db2
#'
#' svm_t1 <- fit_svm(
#'   data = abies_db2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = c("max_sens_spec", "equal_sens_spec", "mas_sorensen"),
#'   fit_formula = NULL
#' )
#'
#' svm_t1$model
#' svm_t1$performance
#' svm_t1$selected_thresholds
#' svm_t1$all_thresholds
#'
#' # Using bootstrap partition method and only with presence-absence
#' # and get performance for several method
#' abies_db2 <- data_part(
#'   data = abies_db,
#'   pr_ab = "pr_ab",
#'   method = c(method = "boot", replicates = 10, proportion = 0.7)
#' )
#' abies_db2
#'
#' svm_t2 <- fit_svm(
#'   data = abies_db2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = c("max_sens_spec", "equal_sens_spec", "mas_sorensen"),
#'   fit_formula = NULL
#' )
#' svm_t2
#' }
#'
fit_svm <- function(data,
                    response,
                    predictors,
                    predictors_f = NULL,
                    fit_formula = NULL,
                    partition,
                    thr = NULL,
                    sigma = "automatic",
                    C = 1,
                    ...) {
  variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))

  data <- data.frame(data)

  # Transform response variable as factor
  data[, response] <- as.factor(data[, response])

  if (is.null(predictors_f)) {
    data <- data %>%
      dplyr::select(
        dplyr::all_of(response),
        dplyr::all_of(predictors),
        dplyr::starts_with(partition)
      )
    data <- data.frame(data)
  } else {
    data <- data %>%
      dplyr::select(
        dplyr::all_of(response),
        dplyr::all_of(predictors),
        dplyr::all_of(predictors_f),
        dplyr::starts_with(partition)
      )
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
      paste(c(
        predictors, predictors_f
      ), collapse = " + ")
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
  p_names <-
    names(data %>% dplyr::select(dplyr::starts_with(partition)))
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

    eval_partial <- list()
    pred_test <- list()
    mod <- list()

    for (i in 1:np2) {
      message("Partition number: ", i, "/", np2)
      tryCatch(
        {
          if (sigma == "automatic") {
            kpar_ <- "automatic"
          } else {
            kpar_ <- list(sigma = sigma)
          }

          set.seed(1)
          mod[[i]] <-
            suppressMessages(
              kernlab::ksvm(
                formula1,
                data = train[[i]],
                type = "C-svc",
                kernel = "rbfdot",
                kpar = kpar_,
                C = C,
                prob.model = TRUE
              )
            )

          pred_test <- data.frame(
            pr_ab = test[[i]][, response],
            pred = suppressMessages(kernlab::predict(
              mod[[i]],
              newdata = test[[i]],
              type = "prob",
            )[, 2])
          )

          pred_test_ens[[h]][[i]] <- pred_test %>%
            dplyr::mutate(rnames = rownames(.))

          # Validation of model
          eval <-
            sdm_eval(
              p = pred_test$pred[pred_test$pr_ab == 1],
              a = pred_test$pred[pred_test$pr_ab == 0],
              thr = thr
            )

          eval_partial[[i]] <- dplyr::tibble(model = "svm", eval)

          names(eval_partial) <- i
        },
        error = function(cond) {
          message("It was not possible to fit this model")
        }
      )
    }

    # Create final database with parameter performance
    eval_partial <- eval_partial %>%
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
  if (sigma == "automatic") {
    kpar_ <- "automatic"
  } else {
    kpar_ <- list(sigma = sigma)
  }
  set.seed(1)
  suppressMessages(
    mod <-
      kernlab::ksvm(
        formula1,
        data = data,
        type = "C-svc",
        kernel = "rbfdot",
        kpar = kpar_,
        C = C,
        prob.model = TRUE
      )
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
