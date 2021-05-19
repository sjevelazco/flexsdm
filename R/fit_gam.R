#' Fit and validate Generalized Additive Models
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous or discrete variables).
#' Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative
#' predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param partition character. Column name with training and validation partition groups.
##' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The next threshold area available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission. Usage thr=c(type='lpt').
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal (aka threshold that maximizes the TSS).
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest.
#'   Usage thr=c(type='max_sens_spec').
#'   \item max_kappa: The threshold at which Kappa is the highest ("max kappa"). Usage thr=c(type='max_kappa').
#'   \item max_jaccard: The threshold at which Jaccard is the highest. Usage thr=c(type='max_jaccard').
#'   \item max_sorensen: The threshold at which Sorensen is highest. Usage thr=c(type='max_sorensen').
#'   \item max_fpb: The threshold at which FPB is highest. Usage thr=c(type='max_fpb').
#'   \item specific: A threshold value specified by user. Usage thr=c(type='specific', sens='0.6'). 'sens' refers to models will be binarized using this suitability value.
#'   }
#' @param fit_formula formula. A formula object with response and predictor
#' variables (e.g. formula(pr_ab ~ aet + ppt_jja + pH + awc + depth + landform)).
#' Note that the variables used here must be consistent with those used in
#' response, predictors, and predictors_f arguments
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "Gam" "glm" "lm"  class object. This object can be used for predicting.
#' \item predictors: A character with quantitative (elements names with c) and qualitative (elements names with f) variables use for modeling.
#' \item performance: Performance metric (see \code{\link{sdm_eval}}).
#' Those threshold dependent metric are calculated based on the threshold specified in thr argument .
#' \item selected_thresholds: Value of the threshold selected.
#' \item threshold_table: Value of all threshold.
#' }
#'
#' @export
#'
#' @importFrom dplyr select all_of starts_with bind_rows summarise across everything
#' @importFrom gam gam predict.Gam s
#' @importFrom stats formula sd
#'
#' @examples
#' \dontrun{
#' data("abies_db")
#' require(gam)
#'
#' # Using k-fold partition method
#' abies_db2 <- data_part(
#'   data = abies_db,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 10)
#' )
#' abies_db2
#'
#' gam_t1 <- fit_gam(
#'   data = abies_db2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = "max_sens_spec"
#' )
#' gam_t1$model
#' gam_t1$performance
#' gam_t1$selected_thresholds
#' gam_t1$all_thresholds
#'
#' # Using our own formula
#' gam_t2 <- fit_gam(
#'   data = abies_db2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = "max_sens_spec",
#'   fit_formula = stats::formula(pr_ab ~ s(aet, df = 4) +
#'     s(ppt_jja, df = 3) +
#'     s(pH, df = 3) + landform)
#' )
#'
#' gam_t2$model
#' gam_t2$performance %>% dplyr::select(ends_with("_mean"))
#' gam_t2$selected_thresholds
#' gam_t2$threshold_table
#'
#' # Using repeated k-fold partition method
#' abies_db2 <- data_part(
#'   data = abies_db,
#'   pr_ab = "pr_ab",
#'   method = c(method = "rep_kfold", folds = 10, replicates = 10)
#' )
#' abies_db2
#'
#' gam_t3 <- fit_gam(
#'   data = abies_db2,
#'   response = "pr_ab",
#'   predictors = c("ppt_jja", "pH", "awc"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = "max_sens_spec",
#'   poly = 3,
#'   inter_order = 2
#' )
#' gam_t3$model
#' }
#'
fit_gam <- function(data,
                    response,
                    predictors,
                    predictors_f = NULL,
                    partition,
                    thr = NULL,
                    fit_formula = NULL,
                    ...) {
  variables <- c(c = predictors, f = predictors_f)

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
  data <- rm_na(x = data)

  # Formula
  if (is.null(fit_formula)) {
    formula1 <-
      paste(c(
        paste("s(", predictors, ")", collapse = " + ", sep = ""),
        predictors_f
      ), collapse = " + ")
    formula1 <- stats::formula(paste(
      response, "~", formula1
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
    lapply(., as.list)

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
      try(suppressWarnings(mod[[i]] <-
        gam::gam(formula1,
          data = train[[i]],
          family = "binomial"
        )))


      # Predict for presences absences data
      if (!is.null(predictors_f)) {
        for (fi in 1:length(predictors_f)) {
          lev <- as.character(unique(mod[[i]]$data[, predictors_f[fi]]))
          lev_filt <- test[[i]][, predictors_f[fi]] %in% lev
          test[[i]] <- test[[i]][lev_filt, ]
        }
      }

      pred_test <- try(data.frame(
        pr_ab = test[[i]][, response],
        pred = suppressWarnings(
          gam::predict.Gam(
            mod[[i]],
            newdata = test[[i]],
            type = "response",
            se.fit = FALSE
          )
        )
      ))

      pred_test_ens[[h]][[i]] <- pred_test %>%
        dplyr::mutate(rnames=rownames(.))

      # Validation of model
      eval <-
        sdm_eval(
          p = pred_test$pred[pred_test$pr_ab == 1],
          a = pred_test$pred[pred_test$pr_ab == 0],
          thr = thr
        )
      if (is.null(thr)) {
        eval_partial[[i]] <- eval$all_thresholds
      } else {
        eval_partial[[i]] <- eval$selected_thresholds
      }
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

  # Bind data for ensemble
  pred_test_ens <-
    lapply(pred_test_ens, function(x)
      bind_rows(x, .id = 'part')) %>%
    bind_rows(., .id = 'replicates') %>% dplyr::tibble() %>%
    dplyr::relocate(rnames)


  # Fit final models with best settings
  suppressWarnings(mod <-
    gam::gam(formula1,
      data = data,
      family = "binomial"
    ))

  pred_test <- data.frame(
    pr_ab = data[, response],
    pred = suppressMessages(gam::predict.Gam(
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

  if (!is.null(thr)) {
    st <- threshold$selected_thresholds
  } else {
    st <- threshold$all_thresholds
  }

  result <- list(
    model = mod,
    predictors = variables,
    performance = eval_final,
    selected_thresholds = st %>% dplyr::select(threshold:values),
    all_thresholds = threshold$all_thresholds %>% dplyr::select(threshold:values),
    data_ens = pred_test_ens
  )
  return(result)
}
