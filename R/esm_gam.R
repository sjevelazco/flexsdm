#' Fit and validate Generalized Additive Models based Ensemble of Small Model approach
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous variables).
#' Usage predictors = c("aet", "cwd", "tmin"). This function only can construct models with
#' continuous variables.
#' @param partition character. Column name with training and validation partition groups.
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The next threshold area available:
#' \itemize{
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which Jaccard is the highest.
#'   \item max_sorensen: The threshold at which Sorensen is highest.
#'   \item max_fpb: The threshold at which FPB is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to sensitivity value. If it is not specified a sensitivity values, function will use by default 0.9
#'   }
#' In the case of use more than one threshold type it is necessary concatenate threshold types, e.g., thr=c('max_sens_spec', 'max_jaccard'), or thr=c('max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('max_sens_spec', 'sensitivity'). Function will use all thresholds if no threshold is specified
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A list with "Gam" class object for each bivariate model. This object can be used for predicting ensemble of small model with \code{\link{sdm_predict}} function.
#' \item predictors: A tibble with variables use for modeling.
#' \item performance: Performance metric (see \code{\link{sdm_eval}}).
#' Those threshold dependent metric are calculated based on the threshold specified in thr argument.
#' }
#'
#' @export
#'
#' @importFrom dplyr bind_rows filter group_by distinct pull mutate inner_join select starts_with bind_cols summarise across left_join relocate
#' @importFrom stats sd
#' @importFrom utils combn txtProgressBar setTxtProgressBar
#'
#' @examples
#' \dontrun{
#' data("abies_db")
#' require(dplyr)
#'
#' # Using k-fold partition method
#' set.seed(10)
#' abies_db2 <- abies_db %>% na.omit %>%
#'   group_by(pr_ab) %>%
#'   dplyr::slice_sample(n = 10) %>%
#'   group_by()
#'
#' abies_db2 <- part(
#'   data = abies_db2,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 3)
#' )
#' abies_db2
#'
#' # Without thrshold specification and with kfold
#' esm_gam_t1 <- esm_gam(
#'   data = abies_db2,
#'   response = "pr_ab",
#'   predictors = c("aet", "cwd", "tmin", "ppt_djf", "ppt_jja", "pH", "awc", "depth", "percent_clay"),
#'   partition = ".part",
#'   thr = NULL
#' )
#'
#' esm_gam_t1$model
#' esm_gam_t1$predictors
#' esm_gam_t1$performance
#' }
#'
esm_gam <- function(data,
                    response,
                    predictors,
                    partition,
                    thr = NULL) {
  . <- model <- TPR <- IMAE <- rnames <- thr_value <- n_presences <- n_absences <- NULL
  variables <- dplyr::bind_rows(c(c = predictors))

  # Formula
  formula1 <- utils::combn(variables, 2)
  nms <- apply(utils::combn(variables, 2), 2, function(x) paste(x, collapse = "_")) %>%
    paste0(".", .)

  # Fit models
  eval_esm <- list()
  list_esm <- list()
  pb <- utils::txtProgressBar(min = 0, max = ncol(formula1), style = 3)
  for (f in 1:ncol(formula1)) {
    suppressMessages(
      list_esm[[f]] <- fit_gam(
        data = data,
        response = response,
        predictors = unlist(formula1[, f]),
        predictors_f = NULL,
        partition = partition,
        thr = thr
      )
    )
    utils::setTxtProgressBar(pb, which(1:ncol(formula1) == f))
  }
  close(pb)

  # Extract performance
  eval_esm <- lapply(list_esm, function(x) {
    x <- x$performance
    x$model <- "esm_gam"
    x
  })
  names(eval_esm) <- nms
  eval_esm <- eval_esm %>%
    dplyr::bind_rows(., .id = "variables") %>%
    dplyr::filter(threshold != "lpt")

  # Calculate Somers' metric and remove small models with bad performance (AUC<0.5)
  mtrc <- eval_esm %>%
    dplyr::group_by(variables) %>%
    dplyr::distinct(AUC_mean) %>%
    dplyr::pull()
  D <- 2 * (mtrc - 0.5) # Somers'D
  filt <- mtrc >= 0.5

  # Filter data based on Somers<0.5
  D <- D[filt]
  list_esm <- list_esm[filt]
  nms <- nms[filt]

  # Perform weighted ensemble
  data_ens <- sapply(list_esm, function(x) {
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
          data_ens[[i]],
          by = c("rnames", "replicates", "part", "pr_ab")
        )
    }
  }
  rm(data_ens)

  #### Extract predicted suitability of each model
  values <- data_ens2 %>%
    dplyr::select(dplyr::starts_with("."))

  #### Remove suitability values from data_ens2
  data_ens2 <- data_ens2 %>% dplyr::select(-dplyr::starts_with("."))

  pred <- mapply(function(x, v) {
    (x * v)
  }, values, D, SIMPLIFY = TRUE) %>%
    apply(., 1, function(x) {
      sum(x, na.rm = TRUE)
    }) / sum(D)

  pred_test <- dplyr::bind_cols(data_ens2, pred = pred)

  threshold <- sdm_eval(
    p = pred_test$pred[pred_test$pr_ab == 1],
    a = pred_test$pred[pred_test$pr_ab == 0],
    thr = thr
  )

  eval_esm <- eval_esm %>% dplyr::select(model:n_absences, dplyr::ends_with("_mean"))
  names(eval_esm) <- gsub("_mean", "", names(eval_esm))

  eval_final <- eval_esm %>%
    dplyr::group_by(model, threshold) %>%
    dplyr::summarise(dplyr::across(
      TPR:IMAE,
      list(mean = mean, sd = stats::sd)
    ), .groups = "drop")

  mod <- lapply(list_esm, function(x) x$mode)
  names(mod) <- gsub("[$.]", "", nms)

  result <- list(
    model = mod,
    predictors = variables,
    performance = dplyr::left_join(eval_final, threshold[1:4], by = "threshold") %>%
      dplyr::relocate(model, threshold, thr_value, n_presences, n_absences)
  )

  return(result)
}
