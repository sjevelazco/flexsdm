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
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "Gam" "glm" "lm"  class object. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c colum names) and qualitative (f colum names) variables use for modeling.
#' \item performance: Performance metric (see \code{\link{sdm_eval}}).
#' Those threshold dependent metric are calculated based on the threshold specified in thr argument .
#' \item data_ens: Predicted suitability for each test partition. This database is used in \code{\link{fit_ensemble}}
#' }
#'
#' @export
#'
#' @importFrom dplyr %>% select all_of starts_with bind_rows summarise across everything
#' @importFrom gam gam predict.Gam s
#' @importFrom stats formula sd
#'
#' @examples
#' \dontrun{
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
  nms <- apply(utils::combn(variables, 2), 2, function(x) paste(x, collapse = "_"))

  # Fit models
  eval_esm <- list()
  list_esm <- list()
  for (f in 1:ncol(formula1)) {
    message("Small model number: ", f)
    list_esm[[f]] <- fit_gam(
      data = data,
      response = response,
      predictors = unlist(formula1[, f]),
      predictors_f = NULL,
      partition = partition,
      thr = thr
    )
  }

  # Extract performance
  eval_esm <- lapply(list_esm, function(x) {
    x <- x$performance
    x$model <- 'esm_gam'
    x
  }) %>%
    dplyr::bind_rows()

  # Extract final models
  mod <- lapply(list_esm, function(x) x$mode)

  # Extract prediction for used occurrences database
  pred_test_ens <- lapply(list_esm, function(x) x$data_ens)

  # Calculate Somers' metric and remove small models with bad performance (AUC<0.5)
  mtrc <- eval_esm[, "AUC_mean"][[1]]
  D <- 2 * (mtrc - 0.5) #Somers'D
  filt <- mtrc>=0.5
  D <- D[filt]
  pred_test_ens <- pred_test_ens[filt]
  mod <- mod[filt]
  names(mod) <- nms[filt]
  pred_test_ens <-

  # Perform weighted ensemble
    data_ens <- sapply(mod, function(x) {
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


  data_ens <- sapply(mod, function(x) {
    x$data_ens
  })

  data_ens <- mapply(function(x, nms) {
    colnames(x)[colnames(x) %in% "pred"] <- nms
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

  data_ens <- sapply(pred_test_ens, function(x) {
    x["pred"]
  })
  names(data_ens) <- nms
  data_ens <- bind_cols(data_ens)
  pred <- mapply(function(x, v) {
    x * v
  }, data_ens, D) %>%
    apply(., 1, function(x) {
      sum(x, na.rm = TRUE)
    }) / sum(D)


  pred_test <- tibble(pred_test[[1]][1], pred)

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
