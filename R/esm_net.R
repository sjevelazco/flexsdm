#' Fit and validate Neural Networks based on Ensembles of Small of Models approach
#'
#' @description This function constructs Neural Networks using the
#' Ensembles of Small Models (ESM) approach (Breiner et al., 2015, 2018).
#'
#' @param data data.frame. Database with the response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1)
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous variables). This function only can construct models with continuous variables and does not allow categorical variables.
#' Usage predictors = c("aet", "cwd", "tmin").
#' @param partition character. Column name with training and validation partition groups.
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The following threshold criteria are available:
#' \itemize{
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which the Jaccard index is the highest.
#'   \item max_sorensen: The threshold at which the Sorensen index is highest.
#'   \item max_fpb: The threshold at which FPB (F-measure on presence-background data) is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to sensitivity value. If a sensitivity values is not specified, the default used is 0.9
#'   }
#' If the user wants to include more than one threshold type, it is necessary concatenate threshold types, e.g., thr=c('max_sens_spec', 'max_jaccard'), or thr=c('max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('max_sens_spec', 'sensitivity'). Function will use all thresholds if no threshold is specified
#' @param size numeric. Number of units in the hidden layer. Can be zero if there are skip-layer units. Default 2.
#' @param decay numeric. Parameter for weight decay. Default 0.
#'
#' @details This method consists of creating bivariate models with all the pair-wise combinations
#' of predictors and perform an ensemble based on the average of suitability weighted by
#' Somers' D metric (D = 2 x (AUC -0.5)). ESM is recommended for modeling species with few occurrences.
#' This function does not allow categorical variables because the use of these types of variables
#' could be problematic when using with few occurrences. Further detail see
#' Breiner et al. (2015, 2018).
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item esm_model: A list with "nnet" class object from nnet package for each bivariate model. This object can be used
#' for predicting ensemble of small model with \code{\link{sdm_predict}} function.
#' \item predictors: A tibble with variables use for modeling.
#' \item performance: Performance metric (see \code{\link{sdm_eval}}).
#' Those threshold dependent metric are calculated based on the threshold specified in thr argument.
#' }
#'
#' @seealso \code{\link{esm_gam}}, \code{\link{esm_gau}}, \code{\link{esm_gbm}},
#' \code{\link{esm_glm}}, \code{\link{esm_max}}, and \code{\link{esm_svm}}.
#'
#' @export
#'
#' @references
#' \itemize{
#' \item Breiner, F. T., Guisan, A., Bergamini, A., & Nobis, M. P. (2015). Overcoming limitations of modelling rare species by using ensembles of small models. Methods in Ecology and Evolution, 6(10), 1210-218. https://doi.org/10.1111/2041-210X.12403
#' \item Breiner, F. T., Nobis, M. P., Bergamini, A., & Guisan, A. (2018). Optimizing ensembles of small models for predicting the distribution of species with few occurrences. Methods in Ecology and Evolution, 9(4), 802-808. https://doi.org/10.1111/2041-210X.12957
#' }
#'
#' @importFrom dplyr bind_rows filter group_by distinct pull mutate inner_join select starts_with bind_cols summarise across left_join relocate
#' @importFrom stats sd
#' @importFrom utils combn txtProgressBar setTxtProgressBar
#'
#' @examples
#' \dontrun{
#' data("abies")
#' require(dplyr)
#'
#' # Using k-fold partition method
#' set.seed(10)
#' abies2 <- abies %>%
#'   na.omit() %>%
#'   group_by(pr_ab) %>%
#'   dplyr::slice_sample(n = 10) %>%
#'   group_by()
#'
#' abies2 <- part_random(
#'   data = abies2,
#'   pr_ab = "pr_ab",
#'   method = c(method = "rep_kfold", folds = 3, replicates = 5)
#' )
#' abies2
#'
#' # Without threshold specification and with kfold
#' esm_net_t1 <- esm_net(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "cwd", "tmin", "ppt_djf", "ppt_jja", "pH", "awc", "depth"),
#'   partition = ".part",
#'   thr = NULL
#' )
#'
#' esm_net_t1$esm_model # bivariate model
#' esm_net_t1$predictors
#' esm_net_t1$performance
#' }
esm_net <- function(data,
                    response,
                    predictors,
                    partition,
                    thr = NULL,
                    size = 2, # TODO find a formula to calculate default value for this argument
                    decay = 0) {
  . <- part <- model <- TPR <- IMAE <- rnames <- thr_value <- n_presences <- n_absences <- AUC_mean <- pr_ab <- NULL
  variables <- dplyr::bind_rows(c(c = predictors))

  # N of predictor requirement
  if (length(predictors) <= 2) {
    stop("The 'esm_' family function should be used to build models with more than 2 predictors, use the 'fit_' or 'tune_' family functions instead")
  }

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
      list_esm[[f]] <- fit_net(
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
    x$model <- "esm_net"
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
  filt <- D > 0

  if (sum(filt) == 0) {
    message("None bivariate model had Somer's D > 0.5. Try with another esm_* function. NA will be returned")
    return(NA)
  }

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

  if (length(data_ens) > 1) {
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
  } else {
    data_ens2 <- data_ens[[1]]
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

  pred_test <- dplyr::bind_cols(data_ens2, pred = pred) # This dataset will be use to calculate
  # esm peformance

  # Validate ensemble base on weighted average suitability, split data by
  # replicates and partition
  testlist <- pred_test %>% dplyr::distinct(replicates, part)
  replicates <- as.list(unique(pred_test$replicates))
  names(replicates) <- unique(pred_test$replicates)

  for (r in unique(pred_test$replicates)) {
    x0 <- pred_test %>% dplyr::filter(replicates == r)
    replicates[[r]] <- lapply(
      split(x0, x0$part),
      function(x) {
        sdm_eval(
          p = x$pred[x$pr_ab == 1],
          a = x$pred[x$pr_ab == 0],
          thr = thr
        )
      }
    ) %>%
      dplyr::bind_rows(., .id = "part")
  }

  eval_esm <- dplyr::bind_rows(replicates, .id = "replicates")


  eval_final <- eval_esm %>%
    dplyr::group_by(threshold) %>%
    dplyr::summarise(dplyr::across(
      TPR:IMAE,
      list(mean = mean, sd = stats::sd)
    ), .groups = "drop")

  # Calculate final threshold
  threshold <- sdm_eval(
    p = pred_test$pred[pred_test$pr_ab == 1],
    a = pred_test$pred[pred_test$pr_ab == 0],
    thr = thr
  )

  # List of models used for prediction
  mod <- lapply(list_esm, function(x) x$mode)
  names(mod) <- D

  result <- list(
    esm_model = mod,
    predictors = variables,
    performance = dplyr::left_join(tibble(model = "esm_net", eval_final),
      threshold[1:4],
      by = "threshold"
    ) %>%
      dplyr::relocate(model, threshold, thr_value, n_presences, n_absences)
  )

  return(result)
}
