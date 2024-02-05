#' Fit and validate Maximum Entropy Models based on Ensemble of Small of Model approach
#'
#' @description This function constructs Maxent Models using the
#' Ensemble of Small Model (ESM) approach (Breiner et al., 2015, 2018).
#'
#' @param data data.frame. Database with the response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1)
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous variables). This function can only construct models with continuous variables, and does not allow categorical variables
#' Usage predictors = c("aet", "cwd", "tmin").
#' @param partition character. Column name with training and validation partition groups.
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The following threshold criteria are available:
#' \itemize{
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which Jaccard is the highest.
#'   \item max_sorensen: The threshold at which Sorensen is highest.
#'   \item max_fpb: The threshold at which FPB (F-measure on presence-background data) is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to sensitivity value. If no sensitivity value is specified, the default is 0.9
#'   }
#' If the user wants to include more than one threshold type, it is necessary concatenate threshold types, e.g., thr=c('max_sens_spec', 'max_jaccard'), or thr=c('max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('max_sens_spec', 'sensitivity'). Function will use all thresholds if no threshold is specified.
#' @param background data.frame. Database with response column only with 0 and predictors variables. All
#' column names must be consistent with data. Default NULL
#' @param clamp logical. It is set with TRUE, predictors and features are restricted to the range seen during model training.
#' @param classes character. A single feature of any combinations of them. Features are symbolized by letters: l (linear), q (quadratic), h (hinge), p (product), and t (threshold). Usage classes = "lpq". Default "default" (see details).
#' @param pred_type character. Type of response required available "link", "exponential", "cloglog" and "logistic". Default "cloglog"
#' @param regmult numeric. A constant to adjust regularization. Because ESM are used
#' for modeling species with few records default value is 2.5
#'
#' @details This method consists of creating bivariate models with all the pair-wise combinations
#' of predictors and perform an ensemble based on the average of suitability weighted by
#' Somers' D metric (D = 2 x (AUC -0.5)). ESM is recommended for modeling species with few occurrences.
#' This function does not allow categorical variables because the use of these types of variables
#' could be problematic when using with few occurrences. For further detail see
#' Breiner et al. (2015, 2018). This function use a default regularization multiplier
#' equal to 2.5 (see  Breiner et al., 2018)
#'
#' When the argument “classes” is set as default MaxEnt will use different features combination
#' depending of the number of presences (np) with the follow rule:
#' if np < 10 classes = "l",
#' if np between 10  and 15 classes = "lq",
#' if np between 15 and 80 classes = "lqh",
#' and if np >= 80 classes = "lqph"
#'
#' When presence-absence (or presence-pseudo-absence) data are used in data argument
#' in addition to background points, the function will fit models with presences and background
#' points and validate with presences and absences. This procedure makes maxent comparable to other
#' presences-absences models (e.g., random forest, support vector machine). If only presences and
#' background points data are used, function will fit and validate model with presences and
#' background data. If only presence-absences are used in data argument and without background,
#' function will fit model with the specified data (not recommended).
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item esm_model: A list with "maxnet" class object from maxnet package for each bivariate model. This object can be used
#' for predicting ensembles of small models with \code{\link{sdm_predict}} function.
#' \item predictors: A tibble with variables use for modeling.
#' \item performance: Performance metrics (see \code{\link{sdm_eval}}).
#' Those threshold dependent metric are calculated based on the threshold specified in the argument.
#' }
#'
#' @seealso \code{\link{esm_gam}}, \code{\link{esm_gau}}, \code{\link{esm_gbm}},
#' \code{\link{esm_glm}}, \code{\link{esm_net}}, and \code{\link{esm_svm}}.
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
#' data("backg")
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
#'   method = c(method = "rep_kfold", folds = 5, replicates = 5)
#' )
#' abies2
#'
#' set.seed(10)
#' backg2 <- backg %>%
#'   na.omit() %>%
#'   group_by(pr_ab) %>%
#'   dplyr::slice_sample(n = 100) %>%
#'   group_by()
#'
#' backg2 <- part_random(
#'   data = backg2,
#'   pr_ab = "pr_ab",
#'   method = c(method = "rep_kfold", folds = 5, replicates = 5)
#' )
#' backg2
#'
#' # Without threshold specification and with kfold
#' esm_max_t1 <- esm_max(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "cwd", "tmin", "ppt_djf", "ppt_jja", "pH", "awc", "depth"),
#'   partition = ".part",
#'   thr = NULL,
#'   background = backg2,
#'   clamp = TRUE,
#'   classes = "default",
#'   pred_type = "cloglog",
#'   regmult = 1
#' )
#'
#' esm_max_t1$esm_model # bivariate model
#' esm_max_t1$predictors
#' esm_max_t1$performance
#' }
esm_max <- function(data,
                    response,
                    predictors,
                    partition,
                    thr = NULL,
                    background = NULL,
                    clamp = TRUE,
                    classes = "default",
                    pred_type = "cloglog",
                    regmult = 2.5) {
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
      list_esm[[f]] <- fit_max(
        data = data,
        response = response,
        predictors = unlist(formula1[, f]),
        predictors_f = NULL,
        fit_formula = NULL,
        partition = partition,
        background = background,
        thr = thr,
        clamp = clamp,
        classes = classes,
        pred_type = pred_type,
        regmult = regmult
      )
    )
    utils::setTxtProgressBar(pb, which(1:ncol(formula1) == f))
  }
  close(pb)

  # Extract performance
  eval_esm <- lapply(list_esm, function(x) {
    x <- x$performance
    x$model <- "esm_max"
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
  filt <- mtrc > 0

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
    performance = dplyr::left_join(tibble(model = "esm_max", eval_final),
      threshold[1:4],
      by = "threshold"
    ) %>%
      dplyr::relocate(model, threshold, thr_value, n_presences, n_absences)
  )

  return(result)
}
