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

  data <- data.frame(data)
  data <- data %>%
    dplyr::select(
      dplyr::all_of(response),
      dplyr::all_of(predictors),
      dplyr::starts_with(partition)
    )
  data <- data.frame(data)

  # Remove NAs
  complete_vec <- stats::complete.cases(data[, c(response, unlist(variables))])
  if (sum(!complete_vec) > 0) {
    message(sum(!complete_vec), " rows were excluded from database because NAs were found")
    data <- data %>% dplyr::filter(complete_vec)
  }
  rm(complete_vec)

  # Formula
  formula1 <- utils::combn(variables, 2) %>%
    apply(., 2, function(x) {
      formula1 <- c(paste("s(", x, ")", collapse = " + ", sep = ""))
      formula1 <- stats::formula(paste(response, "~", formula1))
      return(formula1)
    })

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


  eval_esm <- as.list(rep(NA, length(formula1)))
  names(eval_esm) <- apply(utils::combn(variables, 2), 2, function(x) paste(x[1], x[2], sep='||'))
  for(f in 1:length(formula1)){
    message("Small model number: ", f, "/", np)
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
          suppressWarnings(mod[[i]] <-
                             gam::gam(formula1[[f]],
                                      data = train[[i]],
                                      family = "binomial"))


          # Predict for presences absences data
          if (!is.null(predictors_f)) {
            for (fi in 1:length(predictors_f)) {
              lev <- as.character(unique(mod[[i]]$data[, predictors_f[fi]]))
              lev_filt <- test[[i]][, predictors_f[fi]] %in% lev
              test[[i]] <- test[[i]][lev_filt,]
            }
          }

          pred_test <- data.frame(pr_ab = test[[i]][, response],
                                  pred = suppressWarnings(
                                    gam::predict.Gam(
                                      mod[[i]],
                                      newdata = test[[i]],
                                      type = "response",
                                      se.fit = FALSE
                                    )
                                  ))

          pred_test_ens[[h]][[i]] <- pred_test %>%
            dplyr::mutate(rnames = rownames(.))

          # Validation of model
          eval <-
            sdm_eval(p = pred_test$pred[pred_test$pr_ab == 1],
                     a = pred_test$pred[pred_test$pr_ab == 0],
                     thr = thr)
          # eval_partial[[i]] <- dplyr::tibble(model = "gam", eval)
        },
        error = function(cond) {
          message("It was not possible to fit this model")
        })
      }
      names(eval_partial) <- 1:np2

      # Create final database with parameter performance
      eval_partial <- eval_partial[sapply(eval_partial, function(x) !is.null(dim(x)))] %>%
        dplyr::bind_rows(., .id = "partition")
      eval_partial_list[[h]] <- eval_partial
    }

    eval_esm[[f]] <- eval_partial_list %>%
      dplyr::bind_rows(., .id = "replica")
  }

  eval_esm <- eval_esm %>%
    dplyr::bind_rows(., .id = "esm")

  eval_esm <- eval_esm %>%
    dplyr::group_by(esm) %>%
    dplyr::summarise(dplyr::across(
      TPR:IMAE,
      list(mean = mean, sd = stats::sd)
    ), .groups = "drop") %>%
    dplyr::tibble(model = "esm_gam",.)


  # Fit final models with best settings
  suppressWarnings(mod <- lapply(formula1, function(x)
                     gam::gam(x,
                              data = data,
                              family = "binomial"
                     ))
  )

  pred_test <- lapply(mod, function(x)
    dplyr::tibble(
    pr_ab = data[, response],
    pred = suppressMessages(gam::predict.Gam(
      x,
      newdata = data,
      type = "response"
    ))
  ))

  # Perform weighted ensemble
  fun1 <- function(x1, x2){
    x1[,2]*x2[1,]
  }
  pred_test

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
