#' Fit and validate Gaussian Process models
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous variables).
#' Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative
#' predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param partition character. Column name with training and validation partition groups.
#' If partition = NULL, the model will be validated with the same data used for fitting.
#' @param background data.frame. Database with response column only with 0 and predictors variables. All
#' column names must be consistent with data
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1), useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. A vector must be provided for this argument. The following threshold criteria are available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which the Jaccard index is the highest.
#'   \item max_sorensen: The threshold at which the Sorensen index is the highest.
#'   \item max_fpb: The threshold at which FPB (F-measure on presence-background data) is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to sensitivity value. If a sensitivity value is not specified, the default used is 0.9.
#'   }
#' If more than one threshold type is used they must be concatenated, e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function will use all threshold criteria if none is specified.
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "graf" class object. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: Performance metrics (see \code{\link{sdm_eval}}).
#' Threshold dependent metrics are calculated based on the threshold criteria specified in the argument.
#' \item performance_part: Performance metric for each replica and partition (see \code{\link{sdm_eval}}).
#' \item data_ens: Predicted suitability for each test partition. This database is used in \code{\link{fit_ensemble}}
#' }
#'
#' @seealso \code{\link{fit_gam}}, \code{\link{fit_glm}}, \code{\link{fit_gbm}},
#' \code{\link{fit_max}}, \code{\link{fit_net}}, \code{\link{fit_raf}}, and \code{\link{fit_svm}}.
#'
#' @export
#'
#' @importFrom dplyr %>% select all_of starts_with filter pull bind_rows mutate tibble group_by summarise across relocate left_join
#' @importFrom stats sd dist dnorm pnorm qnorm
#'
#' @examples
#' \dontrun{
#' data("abies")
#'
#' # Using k-fold partition method
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 3)
#' )
#' abies2
#'
#' bg <- abies2
#' bg$pr_ab <- 0
#'
#'
#' gaup_t1 <- fit_gau(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   background = bg,
#'   thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
#' )
#'
#' gaup_t1$model
#' gaup_t1$predictors
#' gaup_t1$performance
#' gaup_t1$performance_part
#' gaup_t1$data_ens
#'
#' # Using bootstrap partition method only with presence-absence
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "boot", replicates = 5, proportion = 0.7)
#' )
#' abies2
#'
#' gaup_t2 <- fit_gau(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("ppt_jja", "pH", "awc"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = c(type = c("lpt", "max_sens_spec", "sensitivity"), sens = "0.8")
#' )
#' gaup_t2
#' }
fit_gau <- function(data,
                    response,
                    predictors,
                    predictors_f = NULL,
                    background = NULL,
                    partition = NULL,
                    thr = NULL) {
  . <- model <- TPR <- IMAE <- rnames <- thr_value <- n_presences <- n_absences <- NULL
  variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))

  data <- data.frame(data)
  if (!is.null(background)) background <- data.frame(background)

  if (is.null(predictors_f)) {
    data <- data %>%
      dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), if (!is.null(partition)) dplyr::starts_with(partition))
    if (!is.null(background)) {
      background <- background %>%
        dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), if (!is.null(partition)) dplyr::starts_with(partition))
    }
  } else {
    data <- data %>%
      dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), dplyr::all_of(predictors_f), if (!is.null(partition)) dplyr::starts_with(partition))
    data <- data.frame(data)
    for (i in predictors_f) {
      data[, i] <- as.factor(data[, i])
    }
    if (!is.null(background)) {
      background <- background %>%
        dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), dplyr::all_of(predictors_f), if (!is.null(partition)) dplyr::starts_with(partition))
      for (i in predictors_f) {
        background[, i] <- as.factor(background[, i])
      }
    }
  }

  if (!is.null(background)) {
    if (!all(table(c(names(background), names(data))) == 2)) {
      stop("Column names of database used in 'data' and 'background' arguments do not match")
    }
  }

  # Remove NAs
  complete_vec <- stats::complete.cases(data[, c(response, unlist(variables))])
  if (sum(!complete_vec) > 0) {
    message(sum(!complete_vec), " rows were excluded from database because NAs were found")
    data <- data %>% dplyr::filter(complete_vec)
  }
  rm(complete_vec)
  if (!is.null(background)) {
    complete_vec <- stats::complete.cases(background[, c(response, unlist(variables))])
    if (sum(!complete_vec) > 0) {
      message(sum(!complete_vec), " rows were excluded from database because NAs were found")
      background <- background %>% dplyr::filter(complete_vec)
    }
    rm(complete_vec)
  }


  # Compare pr_ab and background column names
  if (!is.null(partition)) {
    p_names <- names(data %>% dplyr::select(dplyr::starts_with(partition)))
    for (i in p_names) {
      if (!is.null(background)) {
        Npart_p <- data %>%
          dplyr::filter(!!as.symbol(response) == 1) %>%
          dplyr::pull({{ i }}) %>%
          unique() %>%
          sort()
        Npart_bg <- background %>%
          dplyr::filter(!!as.symbol(response) == 0) %>%
          dplyr::pull({{ i }}) %>%
          unique() %>%
          sort()
        if (!all(table(c(Npart_p, Npart_bg)) == 2)) {
          stop(
            paste(
              "Partition groups between presences and background do not match:\n",
              paste("Part. group presences:", paste(Npart_p, collapse = " "), "\n"),
              paste("Part. group background:", paste(Npart_bg, collapse = " "), "\n")
            )
          )
        }
      }
    }
    rm(i)
  }


  # Fit models
  if (is.null(partition)) {
    suppressWarnings(mod <-
      graf(
        y = data[, response],
        x = data[, c(predictors, predictors_f)],
        opt.l = FALSE,
        method = "Laplace"
      ))

    result <- list(
      model = mod
    )
    return(result)
  } else {
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

      # In the follow code function will substitutes absences by background points
      # only in train database in order to fit gau with presences and background
      # and validate models with presences and absences
      if (!is.null(background)) {
        background2 <- pre_tr_te(background, p_names, h)
        train <- lapply(train, function(x) x[x[, response] == 1, ])
        train <- mapply(dplyr::bind_rows, train, background2$train, SIMPLIFY = FALSE)
        bgt_test <- background2$test
        rm(background2)
      }

      for (i in 1:np2) {
        tryCatch({
          message("Partition number: ", i, "/", np2)
          set.seed(1)
          mod[[i]] <-
            graf(
              y = train[[i]][, response],
              x = train[[i]][, c(predictors, predictors_f)],
              opt.l = FALSE,
              method = "Laplace"
            )

          # Predict for presences absences data
          ## Eliminate factor levels not used in fitting
          if (!is.null(predictors_f)) {
            for (fi in 1:length(predictors_f)) {
              lev <- train[[i]][, predictors_f[fi]] %>%
                unique() %>%
                as.character()
              lev_filt <- test[[i]][, predictors_f[fi]] %in% lev
              test[[i]] <- test[[i]][lev_filt, ]
              if (!is.null(background)) {
                lev_filt <- bgt_test[[i]][, predictors_f[fi]] %in% lev
                bgt_test[[i]] <- bgt_test[[i]][lev_filt, ]
              }
            }
          }

          # Predict for presences absences data
          pred_test <- data.frame(
            pr_ab = test[[i]][, response],
            pred =
              suppressWarnings(
                predict.graf(
                  mod[[i]],
                  newdata = test[[i]][c(predictors, predictors_f)],
                  type = "response",
                  CI = NULL
                )[, 1]
              )
          )

          pred_test_ens[[h]][[i]] <- pred_test %>%
            dplyr::mutate(rnames = rownames(test[[i]]))

          # Predict for background data
          if (!is.null(background)) {
            bgt <-
              data.frame(
                pr_ab = bgt_test[[i]][, response],
                pred = suppressWarnings(
                  predict.graf(
                    mod[[i]],
                    newdata = bgt_test[[i]][c(predictors, predictors_f)],
                    type = "response",
                    CI = NULL
                  )[, 1]
                )
              )
          }

          # Validation of model
          if (is.null(background)) {
            eval <-
              sdm_eval(
                p = pred_test$pred[pred_test$pr_ab == 1],
                a = pred_test$pred[pred_test$pr_ab == 0],
                thr = thr
              )
          } else {
            eval <-
              sdm_eval(
                p = pred_test$pred[pred_test$pr_ab == 1],
                a = pred_test$pred[pred_test$pr_ab == 0],
                thr = thr,
                bg = bgt$pred
              )
          }
          eval_partial[[i]] <- dplyr::tibble(model = "gau", eval)
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
        dplyr::bind_rows(x, .id = "part")
      }) %>%
      dplyr::bind_rows(., .id = "replicates") %>%
      dplyr::tibble() %>%
      dplyr::relocate(rnames)

    # Fit final models with best settings
    set.seed(1)
    suppressMessages(mod <-
      graf(
        y = data[, response],
        x = data[, c(predictors, predictors_f)],
        opt.l = FALSE,
        method = "Laplace"
      ))

    suppressWarnings(pred_test <- data.frame(
      pr_ab = data.frame(data)[, response],
      pred = suppressMessages(
        predict.graf(
          mod,
          newdata = data[, c(predictors, predictors_f)],
          type = "response",
          CI = NULL
        )[, 1]
      )
    ))

    if (is.null(background)) {
      threshold <- sdm_eval(
        p = pred_test$pred[pred_test$pr_ab == 1],
        a = pred_test$pred[pred_test$pr_ab == 0],
        thr = thr
      )
    } else {
      background <- suppressWarnings(predict.graf(
        mod,
        newdata = background[c(predictors, predictors_f)],
        type = "response",
        CI = NULL
      )[, 1])
      threshold <- sdm_eval(
        p = pred_test$pred[pred_test$pr_ab == 1],
        a = pred_test$pred[pred_test$pr_ab == 0],
        thr = thr,
        bg = background
      )
    }

    result <- list(
      model = mod,
      predictors = variables,
      performance = dplyr::left_join(eval_final, threshold[1:4], by = "threshold") %>%
        dplyr::relocate(model, threshold, thr_value, n_presences, n_absences),
      performance_part = eval_partial,
      data_ens = pred_test_ens
    )
    return(result)
  }
}
