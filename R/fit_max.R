#' Fit and validate Maximum Entropy models
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous variables).
#' Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative
#' predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor
#' variables. See maxnet.formula function from maxnet package.
#' Note that the variables used here must be consistent with those used in
#' response, predictors, and predictors_f arguments. Default NULL.
#' @param partition character. Column name with training and validation partition groups.
#' @param background data.frame. Database including only those rows with 0 values in the response column and the predictors variables. All
#' column names must be consistent with data. Default NULL
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1), needed for threshold-dependent performance metrics. More than one threshold type can be used. It is necessary to provide a vector for this argument. The following threshold criteria are available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which the Jaccard index is the highest.
#'   \item max_sorensen: The threshold at which the Sorensen index is highest.
#'   \item max_fpb: The threshold at which FPB (F-measure on presence-background data) is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to sensitivity value. If  a sensitivity values is not specified the default used is 0.9.
#'   }
#' If more than one threshold type is used they must be concatenated, e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function will use all thresholds if no threshold is specified.
#'
#' @param clamp logical. If TRUE, predictors and features are restricted to the range seen during model training.
#' @param classes character. A single feature of any combinations of them. Features are symbolized by letters: l (linear), q (quadratic), h (hinge), p (product), and t (threshold). Usage classes = "lpq". Default "default" (see details).
#' @param pred_type character. Type of response required available "link", "exponential", "cloglog" and "logistic". Default "cloglog"
#' @param regmult numeric. A constant to adjust regularization. Default 1.
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "maxnet" class object from maxnet package. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: Performance metrics (see \code{\link{sdm_eval}}).
#' Threshold dependent metrics are calculated based on the threshold specified in thr argument.
#' \item data_ens: Predicted suitability for each test partition based on the best model. This database is used in \code{\link{fit_ensemble}}
#' }
#' @details
#' When the argument “classes” is set as default MaxEnt will use different features combination
#' depending of the number of presences (np) with the follow rule:
#'  if np < 10 classes = "l",
#'  if np between 10  and 15 classes = "lq",
#'  if np between 15 and 80 classes = "lqh",
#'  and if np >= 80 classes = "lqph"
#'
#'
#' @seealso \code{\link{fit_gam}}, \code{\link{fit_gau}}, \code{\link{fit_gbm}},
#' \code{\link{fit_glm}}, \code{\link{fit_net}}, \code{\link{fit_raf}}, and \code{\link{fit_svm}}.
#'
#' @export
#'
#' @importFrom dplyr %>% select starts_with filter pull bind_rows group_by summarise across everything all_of
#' @importFrom maxnet maxnet maxnet.formula
#' @importFrom stats sd
#'
#' @examples
#' \dontrun{
#' data("abies")
#' data("backg")
#' abies # environmental conditions of presence-absence data
#' backg # environmental conditions of background points
#'
#' # Using k-fold partition method
#' # Note that the partition method, number of folds or replications must
#' # be the same for presence-absence and background points datasets
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 5)
#' )
#' abies2
#'
#' backg2 <- part_random(
#'   data = backg,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 5)
#' )
#' backg2
#'
#' max_t1 <- fit_max(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   background = backg2,
#'   thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
#'   clamp = TRUE,
#'   classes = "default",
#'   pred_type = "cloglog",
#'   regmult = 1
#' )
#' length(max_t1)
#' max_t1$model
#' max_t1$predictors
#' max_t1$performance
#' max_t1$data_ens
#' }
fit_max <- function(data,
                    response,
                    predictors,
                    predictors_f = NULL,
                    fit_formula = NULL,
                    partition,
                    background = NULL,
                    thr = NULL,
                    clamp = TRUE,
                    classes = "default",
                    pred_type = "cloglog",
                    regmult = 1) {
  . <- model <- TPR <- IMAE <- rnames <- thr_value <- n_presences <- n_absences <- NULL
  variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))

  # Test response variable
  r_test <- (data %>% dplyr::pull(response) %>% unique() %>% na.omit())
  if ((!all(r_test %in% c(0, 1)))) {
    stop("values of response variable do not match with 0 and 1")
  }

  data <- data.frame(data)
  if (!is.null(background)) background <- data.frame(background)

  if (is.null(predictors_f)) {
    data <- data %>%
      dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), dplyr::starts_with(partition))
    if (!is.null(background)) {
      background <- background %>%
        dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), dplyr::starts_with(partition))
    }
  } else {
    data <- data %>%
      dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), dplyr::all_of(predictors_f), dplyr::starts_with(partition))
    data <- data.frame(data)
    for (i in predictors_f) {
      data[, i] <- as.factor(data[, i])
    }
    if (!is.null(background)) {
      background <- background %>%
        dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), dplyr::all_of(predictors_f), dplyr::starts_with(partition))
      for (i in predictors_f) {
        background[, i] <- as.factor(background[, i])
      }
    }
  }

  if (!is.null(background)) {
    if (!all(table(c(names(background), names(data))) == 2)) {
      print(table(c(names(background), names(data))))
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

  # Formula
  if (is.null(fit_formula)) {
    formula1 <- maxnet::maxnet.formula(data[response],
      data[predictors],
      classes = classes
    )
  } else {
    formula1 <- fit_formula
  }
  message(
    "Formula used for model fitting:\n",
    Reduce(paste, deparse(formula1)) %>% gsub(paste("  ", "   ", collapse = "|"), " ", .),
    "\n"
  )


  # Compare pr_ab and background column names
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

  # New predictor vectors
  if (!is.null(predictors_f)) {
    predictors <- c(predictors, predictors_f)
  }

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

    # In the follow code function will substitutes absences by background points
    # only in train database in order to fit maxent with presences and background
    # and validate models with presences and absences
    if (!is.null(background)) {
      background2 <- pre_tr_te(background, p_names, h)
      train <- lapply(train, function(x) x[x[, response] == 1, ])
      train <- mapply(dplyr::bind_rows, train, background2$train, SIMPLIFY = FALSE)
      bgt_test <- background2$test
      rm(background2)
    }

    eval_partial <- as.list(rep(NA, np2))
    pred_test <- list()
    mod <- list()


    for (i in 1:np2) {
      message("Partition number: ", i, "/", np2)
      tryCatch({
        sampleback = TRUE
        try(mod[[i]] <-
              suppressMessages(
                maxnet::maxnet(
                  p = train[[i]][, response],
                  data = train[[i]][predictors],
                  f = formula1,
                  regmult = regmult,
                  addsamplestobackground = sampleback
                )
              ))
        if (length(mod) < i) {
          message("Refit with addsamplestobackground = FALSE")
          sampleback = FALSE
          try(mod[[i]] <-
                suppressMessages(
                  maxnet::maxnet(
                    p = train[[i]][, response],
                    data = train[[i]][predictors],
                    f = formula1,
                    regmult = regmult,
                    addsamplestobackground = sampleback
                  )
                ))
        }

        # Predict for presences absences data
        ## Eliminate factor levels not used in fitting
        # if (!is.null(predictors_f)) {
        #   for (fi in 1:length(predictors_f)) {
        #     lev <- as.character(unique(mod[[i]]$x[, predictors_f[fi]]))
        #     lev_filt <- test[[i]][, predictors_f[fi]] %in% lev
        #     test[[i]] <- test[[i]][lev_filt, ]
        #     if (!is.null(background)) {
        #       lev_filt <- bgt_test[[i]][, predictors_f[fi]] %in% lev
        #       bgt_test[[i]] <- bgt_test[[i]][lev_filt, ]
        #     }
        #   }
        # }

        # Predict for presences absences data
        pred_test <- data.frame(
          pr_ab = test[[i]][, response],
          pred = predict_maxnet(
            mod[[i]],
            newdata = test[[i]],
            clamp = clamp,
            type = pred_type,
            addsamplestobackground = sampleback
          )
        )

        pred_test_ens[[h]][[i]] <- pred_test %>%
          dplyr::mutate(rnames = rownames(test[[i]]))

        # Predict for background data
        if (!is.null(background)) {
          bgt <-
            data.frame(
              pr_ab = bgt_test[[i]][, response],
              pred =
                predict_maxnet(
                  mod[[i]],
                  newdata = bgt_test[[i]][c(predictors, predictors_f)],
                  clamp = clamp,
                  type = pred_type,
                  addsamplestobackground = sampleback
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

        eval_partial[[i]] <- dplyr::tibble(model = "max", eval)

        names(eval_partial) <- i
      })
    }

    # Create final database with parameter performance
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
  for (e in 1:length(pred_test_ens)) {
    fitl <- sapply(pred_test_ens[[e]], function(x) !is.null(nrow(x)))
    pred_test_ens[[e]] <- pred_test_ens[[e]][fitl]
  }

  pred_test_ens <-
    lapply(pred_test_ens, function(x) {
      bind_rows(x, .id = "part")
    }) %>%
    bind_rows(., .id = "replicates") %>%
    dplyr::tibble() %>%
    dplyr::relocate(rnames)

  # Fit final models
  suppressMessages(mod <-
    maxnet::maxnet(
      p = data[, response],
      data = data[predictors],
      f = formula1,
      regmult = regmult
    ))

  pred_test <- data.frame(
    "pr_ab" = data[response],
    "pred" = predict_maxnet(
      mod,
      newdata = data,
      clamp = clamp,
      type = pred_type
    )
  )

  if (is.null(background)) {
    threshold <- sdm_eval(
      p = pred_test$pred[pred_test$pr_ab == 1],
      a = pred_test$pred[pred_test$pr_ab == 0],
      thr = thr
    )
  } else {
    background <- predict_maxnet(
      mod,
      newdata = background[c(predictors, predictors_f)],
      clamp = clamp,
      type = pred_type
    )
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
    performance = dplyr::left_join(eval_final, threshold[1:4], by = "threshold") %>% dplyr::relocate(model, threshold, thr_value, n_presences, n_absences),
    data_ens = pred_test_ens
  )
  return(result)
}
