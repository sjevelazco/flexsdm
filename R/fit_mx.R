#' Fit and validate Maximum Entropy models
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous or discrete variables).
#' Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative
#' predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor
#' variables see maxnet.formula function from maxnet package.
#' Note that the variables used here must be consistent with those used in
#' response, predictors, and predictors_f arguments. Default NULL.
#' @param partition character. Column name with training and validation partition groups.
#' @param background data.frame. Database with response column only with 0 and predictors variables. All
#' column names must be consistent with data
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The next threshold area available:
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
#' @param clamp logical. It is set with TRUE, predictors and features are restricted to the range seen during model training.
#' @param classes character. A single feature of any combinations of them. Features are symbolized by letters: l (linear), q (quadratic), h (hinge), p (product), and t (threshold). Usage classes = "lpq". Default "default" (see details).
#' @param pred_type character. Type of response required available "link", "exponential", "cloglog" and "logistic". Default "cloglog"
#' @param regmult numeric. A constant to adjust regularization. Default 1.
#' @param ...
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "MaxEnt" class object. This object can be used for predicting.
#' \item predictors: A character with quantitative (elements names with c) and qualitative (elements names with f) variables use for modeling.
#' \item performance: Performance metric (see \code{\link{sdm_eval}}).
#' Those threshold dependent metric are calculated based on the threshold specified in thr argument .
#' \item selected_thresholds: Value of the threshold selected.
#' \item all_thresholds: Value of all threshold.
#' }
#' @details EXPLAIN HERE DEFAUL SELECTION OF FEATURES BASED ON number of occurrences
#' if (np < 10) {
#'   classes <- "l"
#' } else if (np < 15){
#'   classes <- "lq"
#' } else if (np < 80) {
#'   classes <- "lqh"
#' }
#'
#' @export
#'
#' @importFrom dplyr select starts_with filter pull bind_rows group_by summarise across everything all_of
#' @importFrom maxnet maxnet maxnet.formula
#' @importFrom stats sd
#'
#' @examples
#' \dontrun{
#'
#' }
#'
fit_mx <- function(data,
                   response,
                   predictors,
                   predictors_f = NULL,
                   fit_formula = NULL,
                   partition,
                   background,
                   thr = NULL,
                   clamp = TRUE,
                   classes = "default",
                   pred_type = "cloglog",
                   regmult = 1,
                   ...) {
  variables <- c(c = predictors, f = predictors_f)

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
  data <- rm_na(x = data)
  if (!is.null(background)) {
    background <- rm_na(x = background)
  }

  # Formula
  if(is.null(fit_formula)){
    formula1 <- maxnet::maxnet.formula(data[response],
                                       data[predictors],
                                       classes = classes)
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

    eval_partial <- list()
    pred_test <- list()
    mod <- list()


    for (i in 1:np2) {
      message("Partition number: ", i, "/", np2)
      tryCatch({
      set.seed(1)
        mod[[i]] <-
          suppressMessages(
          maxnet::maxnet(
          p = train[[i]][, response],
          data = train[[i]][predictors],
          f = formula1,
          regmult = regmult
        )
        )


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
        pred =
          maxnet:::predict.maxnet(
            mod[[i]],
            newdata = test[[i]],
            clamp = clamp,
            type = pred_type
          )
      )

      pred_test_ens[[h]][[i]] <- pred_test %>%
        dplyr::mutate(rnames=rownames(.))

      # Predict for background data
      if (!is.null(background)) {
        bgt <-
          data.frame(
            pr_ab = bgt_test[[i]][, response],
            pred =
              maxnet:::predict.maxnet(
                mod[[i]],
                newdata = bgt_test[[i]][c(predictors, predictors_f)],
                clamp = clamp,
                type = pred_type
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

      if (is.null(thr)) {
        eval_partial[[i]] <- eval$all_thresholds
      } else {
        eval_partial[[i]] <- eval$selected_thresholds
      }

      names(eval_partial) <- i

      }
      ,
      error=function(cond) {
        message("It was not possible to fit this model")
      })

    }

    # Create final database with parameter performance
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
  for(e in 1:length(pred_test_ens)){
    fitl <- sapply(pred_test_ens[[e]], function(x) !is.null(nrow(x)))
    pred_test_ens[[e]] <- pred_test_ens[[e]][fitl]
  }

  pred_test_ens <-
    lapply(pred_test_ens, function(x)
      bind_rows(x, .id = 'part')) %>%
    bind_rows(., .id = 'replicates') %>% dplyr::tibble() %>%
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
    "pred" = maxnet:::predict.maxnet(
      mod,
      newdata = data,
      clamp = TRUE,
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
    background <- maxnet:::predict.maxnet(
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
    all_thresholds = threshold[[2]] %>% dplyr::select(threshold:values),
    data_ens = pred_test_ens
  )
  return(result)
}
