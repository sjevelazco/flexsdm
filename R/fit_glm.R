#' Fit and validate Generalized Linear Models
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous variables).
#' Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative
#' predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param select_pred logical. Perform predictor selection.
#' If TRUE predictors will be selected based on backward step wise approach. Default FALSE.
#' @param partition character. Column name with training and validation partition groups.
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1), needed for threshold-dependent performance metrics. More than one threshold type can be used. It is necessary to provide a vector for this argument. The following threshold criteria are available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which the Jaccard index is the highest.
#'   \item max_sorensen: The threshold at which the Sorensen index is highest.
#'   \item max_fpb: The threshold at which FPB (F-measure on presence-background data) is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to sensitivity value. If it is not specified a sensitivity values, function will use by default 0.9
#'   }
#' If more than one threshold type is used they must be concatenated, e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function will use all thresholds if no threshold is specified.
#'
#' @param fit_formula formula. A formula object with response and predictor
#' variables (e.g. formula(pr_ab ~ aet + ppt_jja + pH + awc + depth + landform)).
#' Note that the variables used here must be consistent with those used in
#' response, predictors, and predictors_f arguments
#' @param poly integer >= 2. If used with values >= 2 model will use polynomials
#' for those continuous variables (i.e. used in predictors argument). Default is 0.
#' @param inter_order integer >= 0. The interaction order between explanatory variables. Default is 0.
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "glm" class object from stats package. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: Performance metrics (see \code{\link{sdm_eval}}).
#' Threshold dependent metric are calculated based on the threshold specified in thr argument .
#' \item data_ens: Predicted suitability for each test partition. This database is used in \code{\link{fit_ensemble}}
#' }
#'
#' @seealso \code{\link{fit_gam}}, \code{\link{fit_gau}}, \code{\link{fit_gbm}},
#' \code{\link{fit_max}}, \code{\link{fit_net}}, \code{\link{fit_raf}}, and \code{\link{fit_svm}}.
#'
#' @export
#'
#' @importFrom dplyr bind_rows select all_of starts_with filter mutate tibble group_by summarise across relocate left_join
#' @importFrom stats complete.cases formula glm step na.exclude predict.glm sd
#'
#' @examples
#' \dontrun{
#' data("abies")
#' abies
#'
#' # Using k-fold partition method
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 5)
#' )
#' abies2
#'
#' glm_t1 <- fit_glm(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   select_pred = FALSE,
#'   partition = ".part",
#'   thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
#'   poly = 0,
#'   inter_order = 0
#' )
#' glm_t1$model
#' glm_t1$predictors
#' glm_t1$performance
#' glm_t1$data_ens
#'
#' # Using second order polynomial terms and first-order interaction terms
#' glm_t2 <- fit_glm(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   select_pred = FALSE,
#'   partition = ".part",
#'   thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
#'   poly = 2,
#'   inter_order = 1
#' )
#'
#' # Using repeated k-fold partition method
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "rep_kfold", folds = 3, replicates = 5)
#' )
#' abies2
#'
#' # Using third order polynomial terms and second-order interaction terms
#' glm_t3 <- fit_glm(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("ppt_jja", "pH", "awc"),
#'   predictors_f = c("landform"),
#'   select_pred = FALSE,
#'   partition = ".part",
#'   thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen"),
#'   poly = 3,
#'   inter_order = 2
#' )
#' }
fit_glm <- function(data,
                    response,
                    predictors,
                    predictors_f = NULL,
                    select_pred = FALSE,
                    partition,
                    thr = NULL,
                    fit_formula = NULL,
                    poly = 2,
                    inter_order = 0) {
  . <- model <- TPR <- IMAE <- rnames <- thr_value <- n_presences <- n_absences <- NULL
  variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))

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
  complete_vec <- stats::complete.cases(data[, c(response, unlist(variables))])
  if (sum(!complete_vec) > 0) {
    message(sum(!complete_vec), " rows were excluded from database because NAs were found")
    data <- data %>% dplyr::filter(complete_vec)
  }
  rm(complete_vec)

  # Formula
  if (is.null(fit_formula)) {
    if (poly >= 2) {
      forpoly <- lapply(2:poly, function(x) {
        paste("I(", predictors, "^", x, ")",
          sep = "", collapse = " + "
        )
      }) %>% paste(collapse = " + ")
      formula1 <- paste(c(predictors, predictors_f), collapse = " + ") %>%
        paste(., forpoly, sep = " + ")
    } else {
      formula1 <-
        paste(c(predictors, predictors_f), collapse = " + ")
    }

    if (inter_order > 0) {
      forinter <- c(predictors, predictors_f)
      if (inter_order > length(forinter)) {
        stop("value of inter_order is higher than number of predicors ", "(", length(forinter), ")")
      }
      forinter_l <- list()

      for (i in 1:inter_order) {
        forinter_l[[i]] <- do.call(
          "expand.grid",
          c(lapply(1:(i + 1), function(x) {
            forinter
          }), stringsAsFactors = FALSE)
        )
        forinter_l[[i]] <- apply(forinter_l[[i]], 1, function(x) {
          x <- unique(sort(x))
          if (length(x) > i) {
            paste(x, collapse = ":")
          }
        }) %>%
          unlist() %>%
          unique()
      }
      forinter <- sapply(forinter_l, paste, collapse = " + ")
      forinter <- do.call("paste", c(as.list(forinter), sep = " + "))
    }

    if (exists("forinter")) {
      formula1 <- paste(formula1, forinter, sep = " + ")
      formula1 <- stats::formula(paste(
        response, "~", formula1
      ))
    } else {
      formula1 <- stats::formula(paste(
        response, "~", formula1
      ))
    }
  } else {
    formula1 <- fit_formula
  }
  message(
    "Formula used for model fitting:\n",
    Reduce(paste, deparse(formula1)) %>% gsub(paste("  ", "   ", collapse = "|"), " ", .),
    "\n"
  )

  # Selection predictor
  if (select_pred) {
    message("Selecting predictors")

    var_selected <- suppressWarnings(stats::glm(formula1, data = data, family = "binomial"))
    suppressWarnings(var_selected <- stats::step(var_selected, direction = "backward", trace = 0))
    formula1 <- formula(var_selected)

    var_selected <- as.character(formula1)[3] %>%
      strsplit(., "[ + ]")
    var_selected <- var_selected[[1]]
    var_selected <- var_selected[var_selected != ""]
    if (poly < 2) {
      poly <- NULL
    }

    var_selected <- var_selected %>%
      gsub(paste0("[", poly, ")]"), "", .) %>%
      gsub("[/^]", "", .) %>%
      gsub("[I(]", "", .)
    variables <- variables[variables %in% var_selected]
    predictors <- predictors[predictors %in% var_selected]
    predictors_f <- predictors_f[predictors_f %in% var_selected]
    if (length(predictors_f) == 0) {
      predictors_f <- NULL
    }


    message(
      "Formula used for model fitting:\n",
      Reduce(paste, deparse(formula1)) %>% gsub(paste("  ", "   ", collapse = "|"), " ", .),
      "\n"
    )
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

    eval_partial <- as.list(rep(NA, np2))
    pred_test <- list()
    mod <- list()

    for (i in 1:np2) {
      tryCatch({
        message("Partition number: ", i, "/", np2)
        suppressWarnings(mod[[i]] <-
          stats::glm(formula1,
            data = train[[i]],
            family = "binomial"
          ))


        # Predict for presences absences data
        ## Eliminate factor levels not used for fitting
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
            stats::predict.glm(
              mod[[i]],
              newdata = test[[i]],
              type = "response",
              se.fit = FALSE
            )
          )
        ))

        pred_test_ens[[h]][[i]] <- pred_test %>%
          dplyr::mutate(rnames = rownames(test[[i]]))

        # Validation of model
        eval <-
          sdm_eval(
            p = pred_test$pred[pred_test$pr_ab == 1],
            a = pred_test$pred[pred_test$pr_ab == 0],
            thr = thr
          )
        eval_partial[[i]] <- dplyr::tibble(model = "glm", eval)
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
      bind_rows(x, .id = "part")
    }) %>%
    bind_rows(., .id = "replicates") %>%
    dplyr::tibble() %>%
    dplyr::relocate(rnames)

  # Fit final models with best settings
  suppressWarnings(mod <-
    stats::glm(formula1,
      data = data,
      family = "binomial"
    ))

  pred_test <- data.frame(
    pr_ab = data.frame(data)[, response],
    pred = suppressMessages(stats::predict.glm(
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

  result <- list(
    model = mod,
    predictors = variables,
    performance = dplyr::left_join(eval_final, threshold[1:4], by = "threshold") %>% dplyr::relocate(model, threshold, thr_value, n_presences, n_absences),
    data_ens = pred_test_ens
  )
  return(result)
}
