#' Fit and validate Generalized Linear Models
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous or discrete variables).
#' Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative
#' predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param partition character. Column name with training and validation partition groups.
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1).
#' It is useful for threshold-dependent performance metrics.
#' It is possible to use more than one threshold type. It is necessary to provide a
#' vector for this argument. The next threshold area available:
#' \itemize{
#'   \item LPT: The highest threshold at which there is no omission. Usage thr=c(type='LPT').
#'   \item EQUAL_SENS_SPEC: Threshold at which the sum of the sensitivity and specificity is the highest.
#'   \item MAX_TSS: Threshold at which the sensitivity and specificity are equal.
#'   Usage thr=c(type='MAX_TSS').
#'   \item MAX_KAPPA: The threshold at which kappa is the highest ("max kappa"). Usage thr=c(type='MAX_KAPPA').
#'   \item MAX_JACCARD: The threshold at which Jaccard is the highest. Usage thr=c(type='MAX_JACCARD').
#'   \item MAX_SORENSEN: The threshold at which Sorensen is highest. Usage thr=c(type='MAX_SORENSEN').
#'   \item MAX_FPB: The threshold at which Fpb is highest. Usage thr=c(type='MAX_FPB').
#'   \item SENSITIVITY: A threshold value specified by user. Usage thr=c(type='SENSITIVITY', sens='0.6'). 'sens' refers to models will be binarized using this suitability value.
#'   }
#' @param fit_formula formula. A formula object with response and predictor
#' variables (e.g. forumla(pr_ab ~ aet + ppt_jja + pH + awc + depth + landform)).
#' Note that the variables used here must be consistent with those used in
#' response, predictors, and predictors_f arguments
#' @param poly interger >= 2. If used with values >= 2 model will use polinomius
#' for those continuous variables (i.e. used in predictors argument)
#' @param inter_order interger >= 0. The interaction order between explanatory variables.
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "glm" class object. This object can be used for predicting.
#' \item performance: Performance metric (see \code{\link{enm_eval}}).
#' Those threshold dependent metric are calculated based on the threshold specified in thr argument .
#' \item selected_threshold: Value of the threshold selected.
#' \item threshold_table: Value of all threshold.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("abies_db")
#' abies_db
#'
#' # Using KFOLD partition method
#' abies_db2 <- data_part(
#'   data = abies_db,
#'   p_a = "pr_ab",
#'   method = c(method = "KFOLD", folds = 10)
#' )
#' abies_db2
#'
#' glm_t1 <- fit_glm(
#'   data = abies_db2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = "MAX_TSS",
#'   poly = 0,
#'   inter_order = 0
#' )
#' glm_t1$model
#' glm_t1$performance
#' glm_t1$selected_threshold
#' glm_t1$threshold_table
#'
#'
#' glm_t2 <- fit_glm(
#'   data = abies_db2,
#'   response = "pr_ab",
#'   predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = "MAX_TSS",
#'   poly = 2,
#'   inter_order = 1
#' )
#'
#' # Using REP_KFOLD partition method
#' abies_db2 <- data_part(
#'   data = abies_db,
#'   p_a = "pr_ab",
#'   method = c(method = "REP_KFOLD", folds = 10, replicates = 10)
#' )
#' abies_db2
#'
#' glm_t3 <- fit_glm(
#'   data = abies_db2,
#'   response = "pr_ab",
#'   predictors = c("ppt_jja", "pH", "awc"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   thr = "MAX_TSS",
#'   poly = 3,
#'   inter_order = 2
#' )
#' }
#'
fit_glm <- function(data,
                    response,
                    predictors,
                    predictors_f = NULL,
                    partition,
                    thr = NULL,
                    fit_formula = NULL,
                    poly = 0,
                    inter_order = 0,
                    ...) {
  data <- data.frame(data)

  if (is.null(predictors_f)) {
    data <- data %>%
      dplyr::select(all_of(response), all_of(predictors), dplyr::starts_with(partition))
    data <- data.frame(data)
  } else {
    data <- data %>%
      dplyr::select(all_of(response), all_of(predictors), all_of(predictors_f), dplyr::starts_with(partition))
    data <- data.frame(data)
    for (i in predictors_f) {
      data[, i] <- as.factor(data[, i])
    }
  }

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
    message(
      "Formula used for model fitting:\n",
      Reduce(paste, deparse(formula1)) %>% gsub(paste("  ", "   ", collapse = "|"), " ", .)
    )
  } else {
    formula1 <- fit_formula
  }


  # Fit models
  np <- ncol(data %>% dplyr::select(dplyr::starts_with(partition)))
  p_names <- names(data %>% dplyr::select(dplyr::starts_with(partition)))
  eval_partial_list <- list()
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
        stats::glm(formula1,
          data = train[[i]],
          family = binomial
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
          stats::predict.glm(
            mod[[i]],
            newdata = test[[i]],
            type = "response",
            se.fit = FALSE
          )
        )
      ))

      # Validation of model
      eval <-
        enm_eval(
          p = pred_test$pred[pred_test$pr_ab == 1],
          a = pred_test$pred[pred_test$pr_ab == 0],
          thr = thr
        )
      eval_partial[[i]] <- eval$selected_threshold
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
    dplyr::select(-c(replica:n_absences)) %>%
    dplyr::summarise(dplyr::across(
      dplyr::everything(),
      list(mean = mean, sd = sd)
    ), .groups = "drop")

  # Fit final models with best settings
  suppressWarnings(mod <-
    stats::glm(formula1,
      data = data,
      family = binomial
    ))

  pred_test <- data.frame(
    pr_ab = data[, response],
    pred = suppressMessages(stats::predict.glm(
      mod,
      newdata = data,
      type = "response"
    ))
  )

  threshold <- enm_eval(
    p = pred_test$pred[pred_test$pr_ab == 1],
    a = pred_test$pred[pred_test$pr_ab == 0],
    thr = thr
  )

  result <- list(
    model = mod,
    performance = eval_final,
    selected_threshold = threshold[[1]] %>% dplyr::select(threshold:TNR),
    threshold_table = threshold[[2]] %>% dplyr::select(threshold:TNR)
  )
  return(result)
}
