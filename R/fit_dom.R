#' Fit and validate Domain models
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous variables).
#' Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative
#' predictor variables (i.e. ordinal or nominal variables; factors). Usage predictors_f = c("landform")
#' @param partition character. Column name with training and validation partition groups.
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). This is useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The following threshold criteria are available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which the Jaccard index is the highest.
#'   \item max_sorensen: The threshold at which the Sorensen index is highest.
#'   \item max_fpb: The threshold at which FPB (F-measure on presence-background data) is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to sensitivity value. If a sensitivity value is not specified, the default used is 0.9.
#'   }
#' If more than one threshold type is used they must be concatenated, e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function will use all threshold types if none is specified.
#'
#' @param fit_formula formula. A formula object with response and predictor
#' variables (e.g. formula(pr_ab ~ aet + ppt_jja + pH + awc + depth + landform)).
#' Note that the variables used here must be consistent with those used in
#' response, predictors, and predictors_f arguments
#'
#' @details
#' This function fits and validates Domain models. The Domain model is a simple model that uses the Gower distance to
#' calculate environmental similarity between the presence data and test data.
#' Gower range of values area based on presences data. Gower distance are transformed to max(0, 1 - Gower).
#' This involves subtracting the distance from 1 and then ensuring the result is not negative (clamping it at zero).
#' Gower distance is calculated with \code{\link{map_env_dist}} funcion
#'
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A tibble with presences. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: Performance metric (see \code{\link{sdm_eval}}).
#' Threshold dependent metrics are calculated based on the threshold specified in the argument.
#' \item data_ens: Predicted suitability for each test partition. This database is used in \code{\link{fit_ensemble}}
#' }
#'
#' @details
#' This function fit and validate Domain models. The Domain model is a simple model that uses the Gower distance to calculate
#' the similarity between the presences training and presence-absences test data.
#'
#' @seealso \code{\link{fit_gam}}, \code{\link{fit_gau}}, \code{\link{fit_gbm}}, \code{\link{fit_glm}},
#' \code{\link{fit_max}}, \code{\link{fit_net}}, \code{\link{fit_raf}}, and \code{\link{fit_svm}}.
#'
#' @export
#'
#' @importFrom dplyr bind_rows select all_of starts_with filter reframe across everything mutate tibble group_by summarise relocate left_join
#' @importFrom stats complete.cases formula na.exclude sd
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' require(terra)
#'
#' data("spp")
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' # Extract data
#' some_sp <- spp %>%
#'   filter(species == "sp2")
#'
#' some_sp <-
#'   sdm_extract(
#'     data = some_sp,
#'     x = "x",
#'     y = "y",
#'     env_layer = somevar
#'   )
#'
#' # Partition
#' some_sp <- part_random(
#'   data = some_sp,
#'   pr_ab = "pr_ab",
#'   method = c(method = "rep_kfold", folds = 3, replicates = 5)
#' )
#'
#'
#' ##%######################################################%##
#' #                                                          #
#' ####                Fit a Domain model                  ####
#' #                                                          #
#' ##%######################################################%##
#' # Fit some models
#' mdom <- fit_dom(
#'   data = some_sp,
#'   response = "pr_ab",
#'   predictors = c("CFP_1", "CFP_2", "CFP_3", "CFP_4"),
#'   predictors_f = NULL,
#'   fit_formula = NULL,
#'   partition = ".part",
#'   thr = c("max_sens_spec")
#' )
#'
#' mdom
#'
#' # Predict model
#' ind_p <- sdm_predict(
#'   models = mdom,
#'   pred = somevar,
#'   thr = "max_sens_spec",
#'   con_thr = TRUE,
#'   predict_area = NULL
#' )
#' plot(ind_p$dom)
#'
#' ##%######################################################%##
#' #                                                          #
#' ####             Explore Domain suitabiltiy             ####
#' ####             in the environmental space             ####
#' #                                                          #
#' ##%######################################################%##
#'
#' p_extra(
#'   training_data = some_sp %>% dplyr::filter(pr_ab == 1), #select only presences
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   extra_suit_data = ind_p$dom$dom,
#'   projection_data = somevar,
#'   geo_space = FALSE,
#'   prop_points = 0.3,
#'   alpha_p = 0.8,
#'   color_p = "black",
#'   color_gradient = c("#000033", "#1400FF", "#C729D6", "#FF9C63", "#FFFF60")
#' )
#'
#' p_extra(
#'   training_data = some_sp %>% dplyr::filter(pr_ab == 1), #select only presences
#'   x = "x",
#'   y = "y",
#'   pr_ab = "pr_ab",
#'   predictors = c("CFP_1", "CFP_2"), # Just the first two predictors
#'   extra_suit_data = ind_p$dom$dom > 0.96, # a binary map
#'   projection_data = somevar,
#'   geo_space = TRUE,
#'   prop_points = 0.4,
#'   alpha_p = 0.8,
#'   color_p = "black",
#'   color_gradient = c("#1400FF", "#C729D6")
#' )
#'
#' }
fit_dom <- function(data,
                    response,
                    predictors,
                    predictors_f = NULL,
                    partition,
                    thr = NULL,
                    fit_formula = NULL) {
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
    formula1 <-
      paste(c(
        paste(predictors, collapse = " + ", sep = ""),
        predictors_f
      ), collapse = " + ")
    formula1 <- stats::formula(paste(
      response, "~", formula1
    ))
  } else {
    formula1 <- fit_formula
  }
  message(
    "Formula used for model fitting:\n",
    Reduce(paste, deparse(formula1)) %>% gsub(paste("  ", "   ", collapse = "|"), " ", .),
    "\n"
  )


  # Calculate range for each column
  # range_var <- data[, predictors] %>% dplyr::reframe(dplyr::across(dplyr::everything(), range))

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
    train <- lapply(train, function(x) x[x[, response] == 1, ]) # this algorithm only works with presence data
    test <- out$test
    np2 <- out$np2
    rm(out)

    eval_partial <- as.list(rep(NA, np2))
    pred_test <- list()
    mod <- list()

    for (i in 1:np2) {
      message("Partition number: ", i, "/", np2)
      tryCatch(
        {
          mod[[i]] <- min_gower_rcpp(
            train[[i]][, c(predictors, predictors_f)],
            test[[i]][, c(predictors, predictors_f)]
          )

          pred_test <- data.frame(
            pr_ab = test[[i]][, response],
            pred = mod[[i]]
          )

          pred_test_ens[[h]][[i]] <- pred_test %>%
            dplyr::mutate(rnames = rownames(test[[i]]))

          # Validation of model
          eval <-
            sdm_eval(
              p = pred_test$pred[pred_test$pr_ab == 1],
              a = pred_test$pred[pred_test$pr_ab == 0],
              thr = thr
            )
          eval_partial[[i]] <- dplyr::tibble(model = "dom", eval)
        },
        error = function(cond) {
          message("Sorry, but it was not possible to fit the model with this data")
        }
      )
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
  # mod <- min_gower_rcpp(
  #   data[data[, response] == 1, c(predictors, predictors_f)],
  #   data[, c(predictors, predictors_f)]
  # )
  #
  # pred_test <- data.frame(
  #   pr_ab = data[, response],
  #   pred = mod
  # )
  #
  # threshold <- sdm_eval(
  #   p = pred_test$pred[pred_test$pr_ab == 1],
  #   a = pred_test$pred[pred_test$pr_ab == 0],
  #   thr = thr
  # )

  da_pres <- data[data[, response] == 1, c(predictors, predictors_f)]
  nnn <- ifelse(nrow(da_pres)>=100, 50, nrow(da_pres))

  threshold <- list()
  for(i in 1:nnn){
    if(nnn==50){
      set.seed(i)
      mod <- min_gower_rcpp(da_pres %>% dplyr::slice_sample(prop = 0.50),
                             data[, c(predictors, predictors_f)])
    } else {
      mod <- min_gower_rcpp(da_pres[-i,],
                             data[, c(predictors, predictors_f)])
    }

    pred_test <- data.frame(
      pr_ab = data[, response],
      pred = mod)

    threshold[[i]] <- sdm_eval(
      p = pred_test$pred[pred_test$pr_ab == 1],
      a = pred_test$pred[pred_test$pr_ab == 0],
      thr = thr)
  }

  threshold <- threshold %>%
    dplyr::bind_rows(threshold) %>%
    dplyr::group_by(threshold) %>%
    dplyr::filter(!is.na(n_presences)) %>%
    dplyr::summarise(thr_value=mean(thr_value))
  threshold <- dplyr::bind_cols(threshold, sdm_eval(
    p = pred_test$pred[pred_test$pr_ab == 1],
    a = pred_test$pred[pred_test$pr_ab == 0],
    thr = thr)[c("n_presences", "n_absences")]
  )

  # create a new object similar than data.frame
  result <- list(
    model = list(domain = data[data[, response] == 1, c(predictors, predictors_f)]),
    predictors = variables,
    performance = dplyr::left_join(eval_final, threshold[1:4], by = "threshold") %>%
      dplyr::relocate(model, threshold, thr_value, n_presences, n_absences),
    data_ens = pred_test_ens
  )
  return(result)
}
