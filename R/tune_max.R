#' Fit and validate Maximum Entropy models with exploration of hyper-parameters that optimize performance
#'
#'
#' @param data data.frame. Database with response (0,1) and predictors values.
#' @param response character. Column name with species absence-presence data (0,1).
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous variables).
#' Usage predictors = c("aet", "cwd", "tmin")
#' @param predictors_f character. Vector with the column names of qualitative
#' predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param partition character. Column name with training and validation partition groups.
#' @param background data.frame. Database with response variable column only containing 0 values,
#' and predictors variables. All
#' column names must be consistent with data
#' @param grid data.frame. A data frame object with algorithm hyper-parameters values to be tested.
#' It is recommended to generate this data.frame with the grid() function. Hyper-parameters needed
#' for tuning are 'regmult' and 'classes' (any combination of following letters l -linear-, q
#'  -quadratic-, h -hinge-, p -product-, and t -threshold-).
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1)., needed for
#'  threshold-dependent performance metrics. More than one threshold type can be used. It is
#'  necessary to provide a vector for this argument. The following threshold types are available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the
#'   highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which the Jaccard index is the highest.
#'   \item max_sorensen: The threshold at which the Sorensen index is highest.
#'   \item max_fpb: The threshold at which # FPB (F-measure on presence-background data) is
#'   highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value.
#'   Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to
#'   sensitivity value. If  a sensitivity value is not specified, a default of 0.9 will be used.
#'   }
#' If more than one threshold type is used, concatenate them, e.g., thr=c('lpt', 'max_sens_spec',
#'  'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt',
#'   'max_sens_spec', 'sensitivity'). Function will use all thresholds if no threshold is specified.
#'
#' @param metric character. Performance metric used for selecting the best combination of hyper
#' -parameter values. One of the following metrics can be used: SORENSEN, JACCARD, FPB, TSS, KAPPA,
#' AUC, and BOYCE. TSS is used as default.
#' @param clamp logical. If TRUE, predictors and features are restricted to the range seen during
#'  model training.
#' @param pred_type character. Type of response required available "link", "exponential", "cloglog"
#' and "logistic". Default "cloglog"
#' @param n_cores numeric. Number of cores use for parallelization. Default 1
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "maxnet" "lognet" "glmnet" class object. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names)
#'  variables use for modeling.
#' \item performance: Hyper-parameters values and performance metrics (see \code{\link{sdm_eval}})
#' for the best hyper-parameters combination.
#' \item hyper_performance: Performance metrics (see \code{\link{sdm_eval}}) for each combination
#' of the hyper-parameters.
#' \item data_ens: Predicted suitability for each test partition based on the best model. This
#' database is used in \code{\link{fit_ensemble}}
#' }
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr bind_rows pull select all_of starts_with filter tibble left_join mutate group_by_at summarise across relocate
#' @importFrom foreach foreach
#' @importFrom maxnet maxnet maxnet.formula
#' @importFrom parallel makeCluster stopCluster
#' @importFrom stats complete.cases
#'
#' @seealso \code{\link{tune_gbm}}, \code{\link{tune_net}}, \code{\link{tune_raf}}, and \code{\link{tune_svm}}.
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("abies")
#' data("backg")
#' abies # environmental conditions of presence-absence data
#' backg # environmental conditions of background points
#'
#' # Using k-fold partition method
#' # Remember that the partition method, number of folds or replications must
#' # be the same for presence-absence and background points datasets
#' abies2 <- part_random(
#'   data = abies,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 3)
#' )
#' abies2
#'
#' set.seed(1)
#' backg <- dplyr::sample_n(backg, size = 2000, replace = FALSE)
#' backg2 <- part_random(
#'   data = backg,
#'   pr_ab = "pr_ab",
#'   method = c(method = "kfold", folds = 3)
#' )
#' backg
#'
#'
#' gridtest <-
#'   expand.grid(
#'     regmult = seq(0.1, 3, 0.5),
#'     classes = c("l", "lq", "lqh")
#'   )
#'
#' max_t1 <- tune_max(
#'   data = abies2,
#'   response = "pr_ab",
#'   predictors = c("aet", "pH", "awc", "depth"),
#'   predictors_f = c("landform"),
#'   partition = ".part",
#'   background = backg2,
#'   grid = gridtest,
#'   thr = "max_sens_spec",
#'   metric = "TSS",
#'   clamp = TRUE,
#'   pred_type = "cloglog",
#'   n_cores = 2 # activate two cores for speed up this process
#' )
#'
#' length(max_t1)
#' max_t1$model
#' max_t1$predictors
#' max_t1$performance
#' max_t1$data_ens
#' }
tune_max <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           background = NULL,
           partition,
           grid = NULL,
           thr = NULL,
           metric = "TSS",
           clamp = TRUE,
           pred_type = "cloglog",
           n_cores = 1) {
    . <- model <- TPR <- IMAE <- thr_value <- n_presences <- n_absences <- NULL
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
        stop("Column names of database used in 'data' and 'background' arguments do not match")
        print(table(c(names(background), names(data))))
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

    # Prepare grid when grid=default or NULL
    if (is.null(grid)) {
      grid <- expand.grid(
        regmult = seq(0.1, 3, 0.5),
        classes = c("l", "lq", "lqh", "lqhp", "lqhpt")
      )
      message("Hyper-parameter values were not provided, default values will be used")
      message("regmult = seq(0.1, 3, 0.5)")
      message('classes = c("l", "lq", "lqh", "lqhp", "lqhpt")')
    }

    # Test hyperparameter names
    hyperp <- names(grid)
    if (!all(c("regmult", "classes") %in% hyperp)) {
      stop("Database used in 'grid' argument has to contain this columns for tunning: 'regmult', 'classes'")
    }

    grid$tune <- 1:nrow(grid)


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

    # New predictor vectors
    if (!is.null(predictors_f)) {
      predictors <- c(predictors, predictors_f)
    }

    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

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

      eval_partial <- foreach::foreach(i = 1:np2, .export=c('sdm_eval', 'boyce'), .packages = c("dplyr")) %dopar%{
        # message("Partition number: ", i, "/", np2)
        mod <- as.list(rep(NA, nrow(grid)))
        names(mod) <- 1:nrow(grid)
        for (ii in 1:nrow(grid)) {
          set.seed(1)
          try(mod[[ii]] <-
            maxnet::maxnet(
              p = train[[i]][, response],
              data = train[[i]][predictors],
              f = maxnet::maxnet.formula(train[[i]][response],
                                         train[[i]][predictors],
                                         classes = grid$classes[ii]
              ),
              regmult = grid$regmult[ii]
            ))
        }

        filt <- sapply(mod, function(x) length(class(x)) == 3)
        filt <- filt&!is.na(sapply(mod, function(x) x$entropy))
        mod <- mod[filt]
        grid2 <- grid[filt, ]
        tnames <- apply(grid2, 1, function(x) paste(x, collapse = "_"))

        # Predict for presences absences data
        pred_test <-
          lapply(mod, function(x) {
            data.frame(
              pr_ab = test[[i]][, response],
              pred = predict_maxnet(
                x,
                newdata = test[[i]],
                clamp = clamp,
                type = pred_type
              )
            )
          })

        # Predict for background data
        if (!is.null(background)) {
          bgt <-
            lapply(mod, function(x) {
              data.frame(
                pr_ab = bgt_test[[i]][, response],
                pred = predict_maxnet(
                  x,
                  newdata = bgt_test[[i]],
                  clamp = clamp,
                  type = pred_type
                )
              )
            })
        }

        eval <- list()
        for (ii in 1:length(pred_test)) {
          if (is.null(background)) {
            eval[[ii]] <-
              sdm_eval(
                p = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 1],
                a = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 0],
                thr = thr
              ) %>% dplyr::tibble(model = "max", .)
          } else {
            eval[[ii]] <-
              sdm_eval(
                p = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 1],
                a = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 0],
                thr = thr,
                bg = bgt[[ii]]$pred
              ) %>% dplyr::tibble(model = "max", .)
          }
        }

        names(eval) <- tnames
        eval <- dplyr::bind_rows(eval, .id = "tnames")

        eval <-
          dplyr::tibble(dplyr::left_join(dplyr::mutate(grid2, tnames),
                                         eval,
                                         by = "tnames"
          )) %>%
          dplyr::select(-tnames)
        eval
      }
      parallel::stopCluster(cl)

      # Create final database with parameter performance 1
      names(eval_partial) <- 1:np2
      eval_partial <- eval_partial[sapply(eval_partial, function(x) !is.null(dim(x)))] %>%
        dplyr::bind_rows(., .id = "partition")
      eval_partial_list[[h]] <- eval_partial
    }

    # Create final database with parameter performance 2
    eval_partial <- eval_partial_list %>%
      dplyr::bind_rows(., .id = "replica")

    eval_final <- eval_partial %>%
      dplyr::group_by_at(c(hyperp, "model", "threshold")) %>%
      dplyr::summarise(dplyr::across(
        TPR:IMAE,
        list(mean = mean, sd = sd)
      ), .groups = "drop")

    filt <- eval_final %>% dplyr::pull(paste0(metric, "_mean"))
    filt <- which.max(filt)
    best_tune <- eval_final[filt, ]
    best_hyperp <- eval_final[filt, hyperp]


    # Fit final models with best settings
    # Get data for ensemble
    mod <- fit_max(
      data = data,
      response = response,
      predictors = predictors,
      predictors_f = predictors_f,
      partition = partition,
      thr = thr,
      fit_formula = NULL,
      background = background,
      clamp = clamp,
      classes = best_tune$classes,
      pred_type = pred_type,
      regmult = best_tune$regmult
    )
    pred_test_ens <- mod[["data_ens"]]

    pred_test <- data.frame(
      "pr_ab" = data[response],
      "pred" = predict_maxnet(
        mod$model,
        newdata = data,
        clamp = TRUE,
        type = pred_type
      )
    )

    threshold <- sdm_eval(
      p = pred_test$pred[pred_test$pr_ab == 1],
      a = pred_test$pred[pred_test$pr_ab == 0],
      thr = thr
    )

    result <- list(
      model = mod$model,
      predictors = variables,
      performance = dplyr::left_join(best_tune, threshold[1:4], by = "threshold") %>%
        dplyr::relocate(dplyr::all_of(hyperp), model, threshold, thr_value, n_presences, n_absences),
      hyper_performance = eval_final,
      data_ens = pred_test_ens
    )
    return(result)
  }
