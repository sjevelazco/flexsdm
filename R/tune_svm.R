#' Function for constructing Support Vector Machine with exploration of hyper-parameters
#'
#'
#'
#' @param data
#' @param response
#' @param predictors
#' @param predictors_f
#' @param partition
#' @param grid
#' @param thr
#' @param metric
#'
#' @importFrom dismo predict
#' @importFrom dplyr bind_rows tibble select group_by_at summarise across everything pull
#' @importFrom kernlab ksvm
#' @importFrom stats formula
#'
#' @return
#' @export
#'
#' @examples
tune_svm <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           partition,
           grid = NULL,
           thr = NULL,
           metric = "TSS",
           ...) {
    data <- data.frame(data)

    # Transform response variable as factor
    data[,response] <- as.factor(data[,response])

    if (is.null(predictors_f)) {
      data <- data %>%
        dplyr::select(response, predictors, dplyr::starts_with(partition))
      data <- data.frame(data)
    } else {
      data <- data %>%
        dplyr::select(response, predictors, predictors_f, dplyr::starts_with(partition))
      data <- data.frame(data)
      for (i in predictors_f) {
        data[, i] <- as.factor(data[, i])
      }
    }

    # Formula
    Fmula <- stats::formula(paste(response, '~',
                                  paste(c(predictors, predictors_f), collapse = " + ")))

    # Prepare grid when grid=default or NULL
    if(is.null(grid)){
      grid <- data.frame(C=1, sigma = "automatic")
    }
    if(class(grid)=='character'){
      if(grid=='defalut'){
        grid <-  expand.grid(C = c(1, 2, 4, 8, 16),
                             sigma = c(0.001, 0.01, 0.1, 0.2))
      }
    }

    # Test hyper-parameters names
    hyperp <- names(grid)
    if(!all(c("C", "sigma")%in%hyperp)){
      stop("Database used in 'grid' argument has to contain this columns for tunning: 'C', 'sigma'")
    }

    grid$tune <- 1:nrow(grid)


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

      for (i in 1:np2) {
        message("Partition number: ", i, "/", np2)
        mod <- as.list(rep(NA, nrow(grid)))
        names(mod) <- 1:nrow(grid)
        for (ii in 1:nrow(grid)) {
          set.seed(1)
          try(capture.output({
            mod[[ii]] <-
              kernlab::ksvm(
                Fmula,
                data = train[[i]],
                type = "C-svc",
                kernel = "rbfdot",
                kpar = list(sigma = grid$sigma[ii]),
                C = grid$C[ii],
                prob.model = TRUE
              )
          })
          )
        }


        filt <- sapply(mod, function(x) class(x) == "ksvm")
        mod <- mod[filt]
        grid2 <- grid[filt, ]

        # Predict for presences absences data
        pred_test <-
          lapply(mod, function(x)
            data.frame(
              pr_ab = test[[i]][,response],
              pred = dismo::predict(
                x,
                newdata = test[[i]],
                type = 'prob'
              )[,2]
            ))

        # Validation of parameter combination
        eval <- list()
        for (ii in 1:length(pred_test)) {
          eval[[ii]] <-
            enm_eval(
              p = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 1],
              a = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 0],
              thr = thr
            )
        }

        eval <- dplyr::bind_rows(lapply(eval, function(x) x$selected_threshold))
        eval <- dplyr::tibble(cbind(grid2, eval))
        eval_partial[[i]] <- eval
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
      dplyr::select(-replica, -partition, -c(tune:n_absences)) %>%
      dplyr::group_by_at(hyperp) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(),
                                     list(mean = mean, sd = sd)), .groups = "drop")

    # Find the bets parameter setting
    filt <- eval_final %>% dplyr::pull(paste0(metric, '_mean'))
    filt <- which.max(filt)
    best_tune <- eval_final[filt,]
    best_hyperp <- eval_final[filt,hyperp]


    # Fit final models with best settings
    set.seed(1)
    mod <-
      kernlab::ksvm(
        Fmula,
        data = data,
        type = "C-svc",
        kernel = "rbfdot",
        kpar = list(sigma = best_hyperp$sigma),
        C = best_hyperp$C,
        prob.model = TRUE
      )

    pred_test <- data.frame(
      pr_ab = data[,response],
      pred = dismo::predict(
        mod,
        newdata = data,
        type = "prob"
      )[,2])

    threshold <- enm_eval(p = pred_test$pred[pred_test$pr_ab == 1],
                          a = pred_test$pred[pred_test$pr_ab == 0],
                          thr = thr)

    result <- list(
      model = mod,
      tune_performance = eval_final,
      best_hyper_performance = best_tune,
      selected_threshold = threshold[[1]],
      threshold_table = threshold[[2]]
    )
    return(result)
    }
