#' Function for performing support vector machine model exploring hyper-parameters 
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
    
    require(kernlab)
    require(dplyr)
    
    data <- data.frame(data)
    
    predictors <- c(predictors, predictors_f)
    if (is.null(predictors_f)) {
      data <- data[, c(response, predictors, partition)]
      data <- data.frame(data)
    } else {
      data <- data[, c(response, predictors, predictors_f, partition)]
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
    if(any(!hyperp%in%c("C", "sigma"))){
      stop("Database used in 'grid' argument has to contain this columns for tunning: 'C', 'sigma'")
    }
    
    grid$tune <- 1:nrow(grid)
    
    
    N <- max(data[partition])
    
    train <- list()
    test <- list()
    for (i in 1:N) {
      train[[i]] <- data[data[, partition] == i, ]
      test[[i]] <- data[data[, partition] != i, ]
    }
    
    eval_partial <- list()
    for(i in 1:N){
      message('Partition number: ', i, '/', N)
      mod <- list()
      for (ii in 1:nrow(grid)) {
        try(mod[[ii]] <-
              kernlab::ksvm(
                Fmula,
                data = train[[i]],
                type = "C-svc",
                kernel = "rbfdot",
                kpar = list(sigma = grid$sigma[ii]),
                C = grid$C[ii],
                prob.model = TRUE
              )
        )
      }
      
      filt <- !sapply(mod, is.null)
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
            enm_eval(p = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 1],
                     a = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 0],
                     thr = thr)
      }
      
      eval <- dplyr::bind_rows(lapply(eval, function(x) x$selected_threshold))
      eval <- dplyr::tibble(cbind(grid2, eval)) 
      eval_partial[[i]] <- eval
    }
    
    # Create final database with parameter performance
    names(eval_partial) <- 1:N
    eval_partial <- eval_partial %>% dplyr::bind_rows(.,.id='partition')
    
    eval_final <- eval_partial %>% 
      dplyr::select(-partition, -c(tune:n_absences)) %>%
      dplyr::group_by_at(hyperp) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(),
                                     list(mean = mean, sd = sd)), .groups = "drop")
    
    # Find the bets parameter setting 
    filt <- eval_final %>% dplyr::pull(paste0(metric, '_mean'))
    filt <- which.max(filt)
    best_tune <- eval_final[filt,]
    best_hyperp <- eval_final[filt,hyperp]
    
    
    # Fit final models with best settings 
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
  
    result <- list(model = mod, 
                   # eval_partial,
                   tune_performance=eval_final,
                   best_hyperparameter=best_hyperp,
                   selected_threshold=threshold[[1]],
                   threshold_table=threshold[[2]])
    return(result)
  }
