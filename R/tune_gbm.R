#' Function for constructing Generalized Boosted Regression models with exploration of hyper-parameters
#'
#'
#' @param data data.frame. Database with response (0,1) and predictors values. 
#' @param response character. Column name with species absence-presence data (0,1). 
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous or discrete variables). Usage predictors = c()
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors = c()
#' @param partition character. Column name with training and validation partition groups.
#' @param grid data.frame. Provide a data frame object with algorithm hyperparameters values to be tested. It Is recommended to generate this data.frame with grid() function. In the case this argument is set as NULL. It will not perform the tuning process using the defalut values of the parameters.  
#' @param thr character. Threshold used to get binary suitability values (i.e. 0,1). It is useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. It is necessary to provide a vector for this argument. The next threshold area available:
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
#' @param metric character. Performance metric used for selecting the best combination of hyper-parameter values. Can be used one of the next metrics SORENSEN, JACCARD, FPB, TSS, KAPPA, AUC, and BOYCE. TSS is used as default.
#'   
#' @return
#' @export
#'
#' @examples
tune_gbm <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           partition,
           grid = NULL,
           thr = NULL,
           metric = 'TSS',
           ...) {
    
    
    data <- data.frame(data)
    
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
                                  paste(c(predictors, predictors_f), collapse = ' + ')))
    
    # Prepare grid when grid=default or NULL
    if(is.null(grid)){
      nv <- length(stats::na.omit(c(predictors, predictors_f)))
      grid <- data.frame(n.trees = 100, shrinkage = 0.1, n.minobsinnode = 10)
    } 
    if(class(grid)=='character'){
      nv <- length(stats::na.omit(c(predictors, predictors_f)))
      if(grid=='defalut'){
        grid <- expand.grid(
          n.trees = c(20, 50, 100, 200, 500),
          shrinkage = c(0.01, 0.1, 0.3),
          n.minobsinnode = c(1:10, 20, 30) #try define a best default tuning set
        )
      }
    }
    
    # Test hyper-parameters names
    hyperp <- names(grid)
    if(!all(c('n.trees', 'shrinkage', 'n.minobsinnode')%in%hyperp)){
      stop("Database used in 'grid' argument has to contain this columns for tunning: 'n.trees', 'shrinkage', 'n.minobsinnode'")
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
      mod <- as.list(rep(NA, nrow(grid)))
      names(mod) <- 1:nrow(grid)
      for (ii in 1:nrow(grid)) {
        set.seed(1)
        try(mod[[ii]] <-
              gbm::gbm(
                Fmula,
                data = train[[i]],
                distribution = "bernoulli",
                bag.fraction = 1, #Explore more this parameter
                n.trees = grid$n.trees[ii],
                interaction.depth = 1,
                shrinkage = grid$shrinkage[ii],
                n.minobsinnode = grid$n.minobsinnode[ii]
              )
        )
      }
      
      filt <- sapply(mod, function(x) class(x)=="gbm")
      mod <- mod[filt]
      grid2 <- grid[filt, ]
      
      # Predict for presences absences data
      pred_test <-
        lapply(mod, function(x)
          data.frame(
            pr_ab = test[[i]][,response],
            pred = gbm::predict.gbm(
              x,
              newdata = test[[i]],
              type = "response"
            )
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
      eval[hyperp]
      eval_partial[[i]] <- eval
    }
    
    # Create final database with parameter performance
    names(eval_partial) <- 1:N
    eval_partial <- eval_partial %>% dplyr::bind_rows(.,.id='partition')
    
    eval_final <- eval_partial %>% 
      dplyr::select(-partition, -c(tune:n_absences)) %>%
      dplyr::group_by_at(hyperp) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(),
                                     list(mean = mean, sd = sd)), .groups = 'drop')
    
    # Find the bets parameter setting 
    filt <- eval_final %>% dplyr::pull(paste0(metric, '_mean'))
    filt <- which.max(filt)
    best_tune <- eval_final[filt,]
    best_hyperp <- eval_final[filt,hyperp]
    
    
    # Fit final models with best settings 
    set.seed(1)
    mod <- 
      gbm::gbm(
        Fmula,
        data = data,
        distribution = "bernoulli",
        n.trees = best_hyperp$n.trees,
        interaction.depth = 1,
        shrinkage = best_hyperp$shrinkage,
        n.minobsinnode = best_hyperp$n.minobsinnode
      )
    
    
    pred_test <- data.frame(
      pr_ab = data[,response],
      pred = gbm::predict.gbm(
        mod,
        newdata = data,
        type = "response"
      ))
    
    threshold <- enm_eval(p = pred_test$pred[pred_test$pr_ab == 1],
                          a = pred_test$pred[pred_test$pr_ab == 0],
                          thr = thr)
    
    result <- list(model = mod, 
                   tune_performance=eval_final,
                   best_hyper_performance=best_tune,
                   best_hyper=best_hyperp,
                   selected_threshold=threshold[[1]],
                   threshold_table=threshold[[2]])
    return(result)
  }
