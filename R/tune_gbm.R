#' Function for constructing Random Forest with exploration of hyper-parameters
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
#' @importFrom dplyr bind_rows tibble select group_by_at summarise across everything pull
#' @importFrom gbm gbm predict.gbm
#' @importFrom stats formula na.omit
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
    
    require(gbm)
    require(dplyr)
    
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
