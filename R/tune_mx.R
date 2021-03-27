# source("./R/boyce.R")
# source("./R/enmtml_evaluate.R")

tune_mx <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           partition,
           grid = NULL,
           thr = NULL,
           metric = "TSS") {
    
    require(maxnet)
    require(dplyr)
    
    if (is.null(predictors_f)) {
      data <- data[, c(response, predictors, partition)]
      data <- data.frame(data)
    } else {
      data <- data[, c(response, predictors, predictors_f, partition)]
      data <- data.frame(data)
    }
  
    hyperp <- names(grid)
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
            maxnet::maxnet(
              p = train[[i]][,response],
              data = train[[i]][predictors],
              f = maxnet::maxnet.formula(train[[i]][response], 
                                         train[[i]][predictors], 
                                         classes = grid$classes[ii]),
              regmult = grid$regmult[ii]
            ))
    }
    
    filt <- !sapply(mod, is.null)
    mod <- mod[filt]
    grid <- grid[filt, ]
    
    pred_test <-
      lapply(mod, function(x)
        data.frame(
          'pr_ab' = test[[i]][response],
          'pred' = predict(
            x,
            newdata = test[[i]],
            clamp = TRUE,
            type = 'cloglog'
          )
        ))
    
    eval <- list()
    for(ii in 1:length(pred_test)) {
      eval[[ii]] <-
        enmtml_evaluate(p = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 1],
                        a = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 0],
                        thr = thr)
    }
    
    eval <- dplyr::bind_rows(eval)
    eval <- dplyr::tibble(cbind(grid, eval)) 
    eval_partial[[i]] <- eval
  }
  
  names(eval_partial) <- 1:N
  eval_partial <- eval_partial %>% dplyr::bind_rows(.,.id='partition')
  
  eval_final <- eval_partial %>% 
    dplyr::select(-partition, -c(tune:n_absences)) %>%
    dplyr::group_by_at(hyperp) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(),
                                   list(mean = mean, sd = sd))) %>% 
    dplyr::group_by()
  
  
  filt <- eval_final %>% dplyr::pull(paste0(metric, '_mean'))
  filt <- which.max(filt)
  best_tune <- eval_final[filt,]
  best_hyperp <- eval_final[filt,hyperp]
  
  
  # Fit final models with best settings 
  best_hyperp
  
  mod <-
    maxnet::maxnet(
      p = data[,response],
      data = data[predictors],
      f = maxnet::maxnet.formula(data[response], 
                                 data[predictors], 
                                 classes = best_hyperp$classes),
      regmult = best_hyperp$regmult
    )
  
  pred_test <- data.frame(
    'pr_ab' = data[response],
    'pred' = predict(
      mod,
      newdata = data,
      clamp = TRUE,
      type = 'cloglog'
    ))
    
  threshold <- enmtml_evaluate(p = pred_test$pred[pred_test$pr_ab == 1],
                               a = pred_test$pred[pred_test$pr_ab == 0],
                               thr = thr)
  
  
  result <- list(model = mod, 
       # eval_partial,
       tune_performance=eval_final,
       best_hyperparameter=best_hyperp,
       threshold=threshold)
  return(result)
  }