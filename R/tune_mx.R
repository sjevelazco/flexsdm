#' Function for performing maxent model tuning hyperparameters 
#'
#'
#' @param data 
#' @param response 
#' @param predictors 
#' @param predictors_f 
#' @param background 
#' @param partition 
#' @param grid 
#' @param thr 
#' @param metric 
#' @param clamp 
#' @param pred_type 
#'
#' @importFrom dismo predict
#' @importFrom dplyr filter pull bind_rows tibble select group_by_at summarise across everything
#' @importFrom maxnet maxnet maxnet.formula
#'  
#' @return
#' @export
#'
#' @examples
tune_mx <-
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
           pred_type = "cloglog") {
    
    require(maxnet)
    require(dplyr)
    
    if (is.null(predictors_f)) {
      data <- data[, c(response, predictors, partition)]
      data <- data.frame(data)
      if(!is.null(background)){
        background <- data[, c(response, predictors, partition)]
        background <- data.frame(background)
      }
    } else {
      data <- data[, c(response, predictors, predictors_f, partition)]
      data <- data.frame(data)
      for(i in predictors_f){
        data[,i] <- as.factor(data[,i])
      }
      if(!is.null(background)){
        background <- data[, c(response, predictors, partition)]
        background <- data.frame(background)
        for(i in predictors_f){
          background[,i] <- as.factor(background[,i])
        }
      }
    }
    
    if(!is.null(background)){
      if(!all(names(data)%in%names(background))){
        stop("Column names of database used in 'data' and background arguments do not match")
      }
    }
    
    if(is.null(grid)){
      grid <- data.frame(regmult=1, classes = "default")
    } 
    if(class(grid)=='character'){
      if(grid=='defalut'){
        grid <- expand.grid(regmult = seq(0.1, 3, 0.5),
                            classes = c("l", "lq", "lqh", "lqhp", "lqht"))
      }
    }
    hyperp <- names(grid)
    if(any(!hyperp%in%c("regmult", "classes"))){
      stop("Database used in 'grid' argument has to contain this columns for tunning: 'regmult', 'classes'")
    }
    
    grid$tune <- 1:nrow(grid)
    
    
    N <- max(data[partition])
    if(!is.null(background)){
      Npart_p <- data %>% dplyr::filter(!!as.symbol(response)==1) %>% 
        dplyr::pull({{partition}}) %>% unique %>% sort
      Npart_bg <- background %>% dplyr::filter(!!as.symbol(response)==1) %>% 
        dplyr::pull({{partition}}) %>% unique %>% sort
      if(!all(Npart_p %in% Npart_bg)) {
        stop(
          paste(
            "Partition groups between presences and backgroud do not match:\n",
            paste("Part. group presences:", paste(Npart_p, collapse = ' '), "\n"),
            paste("Part. group background:", paste(Npart_bg, collapse = ' '), "\n")
          )
        )
      }
    }
  
  train <- list()
  test <- list()
  for (i in 1:N) {
    train[[i]] <- data[data[, partition] == i, ]
    test[[i]] <- data[data[, partition] != i, ]
  }
  # In the follow code function will substitutes absences by background points
  # only in train database in order to fit maxent with presences and background
  # and validate models with presences and absences
  if(!is.null(background)){
    train <- lapply(train, function(x) x[x[,response]==1,])
    bgt <- split(background, background[,partition])
    train <- mapply(dplyr::bind_rows, train, bgt, SIMPLIFY = FALSE)
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
    grid2 <- grid[filt, ]
    
    pred_test <-
      lapply(mod, function(x)
        data.frame(
          'pr_ab' = test[[i]][response],
          'pred' = dismo::predict(
            x,
            newdata = test[[i]],
            clamp = clamp,
            type = pred_type
          )
        ))
    
    eval <- list()
    for(ii in 1:length(pred_test)) {
      eval[[ii]] <-
        enmtml_evaluate(p = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 1],
                        a = pred_test[[ii]]$pred[pred_test[[ii]]$pr_ab == 0],
                        thr = thr)
    }
    
    eval <- dplyr::bind_rows(lapply(eval, function(x) x$selected_threshold))
    eval <- dplyr::tibble(cbind(grid2, eval)) 
    eval_partial[[i]] <- eval
  }
  
  names(eval_partial) <- 1:N
  eval_partial <- eval_partial %>% dplyr::bind_rows(.,.id='partition')
  
  # n_part_f <- table(eval_partial$partition)
  # if(length(unique(n_part_f))>1){
  #   n_part_f <- eval_partial %>% dplyr::group_by_at(hyperp) %>% dplyr::count()
  # }
  eval_final <- eval_partial %>% 
    dplyr::select(-partition, -c(tune:n_absences)) %>%
    dplyr::group_by_at(hyperp) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(),
                                   list(mean = mean, sd = sd)), .groups = "drop")
  
  
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
    'pred' = dismo::predict(
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
       selected_threshold=threshold[[1]],
       threshold_table=threshold[[2]])
  return(result)
  }