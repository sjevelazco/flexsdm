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
  if(is.null(fit_formula)){
    if(poly>=2){
      forpoly <- lapply(2:poly, function(x){
        paste("I(", predictors, "^", x, ")",
              sep = "", collapse = " + ")
      }) %>% paste(collapse = " + ")
      formula1 <- paste(c(predictors, predictors_f), collapse= ' + ') %>%
        paste(., forpoly, sep =" + ")
    } else {
      formula1 <-
        paste(c(predictors, predictors_f), collapse = " + ")
    }

    if(inter_order>0){
      forinter <- c(predictors, predictors_f)
      if(inter_order>length(forinter)){
        stop('value of inter_order is higher than number of predicors ', '(', length(forinter), ')')
      }
      forinter_l <- list()

      for(i in 1:inter_order){
        forinter_l[[i]] <- do.call('expand.grid',
                            c(lapply(1:(i + 1), function(x)
                              forinter), stringsAsFactors = FALSE))
        forinter_l[[i]] <- apply(forinter_l[[i]], 1, function(x) {
          x <- unique(sort(x))
          if(length(x)>i){
            paste(x, collapse = ":" )
          }
        }) %>% unlist %>% unique

      }
      forinter <- sapply(forinter_l, paste, collapse=' + ')
      forinter <- do.call('paste', c(as.list(forinter), sep=' + '))
    }

    if(exists('forinter')){
      formula1 <- paste(formula1, forinter, sep = ' + ')
      formula1 <- stats::formula(paste(
        response, "~", formula1
      ))
    } else {
      formula1 <- stats::formula(paste(
        response, "~", formula1
      ))
    }
    message("Formula used for model fitting:\n",
            Reduce(paste, deparse(formula1)) %>% gsub(paste('  ', '   ', collapse = "|"), ' ', .))
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
                                        family = binomial)))


    # Predict for presences absences data
      if(!is.null(predictors_f)){
        for(fi in 1:length(predictors_f)){
          lev <- as.character(unique(mod[[i]]$data[,predictors_f[fi]]))
          lev_filt <- test[[i]][,predictors_f[fi]]%in%lev
          test[[i]] <- test[[i]][lev_filt,]
        }
      }

      pred_test <- try(data.frame(pr_ab = test[[i]][, response],
                                       pred = suppressWarnings(
                                         stats::predict.glm(
                                           mod[[i]],
                                           newdata = test[[i]],
                                           type = "response",
                                           se.fit = FALSE
                                         )
                                       )))

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
                                family = binomial))

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
    selected_threshold = threshold[[1]] %>% dplyr::select( threshold:TNR ),
    threshold_table = threshold[[2]] %>% dplyr::select( threshold:TNR )
  )
  return(result)
}
