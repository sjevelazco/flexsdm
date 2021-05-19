fit_ensemble <- function(models, thr,  metric) {
  #### Performance metric
  metric <- paste0(metric, '_mean')

  #### Models names
  nms <- paste0("m_", 1:length(models))

  #### Model performances
  perf <- sapply(models, function(x) {
    x[["performance"]] %>%
      dplyr::filter(threshold == dplyr::all_of(thr)) %>%
      dplyr::pull(dplyr::all_of(metric))
  })

  #### Model thresholds
  thr <- sapply(models, function(x) {
    x[["selected_thresholds"]] %>%
      dplyr::filter(threshold == dplyr::all_of(thr)) %>%
      dplyr::pull(terra::values)
  })

  names(perf) <- names(thr) <- nms

  #### Extract and merge suitability databases
  data_ens <- sapply(models, function(x)
    x['data_ens'])

  data_ens <- mapply(function(x, cn) {
    colnames(x)[colnames(x) %in% "pred"] <- cn
    x
  }, data_ens, nms, SIMPLIFY = FALSE)

  data_ens <-lapply(data_ens, function(x)
    x %>% dplyr::mutate(pr_ab = pr_ab %>%
                          as.character() %>%
                          as.double()))

  data_ens2 <-
    dplyr::inner_join(data_ens[[1]],
                      data_ens[[2]],
                      by = c("rnames", "replicates", "part", "pr_ab"))
  if (length(data_ens) > 2) {
    for (i in 3:length(data_ens)) {
      data_ens2 <-
        dplyr::inner_join(data_ens2,
                          data_ens[[2]],
                          by = c("rnames", "replicates", "part", "pr_ab"))
    }
  }
  rm(data_ens)

  #### Extract predicted suitability of each model
  values <- data_ens2 %>%
    dplyr::select(dplyr::starts_with('m_'))

  #### Remove suitability values from data_ens2
  data_ens2 <- data_ens2 %>% dplyr::select(-dplyr::starts_with('m_'))

  #### Perform ensembles
  data_ens2 <- data_ens2 %>% dplyr::mutate(

    mean = apply(values, 1, function(x)
      mean(x, na.rm = TRUE)),

    meanw = mapply(function(x, v)
      (x * v), values, perf, SIMPLIFY = TRUE) %>%
      apply(., 1, function(x)
        mean(x, na.rm = TRUE)),

    meansup = apply(values[, perf >= mean(perf)], 1, function(x)
      mean(x, na.rm = TRUE)),

    meanthr = mapply(function(x, v)
      ifelse(x >= v, x, 0), values, thr, SIMPLIFY = TRUE) %>%
      apply(., 1, function(x)
        mean(x, na.rm = TRUE)),

    median = apply(values, 1, function(x)
      median(x, na.rm = TRUE))
  )
  rm(values)

  #### Calculate ensemble performance
  data_ens2

  p_names <- data_ens2 %>% dplyr::pull(replicates) %>% unique()
  np <- length(p_names)
  eval_partial_list <- list()

  for (h in 1:np) {
    message("Replica number: ", h, "/", np)

    pred_test <- data_ens2 %>%
      dplyr::filter(replicates == p_names[h]) %>%
      split(., f = part)
    np2 <- length(test)
    eval_partial <- list()

    for (i in 1:np2) {
      # Validation of model
      eval <-
        sdm_eval(
          p = pred_test$pred[pred_test$pr_ab == 1],
          a = pred_test$pred[pred_test$pr_ab == 0],
          thr = thr
        )
      if (is.null(thr)) {
        eval_partial[[i]] <- eval$all_thresholds
      } else {
        eval_partial[[i]] <- eval$selected_thresholds
      }
    }

    # Create final database with parameter performance
    names(eval_partial) <- 1:np2
    eval_partial <- eval_partial %>%
      dplyr::bind_rows(., .id = "partition")
    eval_partial_list[[h]] <- eval_partial
  }

}



