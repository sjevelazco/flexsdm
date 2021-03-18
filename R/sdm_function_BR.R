##%######################################################%##
#                                                          #
####                SDM Function                        ####
#                                                          #
##%######################################################%##
# Written by Brooke Rose
#' Title
#'
#' @param df  
#' @param pr_ab 
#' @param env_preds 
#' @param sp_area 
#' @param pred_rasters 
#' @param species_name 
#' @param dir_save 
#' @param cores
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
sdms <- function(df, # full data set
                 pr_ab, # response variable (0/1)
                 env_preds, # list of characters for the environmental predictors 
                 sp_area, # polygon for delineating prediction area
                 pred_rasters, # prediction rasters
                 species_name, # character of species name
                 dir_save, # directory for saving model objects and spatial predictions
                 cores) 
{ 
  
  require(raster)
  require(dplyr)
  require(vroom)
  require(ENMTML)
  require(ecospat)
  require(e1071)
  require(tidyverse)
  require(dismo)
  require(sf)
  require(gam)
  require(randomForest)
  require(MASS)
  require(biomod2)
  require(caret)
  require(doParallel)
  require(parallel)
  require(nnet)
  require(rasterVis)
  if (!"devtools"%in%installed.packages()){devtools::install_github("babaknaimi/sdm")}  
  if (!"devtools"%in%installed.packages()){install.packages("devtools")}  
  # devtools::install_github("andrefaa/ENMTML") 
  require(sdm)
  require('ENMTML')
  
  if(!dir.exists(paste0(dir_save, 'models/', sep = ''))) {
    dir.create(paste0(dir_save, 'models/', sep = ''))
  }
  
  df_sp <-
    st_as_sf(
      df,
      coords = c('x_albers', 'y_albers'),
      remove = FALSE,
      crs = crs(env_stack)
    )
  
  df_extract <-
    raster::extract(pred_rasters, df_sp, df = TRUE, sp = TRUE) %>%
    st_as_sf(remove = FALSE)
  
  df_extract$terrain <- as.integer(df_extract$terrain)
  
  # data frames with only response and predictor variables
  df_clean <-
    df_extract %>% dplyr::select(all_of(env_preds), pr_ab) %>%
    st_set_geometry(NULL) %>%
    na.omit()
  
  df_clean_sp <-
    df_extract %>%
    na.omit()
  
  p <- df_extract %>% dplyr::select(all_of(env_preds), pr_ab) %>%
    filter(pr_ab ==1)
  a <- df_extract %>% dplyr::select(all_of(env_preds), pr_ab) %>%
    filter(pr_ab ==0)
  
  # preparing prediction rasters
  pred_crop <- raster::crop(pred_rasters, sp_area)
  pred_mask <- raster::mask(pred_crop, sp_area)
  
  
  message('Number of cores: ', detectCores())
  message('Cores used: ', cores)
  # cl <- makePSOCKcluster(cores)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  
  # seting trainControl function for tuning models with caret
  fit_control <- caret::trainControl(
    # method = "repeatedcv",## 10-fold CV
    method = "cv",
    number = 10, ## number of folds
    # repeats = 5, ## for repeating fold CV 
    selectionFunction = "best",
    classProbs = TRUE, ## Estimate class probabilities
    summaryFunction = caret::twoClassSummary,
  ) 
  
  # Testing and Training data for models 
  set.seed(123)
  df_idx = createDataPartition(df_clean$pr_ab, p = 0.7, list = FALSE)
  df_trn = df_clean[df_idx, ]
  df_tst = df_clean[-df_idx, ]
  
  
  
  ##%######################################################%##
  #                                                          #
  ####             Generalized Linear Models              ####
  #                                                          #
  ##%######################################################%##
  #### Final model built with all data
  glm.formula <-
    makeFormula("pr_ab", df_clean[, env_preds],
                "quadratic", interaction.level = 1)
  
  glm_final <-
    stepAIC(glmStart <- glm(pr_ab ~ 1,
                            data = df_clean,
                            family = binomial),
            glm.formula,
            data = df_clean,
            direction = "both",
            trace = FALSE,
            k = 2,
            control = glm.control(maxit = 100))
  
  glmVarImp <- caret::varImp(glm_final, scale = FALSE)
  
  # training and testing model
  glm_train <-
    stepAIC(glmStart <- glm(pr_ab ~ 1,
                            data = df_trn,
                            family = binomial),
            glm.formula,
            data = df_trn,
            direction = "both",
            trace = FALSE,
            k = 2,
            control = glm.control(maxit = 100))
  
  glm_pred <- predict(glm_train, df_tst, type = "response")
  glm_e <- sdm::evaluates(x = df_tst$pr_ab, p = glm_pred)
  
  ##%######################################################%##
  #                                                          #
  ####            Generalized Additive Models             ####
  #                                                          #
  ##%######################################################%##
  #### Final model built with all data
  # gam.formula <- paste("s(", env_preds, ")",sep="") #",k=3)",
  # WE NEED TO FIND A WAY TO CREAT THE FORMULA AUTOMATICALLY without writing variable names
  
  gam_final <- gam(
    pr_ab ~ s(cwd) + s(aet) + s(tmn) + s(ppt_djf) + s(ppt_jja) +
      s(ph) + s(awc) + s(depth) + s(pct_clay) + terrain,
    data = df_clean,
    family = "binomial"
  )
  
  gamVarImp <- caret::varImp(gam_final, scale = FALSE)
  
  # saveRDS(gam_final, file = file.path(dir_save, paste(species_name, "_gam.rda", sep="")))
  
  
  #### Training and testing model
  gam_train <- gam(
    pr_ab ~ s(cwd) + s(aet) + s(tmn) + s(ppt_djf) + s(ppt_jja) +
      s(ph) + s(awc) + s(depth) + s(pct_clay) + terrain,
    data = df_trn,
    family = "binomial"
  )
  
  gam_pred <- predict(gam_train, df_tst, type = "response")
  gam_e <- sdm::evaluates(x = df_tst$pr_ab, p = gam_pred)
  
  ##%######################################################%##
  #                                                          #
  ####                   Random Forest                    ####
  #                                                          #
  ##%######################################################%##
  
  #### Final model built with all data
  # Find best parameters for final model
  
  tune_grid <- expand.grid(mtry = seq(2, length(env_preds), 1))
  
  set.seed(123)
  best_tune <- caret::train(
    pr_ab ~ .,
    data = dplyr::mutate(df_clean,
                         pr_ab =
                           as.factor(ifelse(
                             pr_ab == 1, 'Pres', 'Abs'
                           ))),
    method = "rf",
    type = 'classification',
    verbose = FALSE,
    metric = "ROC",
    trControl = fit_control,
    tuneGrid = tune_grid
  )
  # ggplot(best_tune)
  
  # model built with all data
  set.seed(123)
  rf_final <- randomForest::randomForest(
    pr_ab ~ .,
    data = dplyr::mutate(df_clean, pr_ab = as.factor(pr_ab)),
    ntree = 1000,
    importance = TRUE,
    type = 'classification',
    mtry = best_tune$bestTune$mtry
  )
  
  #saveRDS(rf_final, file = file.path(dir_save, 'models/rf_final.rda'))
  
  rfVarImp <- randomForest::importance(rf_final)
  
  #### Training and testing model
  # with best tuning 
  set.seed(123)
  rf_train <- randomForest::randomForest(
    x = df_trn[, env_preds],
    y = as.factor(df_trn$pr_ab),
    ntree = 1000,
    importance = TRUE, 
    type='classification',
    mtry=best_tune$bestTune$mtry
  )
  
  rf_pred <- predict(rf_train, df_tst, type = "prob")[, 2]
  rf_e <- sdm::evaluates(x = df_tst$pr_ab, p = rf_pred)
  
  ##%######################################################%##
  #                                                          #
  ####              Boosted Regression Tree               ####
  #                                                          #
  ##%######################################################%##
  #### Final model built with all data
  # Find best parameters for final model
  tune_grid <- expand.grid(
    interaction.depth = seq(2, 10, 20),
    n.trees = c(100, 200, 500, 1000, 2000),
    shrinkage = 0.1,
    n.minobsinnode = 20
  )
  
  set.seed(123)
  best_tune <- caret::train(
    pr_ab ~ .,
    data =
      dplyr::mutate(df_clean,
                    pr_ab =
                      as.factor(ifelse(pr_ab == 1, 'Pres', 'Abs'))),
    method = "gbm",
    verbose = FALSE,
    metric = "ROC",
    trControl = fit_control,
    tuneGrid = tune_grid
  )
  
  # ggplot(best_tune)
  
  set.seed(123)
  brt_final <-
    gbm::gbm(
      pr_ab ~ .,
      data = df_clean,
      distribution = "bernoulli",
      n.trees = best_tune$bestTune$n.trees,
      interaction.depth = best_tune$bestTune$interaction.depth,
      shrinkage = best_tune$bestTune$shrinkage,
      n.minobsinnode = best_tune$bestTune$n.minobsinnode
    )
  
  # saveRDS(brt_final, file = file.path(dir_save, 'models/brt_final.rda'))
  
  # variable importance
  brtVarImp <- caret::varImp(brt_final, scale = FALSE, 
                             numTrees = best_tune$bestTune$n.trees)
  
  #### Training and testing model
  # with best tuning
  set.seed(123)
  brt_train <-
    gbm::gbm(pr_ab ~ .,
             data = df_trn,
             distribution = "bernoulli",
             n.trees = best_tune$bestTune$n.trees,
             interaction.depth = best_tune$bestTune$interaction.depth,
             shrinkage = best_tune$bestTune$shrinkage,
             n.minobsinnode = best_tune$bestTune$n.minobsinnode
    )
  
  
  brt_pred <- predict(brt_train, df_tst, type = 'response')
  brt_e <- sdm::evaluates(x = df_tst$pr_ab, p = brt_pred)
  
  ##%######################################################%##
  #                                                          #
  ####              Support Vector Machines               ####
  #                                                          #
  ##%######################################################%##
  
  #### Final model built with all data
  
  # Combination of parameters values to be tested
  tune_grid <-
    expand.grid(C = c(1, 2, 4, 8, 16),
                sigma = c(0.001, 0.01, 0.1, 0.2))
  
  set.seed(123)
  best_tune <- caret::train(
    pr_ab ~ .,
    data =
      dplyr::mutate(df_clean,
                    pr_ab =
                      as.factor(ifelse(
                        pr_ab == 1, 'Pres', 'Abs'
                      ))),
    method = "svmRadial",# Radial kernel (kernel = "rbfdot")
    metric = "ROC",
    trControl = fit_control,
    tuneGrid = tune_grid
  )
  
  # ggplot(best_tune)
  
  set.seed(123)
  svm_final <- kernlab::ksvm(
    pr_ab ~ .,
    data = df_clean,
    type = "C-svc",
    kernel = "rbfdot",
    kpar = list(sigma = best_tune$bestTune$sigma),
    C = best_tune$bestTune$C,
    prob.model = TRUE
  )
  
  # saveRDS(svm_final, file = file.path(dir_save, 'models/svm_final.rda'))
  
  pred_svm <- predict(svm_final, df_clean, type = 'prob')[, 2]
  svmVarImp <- caret::filterVarImp(df_clean[, env_preds], pred_svm, nonpara = FALSE)
  
  #### Training and testing model
  set.seed(123)
  svm_train <- kernlab::ksvm(
    pr_ab ~ .,
    data = df_tst,
    type = "C-svc",
    kernel = "rbfdot",
    kpar = list(sigma = best_tune$bestTune$sigma),
    C = best_tune$bestTune$C,
    prob.model = TRUE
  )
  
  
  svm_pred <- predict(svm_train, df_tst, type = "prob")[, 2]
  svm_e <- sdm::evaluates(x = df_tst$pr_ab, p = svm_pred)
  
  ##%######################################################%##
  #                                                          #
  ####            Artificial Neural Networks              ####
  #                                                          #
  ##%######################################################%##
  
  #### Final model built with all data
  # Find best parameters for final model
  # Combination of parameters values to be tested
  tune_grid <-
    expand.grid(size = c(2, 4, 6, 8, 10), 
                decay = c(0.001, 0.01, 0.05, 0.1))
  
  cl <- makePSOCKcluster(cores)
  registerDoParallel(cl)
  
  set.seed(123)
  best_tune <- caret::train(
    pr_ab ~ .,
    data =
      dplyr::mutate(df_clean,
                    pr_ab =
                      as.factor(ifelse(
                        pr_ab == 1, 'Pres', 'Abs'
                      ))),
    method = "nnet",
    # Radial kernel
    metric = "ROC",
    trControl = fit_control,
    tuneGrid = tune_grid
  )
  
  #ggplot(best_tune)
  
  # final model
  # cv_nnet_final <- biomod2:::.CV.nnet(Input = df_clean[, env_preds],
  # Target = df_clean$pr_ab)
  set.seed(123)
  nnet_final <- nnet::nnet(
    pr_ab ~ .,
    data = df_clean,
    size = best_tune$bestTune$size,
    rang = 0.1,
    decay = best_tune$bestTune$decay,
    maxit = 200,
    trace = FALSE
  )
  
  #saveRDS(nnet_final, file = file.path(dir_save, 'models/nnet_final.rda'))
  
  # variable importance
  nnetVarImp <- caret::varImp(nnet_final)
  
  #### Training and testing model
  # with best tuning
  nnet_train <- nnet::nnet(
    pr_ab ~ .,
    data = df_trn,
    size = best_tune$bestTune$size,
    rang = 0.1,
    decay = best_tune$bestTune$decay,
    maxit = 200,
    trace = FALSE
  )
  
  
  # model predictions on evaluation data
  nnet_pred <- predict(nnet_train, df_tst[, env_preds])
  
  # model evaluation for model built on evaluation data
  nnet_e <- sdm::evaluates(x = df_tst$pr_ab, p = nnet_pred)
  
  # STOP CLUSTER
  stopCluster(cl)
  
  ##%######################################################%##
  #                                                          #
  ####                  post-processing                   ####
  #                                                          #
  ##%######################################################%##
  
  
  ### AUC based on evaluation model
  auc <- list(
    glm_e@statistics$AUC,
    gam_e@statistics$AUC,
    rf_e@statistics$AUC,
    brt_e@statistics$AUC,
    svm_e@statistics$AUC,
    nnet_e@statistics$AUC
  )
  names(auc) <- c('GLM', 'GAM', 'RF', 'BRT', 'SVM', 'NNET')
  
  ### spatial predictions
  beginCluster()
  
  raw_preds <- raster::brick(
    clusterR(pred_mask, raster::predict, args=list(model = glm_final, type = "response")),
    clusterR(pred_mask, raster::predict, args=list(model = gam_final, type = "response")),
    clusterR(pred_mask, raster::predict, args=list(model = rf_final, type = "prob", index = 2)),
    clusterR(pred_mask, raster::predict, args=list(model = brt_final, n.trees = brt_final$n.trees, type = "response")),
    clusterR(pred_mask, raster::predict, args=list(model = svm_final, type = "prob", index = 2)))
  
  names(raw_preds) <- c('glm', 'gam', 'rf', 'brt', 'svm')
  
  ### Ensembles
  
  ensemble <- 
    raster::brick(
      clusterR(raw_preds, fun = calc, args = list(fun = mean)),
      clusterR(raw_preds, fun = weighted.mean, args = list(w=unlist(auc))),
      clusterR(raw_preds, fun = calc, args = list(fun = sd)))
  
  names(ensemble) <- c('mean', 'weighted average', 'standard deviation')
  
  endCluster()
  
  ensemble_extract <- raster::extract(ensemble, df_clean_sp, df = TRUE, sp = TRUE)
  
  ensemble_data <- st_as_sf(ensemble_extract) %>%
    dplyr::select(pr_ab, mean, weighted.average, standard.deviation) %>%
    na.omit()
  
  mean_e <- sdm::evaluates(x = ensemble_data$pr_ab, p = ensemble_data$mean)
  w_avg_e <- sdm::evaluates(x = ensemble_data$pr_ab, p = ensemble_data$weighted.average)
  
  auc_final <- list(
    glm_e@statistics$AUC,
    gam_e@statistics$AUC,
    rf_e@statistics$AUC,
    brt_e@statistics$AUC,
    svm_e@statistics$AUC,
    #nnet_e@statistics$AUC,
    mean_e@statistics$AUC,
    w_avg_e@statistics$AUC
  )
  
  names(auc_final) <- c('GLM',
                        'GAM',
                        'RF',
                        'BRT',
                        'SVM',
                        #'NNET',
                        'ensemble_mean',
                        'ensemble_weighted_avg')
  
  # saveRDS(auc_final, file = file.path(dir_save, 'models/final_auc.rda'))
  
  
  ## all rasters
  all_raw <- raster::stack(raw_preds, ensemble)
  
  # writeRaster(all_raw, filename = file.path(dir_save, 'models/raw_pred_maps.grd'), overwrite = TRUE)
  
  ### Thresholds for based on final models
  
  # list of threshold data frames
  threshold <- list(
    glm_e@threshold_based,
    gam_e@threshold_based,
    rf_e@threshold_based,
    brt_e@threshold_based,
    svm_e@threshold_based,
    #nnet_e@threshold_based,
    mean_e@threshold_based,
    w_avg_e@threshold_based
  )
  
  names(threshold) <- c('GLM',
                        'GAM',
                        'RF',
                        'BRT',
                        'SVM',
                        #'NNET',
                        'ensemble_mean',
                        'ensemble_weighted_avg')
  
  #sapply(names(threshold),
  #      function (x)
  #       utils::write.table(
  #        threshold[[x]],
  #        file.path(dir_save, 'models/', paste(x,"_thresholds.txt")),
  #       sep = "\t",
  #      row.names = F
  #    ))
  
  sens_spec <- list() # initiating list for sensitivity = specificity threshold maps
  
  # function that creates rasters with everything less than threshold = 0 and retain values above 0
  for (i in 1:length(threshold)) {
    sens_spec[i] <-
      calc(
        all_raw[[i]],
        fun = function(x) {
          x[x < threshold[[i]][1, 2]] <- 0
          return(x)
        }
      )
  }
  
  sens_spec <- raster::stack(sens_spec)
  names(sens_spec) <- c('GLM',
                        'GAM',
                        'RF',
                        'BRT',
                        'SVM',
                        #'NNET',
                        'ensemble_mean',
                        'ensemble_weighted_avg')
  
  sp_dir <- file.path(dir_save, paste(species_name))
  dir.create(sp_dir)
  
  writeRaster(
    sens_spec,
    file.path(sp_dir, names(sens_spec)),
    bylayer = TRUE,
    format = 'GTiff'
  )
  
  
  ##### PDF Output #####
  auc_final_df <- as.data.frame(auc_final) %>%
    dplyr::rename(mean = ensemble_mean,
                  w_average = ensemble_weighted_avg) %>%
    pivot_longer(cols = 1:7,
                 names_to = "Model",
                 values_to = "AUC") 
  
  
  auc_plot <-
    ggplot(auc_final_df, aes(
      x = Model,
      y = AUC,
      label = round(AUC, 5),
    )) +
    geom_point(size = 5) + geom_label(size = 7) +
    labs(
      title = paste0(species_name, ": AUC Comparison"),
      x = "Model Type",
      y = "AUC"
    ) +
    theme(text = element_text(size = 17, family = "serif"),
          axis.title.x = element_text(vjust = .25))
  
  presence <- df_clean_sp %>%
    filter(pr_ab == 1)
  
  absence <- df_clean_sp %>%
    filter(pr_ab == 0)
  
  p.points <- as(presence, 'Spatial')
  a.points <- as(absence, 'Spatial')
  cfp.pol <- as(sp_area, 'Spatial')
  
  myTheme <- rasterTheme(region = rev(terrain.colors(7)))
  
  pa_map <- pretty_map_fun(plot_area = cfp,
                           occ_data = df_clean_sp,
                           x = 'x_albers',
                           y = 'y_albers',
                           epsg_code = 3310,
                           fill_att = "pr_ab",
                           title = paste(species_name),
                           subtitle = paste0(
                             nrow(df_clean_sp %>% filter(pr_ab == 1)),
                             " presences; ",
                             nrow(df_clean_sp %>% filter(pr_ab == 0)),
                             " absences"
                           ))
  
  raw_maps <-
    levelplot(
      all_raw,
      main = paste0(species_name, ": Current distribution"),
      par.settings = myTheme,
      layout=c(4, 2)
    ) +
    layer(sp.polygons(cfp.pol, fill = 'transparent', col = 1))
  
  
  threshold_maps <-
    levelplot(
      sens_spec,
      main = paste0(species_name, ": Current distribution with sens=spec threshold"),
      par.settings = myTheme,
      layout=c(4, 2)
    ) +
    layer(sp.polygons(cfp.pol, fill = 'transparent', col = 1))
  
  pdf(file = file.path(sp_dir, 'sdm_outputs.pdf'))
  print(pa_map)
  print(auc_plot)
  print(raw_maps)
  print(threshold_maps)
  print(varImpPlot(rf_final))
  dev.off()
  print(species_name)
}