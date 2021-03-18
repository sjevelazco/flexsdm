##%######################################################%##
#                                                          #
####                SDM Function                        ####
#                                                          #
##%######################################################%##
# Written by Brooke Rose
#' Title
#'
#' @param df 
#' @param eval 
#' @param calib 
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
                 eval, # portion of df for evaluation
                 calib, # portion of df for calibration,
                 pr_ab, # response variable (0/1)
                 env_preds, # list of characters for the environmental predictors 
                 sp_area, # polygon for delineating prediction area
                 pred_rasters, # prediction rasters
                 species_name, # character of species name
                 dir_save,
                 cores=1) # directory for saving model objects and spatial predictions
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
  devtools::install_github("andrefaa/ENMTML") 
  require(sdm)
  require('ENMTML')
  
  if(!dir.exists(paste0(dir_save, 'models/', sep = ''))) {
    dir.create(paste0(dir_save, 'models/', sep = ''))
  }
  
  # data frames with only response and predictor variables
  df_clean <- df %>% dplyr::select(all_of(env_preds), pr_ab)
  eval <- eval %>% dplyr::select(all_of(env_preds), pr_ab)
  calib <- calib %>% dplyr::select(all_of(env_preds), pr_ab)
  
  # preparing prediction rasters
  pred_crop <- raster::crop(pred_rasters, sp_area)
  pred_mask <- raster::mask(pred_crop, sp_area)
  
  # spatial data frame
  spatial_df <- st_as_sf(
    df,
    coords = c('x_albers', 'y_albers'),
    crs = crs(pred_mask),
    remove = FALSE
  )
  
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
  
  
  
##%######################################################%##
#                                                          #
####             Generalized Linear Models              ####
#                                                          #
##%######################################################%##
  #### Final model built with all data
  glm.formula <-
    makeFormula("pr_ab", calib[, env_preds],
                "quadratic", interaction.level = 1)
  
  glm_final <-
      glmStart <- glm(glm.formula,
                      data = df_clean,
                      family = binomial)
  
  saveRDS(glm_final, file = file.path(dir_save, 'models/glm_final.rda'))
  
  # model prediction on evaluation data
  full_pred_glm <- predict(glm_final, df_clean, type = "response")
  
  # model evaluation on model built using all data
  full_eval_glm <- sdm::evaluates(x = df_clean$pr_ab, p = full_pred_glm)
  
  # variable importance
  varImp_full_glm <- varImp(glm_final, scale = FALSE)
  
  
  
  #### Training and testing model
  glm_train <- glm(glm.formula,
                    data = calib,
                    family = binomial)
  
  # saveRDS(glm_train, file = file.path(dir_save, 'models/glm_train.rda'))
  
  # model prediction on evaluation data
  test_pred_glm <- predict(glm_train, eval, type = "response")
  
  # model evaluation on test (eval) data
  test_eval_glm <- sdm::evaluates(x = eval$pr_ab, p = test_pred_glm)
  
  # variable importance
  #### I think it is not necessary calculate variable importance for a model used only for testing. Calculate variable 
  #### importance with the model constructed with all data is enough :^)
  # varImp_glm <- caret::varImp(glm_train, scale = FALSE)
  

##%######################################################%##
#                                                          #
####            Generalized Additive Models             ####
#                                                          #
##%######################################################%##
  #### Final model built with all data
  # gam.formula <- paste("s(", env_preds, ")",sep="") #",k=3)",
  # WE NEED TO FIND A WAY TO CREAT THE FORMULA AUTOMATICALLY without writing variable names
  
  gam_final <- gam(
    pr_ab ~ s(cwd) + s(aet) + s(tmin) + s(ppt_djf) + s(ppt_jja) +
      s(pH) + s(awc) + s(depth) + s(percent_clay) + landform,
    data = df_clean,
    family = "binomial"
  )
  
  saveRDS(gam_final, file = file.path(dir_save, 'models/gam_final.rda'))
  
  # model prediction on all data
  full_pred_gam <- predict(gam_final, df_clean, type = "response")
  
  # model evaluation on model built using all data
  full_eval_gam <- sdm::evaluates(x = df_clean$pr_ab, p = full_pred_gam)
  
  # variable importance
  varImp_full_gam <- caret::varImp(gam_final, scale = FALSE)
  
  
  #### Training and testing model
  gam_train <- gam::gam(
    pr_ab ~ s(cwd) + s(aet) + s(tmin) + s(ppt_djf) + s(ppt_jja) +
      s(pH) + s(awc) + s(depth) + s(percent_clay) + landform,
    data = calib,
    family = "binomial"
  )
  
  # saveRDS(gam_train, file = file.path(dir_save, 'models/gam_train.rda'))
  
  # model prediction on evaluation data
  test_pred_gam <- predict(gam_train, eval, type = "response")
  
  # model evaluation on test (eval) data
  test_eval_gam <- sdm::evaluates(x = eval$pr_ab, p = test_pred_gam)
  
  # variable importance
  # varImp_eval_gam <- varImp(gam_train, scale = FALSE)
  
  
##%######################################################%##
#                                                          #
####                   Random Forest                    ####
#                                                          #
##%######################################################%##
  
  #### Final model built with all data
  # Find best parameters for final model
  
  tune_grid <- expand.grid(mtry = seq(2, length(env_preds), 1))
  
  set.seed(123)
  betst_tune <- caret::train(
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
  # ggplot(betst_tune)
  
  # model built with all data
  set.seed(123)
  rf_final <- randomForest::randomForest(
    pr_ab ~ .,
    data = dplyr::mutate(df_clean, pr_ab = as.factor(pr_ab)),
    ntree = 1000,
    importance = TRUE,
    type = 'classification',
    mtry = betst_tune$bestTune$mtry
  )
  
  saveRDS(rf_final, file = file.path(dir_save, 'models/rf_final.rda'))
  
  # model prediction on all data
  full_pred_rf <- predict(rf_final, df_clean, type = "prob")[, 2]
  
  # model evaluation for model built using evaluation data
  full_eval_rf <- sdm::evaluates(x = df_clean$pr_ab, p = full_pred_rf)
  
  # variable importance
  varImp_full_rf <- caret::varImp(rf_final, scale = FALSE)
  # varImpPlot(rf_final)
  varImp_full_rf <- randomForest::importance(rf_final)
  
  
  #### Training and testing model
  # with best tuning 
  set.seed(123)
  rf_train <- randomForest::randomForest(
    x = calib[, env_preds],
    y = as.factor(calib$pr_ab),
    ntree = 1000,
    importance = TRUE, 
    type='classification',
    mtry=betst_tune$bestTune$mtry
  )
  
  # saveRDS(rf_train, file = file.path(dir_save, 'models/rf_train.rda'))
  
  # model prediction on evaluation data
  test_pred_rf <- predict(rf_train, eval, type = "prob")[, 2]
  
  # model evaluation for model built using evaluation data
  test_eval_rf <- sdm::evaluates(x = eval$pr_ab, p = test_pred_rf)
  
  # variable importance
  # varImp_eval_rf <- varImp(rf_train, scale = FALSE)
  
  
##%######################################################%##
#                                                          #
####              Boosted Regression Tree               ####
#                                                          #
##%######################################################%##
  #### Final model built with all data
  # Find best parameters for final model
  tune_grid <- expand.grid(
    interaction.depth = seq(2, 10, 20),
    n.trees = c(100, 200, 500, 1000, 2000, 3000),
    shrinkage = 0.1,
    n.minobsinnode = 20
  )
  
  set.seed(123)
  betst_tune <- caret::train(
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
  
  # ggplot(betst_tune)
  
  set.seed(123)
  brt_final <-
    gbm::gbm(
      pr_ab ~ .,
      data = df_clean,
      distribution = "bernoulli",
      n.trees = betst_tune$bestTune$n.trees,
      interaction.depth = betst_tune$bestTune$interaction.depth,
      shrinkage = betst_tune$bestTune$shrinkage,
      n.minobsinnode = betst_tune$bestTune$n.minobsinnode
    )
  
  saveRDS(brt_final, file = file.path(dir_save, 'models/brt_final.rda'))
  
  # brt_final_tune = gbm.perf(brt_final, method = "cv", plot.it = F)
  
  # model prediction on all data
  full_pred_brt <- predict(brt_final, df_clean, type = "response", n.trees = brt_final_tune)
  
  # model evaluation for model built using all data
  full_eval_brt <- sdm::evaluates(x = df_clean$pr_ab, p = full_pred_brt)
  
  # variable importance
  varImp_full_brt <- caret::varImp(brt_final, scale = FALSE, 
                            numTrees = betst_tune$bestTune$n.trees)
  
  #### Training and testing model
  # with best tuning
  set.seed(123)
  brt_train <-
    gbm::gbm(pr_ab ~ .,
        data = calib,
        distribution = "bernoulli",
        n.trees = betst_tune$bestTune$n.trees,
        interaction.depth = betst_tune$bestTune$interaction.depth,
        shrinkage = betst_tune$bestTune$shrinkage,
        n.minobsinnode = betst_tune$bestTune$n.minobsinnode
    )
  
  
  # saveRDS(brt_train, file = file.path(dir_save, 'models/brt_train.rda'))
  
  # number of trees to use in predictions
  # brt_train_tune = gbm.perf(brt_train, method = "cv", plot.it = F)
  
  # model prediction on evaluation data
  test_pred_brt <- predict(brt_train, eval, type = "response")
  
  # model evaluation for model built using evaluation data
  test_eval_brt <- sdm::evaluates(x = eval$pr_ab, p = test_pred_brt)
  
  # variable importance
  # varImp_test_brt <- varImp(brt_train, scale = FALSE, numTrees = brt_train_tune)
  
  
##%######################################################%##
#                                                          #
####              Support Vector Machines               ####
#                                                          #
##%######################################################%##
  
  #### Final model built with all data
  # Find best parameters for final model
  
  # Tune final model
  # svm_tune_final <- e1071::tune(
  #   svm,
  #   pr_ab ~ .,
  #   data = df_clean,
  #   ranges = list(gamma = 2 ^ (-1:1), cost = 2 ^ (2:4)),
  #   tunecontrol = tune.control(sampling = "fix")
  # )
  
  # Combination of parameters values to be tested
  tune_grid <-
    expand.grid(C = c(1, 2, 4, 8, 16),
                sigma = c(0.001, 0.01, 0.1, 0.2))
  
  set.seed(123)
  betst_tune <- caret::train(
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
  
  # ggplot(betst_tune)
  
  set.seed(123)
  svm_final <- kernlab::ksvm(
    pr_ab ~ .,
    data = df_clean,
    type = "C-svc",
    kernel = "rbfdot",
    kpar = list(sigma = betst_tune$bestTune$sigma),
    C = betst_tune$bestTune$C,
    prob.model = TRUE
  )
  
  saveRDS(svm_final, file = file.path(dir_save, 'models/svm_final.rda'))
  
  
  # model prediction on all data
  full_pred_svm <- predict(svm_final, df_clean, type = 'prob')[, 2]
  
  # model evaluation for model built using all data
  full_eval_svm <- sdm::evaluates(x = df_clean$pr_ab, p = full_pred_svm)
  
  # variable importance
  varImp_full_svm <- caret::filterVarImp(df_clean[, env_preds], full_pred_svm, nonpara = FALSE)
  # varImp_full_svm <- varImp_full_svm / sum(varImp_full_svm)
  
  #### Training and testing model
  # with best tuning
  # svm_tune_train <- tune(
  #   svm,
  #   pr_ab ~ .,
  #   data = calib,
  #   ranges = list(gamma = 2 ^ (-1:1), cost = 2 ^ (2:4)),
  #   tunecontrol = tune.control(sampling = "fix")
  # )
  
  set.seed(123)
  svm_train <- kernlab::ksvm(
    pr_ab ~ .,
    data = calib,
    type = "C-svc",
    kernel = "rbfdot",
    kpar = list(sigma = betst_tune$bestTune$sigma),
    C = betst_tune$bestTune$C,
    prob.model = TRUE
  )
  
  saveRDS(svm_train, file = file.path(dir_save, 'models/svm_train.rda'))
  
  # model prediction on evaluation data
  test_pred_svm <- predict(svm_train, eval, type = 'prob')[, 2]
  
  # model evaluation for model built using evaluation data
  test_eval_svm <- sdm::evaluates(x = eval$pr_ab, p = test_pred_svm)
  
  
  
##%######################################################%##
#                                                          #
####            Artificial Neural Networks              ####
#                                                          #
##%######################################################%##
  # cv_nnet_train <- biomod2:::.CV.nnet(Input = calib[, env_preds],
  #                                     Target = calib$pr_ab)
  # Tune final model
  
  #### Final model built with all data
  # Find best parameters for final model
  # Combination of parameters values to be tested
  tune_grid <-
    expand.grid(size = c(2, 4, 6, 8, 10), 
                decay = c(0.001, 0.01, 0.05, 0.1))
  
  cl <- makePSOCKcluster(cores)
  registerDoParallel(cl)
  
  set.seed(123)
  betst_tune <- caret::train(
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
  
  ggplot(betst_tune)
  
  # final model
  # cv_nnet_final <- biomod2:::.CV.nnet(Input = df_clean[, env_preds],
                                      # Target = df_clean$pr_ab)
  set.seed(123)
  nnet_final <- nnet::nnet(
    pr_ab ~ .,
    data = df_clean,
    size = betst_tune$bestTune$size,
    rang = 0.1,
    decay = betst_tune$bestTune$decay,
    maxit = 200,
    trace = FALSE
  )
  
  saveRDS(nnet_final, file = file.path(dir_save, 'models/nnet_final.rda'))
  
  # model predictions on all data
  full_pred_nnet <- predict(nnet_final, df_clean[, env_preds])
  
  # model evaluation for model built on evaluation data
  full_eval_nnet <- sdm::evaluates(x = df_clean$pr_ab, p = full_pred_nnet)
  
  # variable importance
  varImp_full_nnet <- caret::varImp(nnet_final)
  
  #### Training and testing model
  # with best tuning
  nnet_train <- nnet::nnet(
    pr_ab ~ .,
    data = calib,
    size = betst_tune$bestTune$size,
    rang = 0.1,
    decay = betst_tune$bestTune$decay,
    maxit = 200,
    trace = FALSE
  )
  
  # saveRDS(nnet_train, file = file.path(dir_save, 'models/nnet_train.rda'))
  
  # model predictions on evaluation data
  test_pred_nnet <- predict(nnet_train, eval[, env_preds])
  
  # model evaluation for model built on evaluation data
  test_eval_nnet <- sdm::evaluates(x = eval$pr_ab, p = test_pred_nnet)
  
  
  # STOP CLUSTER
  stopCluster(cl)
  
##%######################################################%##
#                                                          #
####                  post-processing                   ####
#                                                          #
##%######################################################%##

  
  ### AUC based on evaluation model
  auc <- list(
    test_eval_glm@statistics$AUC,
    test_eval_gam@statistics$AUC,
    test_eval_rf@statistics$AUC,
    test_eval_brt@statistics$AUC,
    test_eval_svm@statistics$AUC,
    test_eval_nnet@statistics$AUC
  )
  names(auc) <- c('GLM', 'GAM', 'RF', 'BRT', 'SVM', 'NNET')
  
  
  saveRDS(auc, file = file.path(dir_save, 'models/train_auc.rda'))
  
  auc <- readRDS(file = file.path(dir_save, 'models/train_auc.rda'))
  
  
  ### spatial predictions
  beginCluster()
  
  raw_preds <- raster::brick(
    clusterR(pred_mask, raster::predict, args=list(model = glm_final, type = "response")),
    clusterR(pred_mask, raster::predict, args=list(model = gam_final, type = "response")),
    clusterR(pred_mask, raster::predict, args=list(model = rf_final, type = "prob", index = 2)),
    clusterR(pred_mask, raster::predict, args=list(model = brt_final, n.trees = brt_final_tune, type = "response")),
    clusterR(pred_mask, raster::predict, args=list(model = svm_final, type = "prob", index = 2)),
    clusterR(pred_mask, raster::predict, args=list(model = nnet_final)))
  
  names(raw_preds) <- c('glm', 'gam', 'rf', 'brt', 'svm', 'nnet')
  
  ### Ensembles
  
  ensemble <- 
    raster::brick(
      clusterR(raw_preds, fun = calc, args = list(fun = mean)),
      clusterR(raw_preds, fun = weighted.mean, args = list(w=unlist(auc))),
      clusterR(raw_preds, fun = calc, args = list(fun = sd)))
  
  names(ensemble) <- c('mean', 'weighted average', 'standard deviation')
  
  endCluster()
  
  ensemble_extract <- raster::extract(ensemble, spatial_df, df = TRUE, sp = TRUE)
  
  ensemble_data <- st_as_sf(ensemble_extract) %>%
    dplyr::select(pr_ab, mean, weighted.average, standard.deviation) %>%
    na.omit()
  
  mean_eval <- evaluates(x = ensemble_data$pr_ab, p = ensemble_data$mean)
  w_avg_eval <- evaluates(x = ensemble_data$pr_ab, p = ensemble_data$weighted.average)
  
  ### auc for models built on all data
  auc_final <- list(
    full_eval_glm@statistics$AUC,
    full_eval_gam@statistics$AUC,
    full_eval_rf@statistics$AUC,
    full_eval_brt@statistics$AUC,
    full_eval_svm@statistics$AUC,
    full_eval_nnet@statistics$AUC,
    mean_eval@statistics$AUC,
    w_avg_eval@statistics$AUC
  )
  
  names(auc_final) <- c('GLM',
                        'GAM',
                        'RF',
                        'BRT',
                        'SVM',
                        'NNET',
                        'ensemble_mean',
                        'ensemble_weighted_avg')
  
  saveRDS(auc_final, file = file.path(dir_save, 'models/final_auc.rda'))

  
  ## all rasters
  all_raw <- raster::stack(raw_preds, ensemble)
  
  writeRaster(all_raw, filename = file.path(dir_save, 'models/raw_pred_maps.grd'), overwrite = TRUE)
  
  ### Thresholds for based on final models
  
  # list of threshold data frames
  threshold <- list(
    test_eval_glm@threshold_based,
    test_eval_gam@threshold_based,
    test_eval_rf@threshold_based,
    test_eval_brt@threshold_based,
    test_eval_svm@threshold_based,
    test_eval_nnet@threshold_based,
    mean_eval@threshold_based,
    w_avg_eval@threshold_based
  )
  
  names(threshold) <- c('GLM',
                        'GAM',
                        'RF',
                        'BRT',
                        'SVM',
                        'NNET',
                        'ensemble_mean',
                        'ensemble_weighted_avg')
  
  sapply(names(threshold),
         function (x)
           utils::write.table(
             threshold[[x]],
             file.path(dir_save, 'models/', paste(x,"_thresholds.txt")),
             sep = "\t",
             row.names = F
           ))
  
  sens_spec <- list() # initiating list for sensitivity = specificity threshold maps
  
  # function that creates rasters with everything less than threshold = 0 and retain values above 0
  for (i in 1:8) {
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
                        'NNET',
                        'ensemble_mean',
                        'ensemble_weighted_avg')
  
  writeRaster(
    sens_spec,
    filename = file.path(dir_save, 'models/sens_spec_maps.grd'),
    overwrite = TRUE)
    
    
    ##### PDF Output #####
    auc_final_df <- as.data.frame(auc_final) %>%
      dplyr::rename(mean = ensemble_mean,
                    w_average = ensemble_weighted_avg) %>%
      pivot_longer(cols = 1:8,
                   names_to = "Model",
                   values_to = "AUC")  %>%
      mutate(model = 'final')
    
    auc_eval_df <- as.data.frame(auc) %>%
      dplyr::mutate(mean = NA,
                    w_average = NA) %>%
      pivot_longer(cols = 1:8,
                   names_to = "Model",
                   values_to = "AUC") %>%
      mutate(model = 'eval')
    
    auc_df <- bind_rows(auc_final_df, auc_eval_df)
    
    auc_plot <-
      ggplot(auc_df, aes(
        x = Model,
        y = AUC,
        label = round(AUC, 3),
        fill = model
      )) +
      geom_point(size = 3) + geom_label(size = 5) +
      labs(
        title = paste0(species_name, ": AUC Comparison"),
        x = "Model Type",
        y = "AUC"
      ) +
      theme(text = element_text(size = 15, family = "serif"),
            axis.title.x = element_text(vjust = .25))
    
    presence <- spatial_df %>%
      filter(pr_ab == 1)
    
    absence <- spatial_df %>%
      filter(pr_ab == 0)
    
    p.points <- as(presence, 'Spatial')
    a.points <- as(absence, 'Spatial')
    cfp.pol <- as(sp_area, 'Spatial')
    
    myTheme <- rasterTheme(region = rev(terrain.colors(7)))
    
    
    raw_maps <-
      levelplot(
        all_raw,
        main = paste0(species_name, ": Current distribution"),
        par.settings = myTheme,
        layout=c(3, 3)
      ) +
      layer(sp.polygons(cfp.pol, fill = 'transparent', col = 1))
    
    
    threshold_maps <-
      levelplot(
        sens_spec,
        main = paste0(species_name, ": Current distribution with sens=spec threshold"),
        par.settings = myTheme,
        layout=c(3, 3)
      ) +
      layer(sp.polygons(cfp.pol, fill = 'transparent', col = 1))
    
    pdf(
      file = paste0(dir_save, 'models/',
        species_name,
        '_sdm_outputs.pdf'
      )
    )
    print(auc_plot)
    print(raw_maps)
    print(threshold_maps)
    print(varImpPlot(rf_final))
    dev.off()
    
}
