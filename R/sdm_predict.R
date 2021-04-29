sdm_predict <- function(models, pred, pred_proj, thr, calibarea = NULL, clamp = TRUE, pred_type = "cloglog") {

  #### Prepare datasets ####
  # Prepare model list
  m <- lapply(models, function(x) x[[1]])
  names(m) <- paste0("m_", 1:length(m))

  # Extract model names object
  clss <- sapply(m, function(x) class(x)[1]) %>%
    tolower() %>%
    gsub(".formula", "", .)

  # Create list for storing raster for current condition
  model_c <- as.list(names(m))
  names(model_c) <- names(m)

  # Transform raster to data.frame
  pred_df <- data.frame(
    ncell = 1:raster::ncell(pred),
    raster::as.data.frame(pred)
  ) %>%
    stats::na.exclude(pred_df)
  rownames(pred_df) <- pred_df$ncell
  pred_df$ncell <- NULL
  factvar <- names(pred_df)[!names(pred_df) %in% names(pred)]
  if (length(factvar)) {
    factvar2 <- gsub("_category", "", factvar)
    for (i in 1:length(factvar2)) {
      df1 <- pred[[factvar2[i]]]@data@attributes[[1]]
      cc <- c("category")
      names(cc) <- factvar
      pred_df[factvar] <-
        dplyr::left_join(pred_df[factvar], df1, by = cc)[, 2]
      names(pred_df)[names(pred_df) == factvar] <- factvar2
      pred_df[, factvar2] <- pred_df[, factvar2] %>% as.factor()
    }
  }

  #### graf models ####
  wm <- which(clss == "graf")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      model_c[[wm]] <- GRaF::predict.graf(m[[wm]], pred, type = "response", CI = NULL)
    }
  }


  #### gam models ####
  wm <- which(clss == "gam")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      r <- pred[[1]]
      values(r) <- NA

      # Test factor levels
      f <- which(sapply(m[[wm]]$data, class) == "factor")
      if (f > 0) {
        for (ii in 1:length(f)) {
          vf <- m[[wm]]$data[, f[ii]] %>% unique()
          vf2 <- pred_df[, names(f[ii])] %>% unique()
          vfilter <- list()
          if (sum(!vf2 %in% vf) > 0) {
            vfilter[[ii]] <- !pred_df[, names(f[ii])] %in% vf
          }
        }
        if (length(vfilter) > 1) {
          vfilter <- vapply(do.call("rbind", vfilter), any, logical(1))
        } else {
          vfilter <- vfilter[[1]]
        }
      }

      if (sum(vfilter) > 0) {
        v <- rep(0, nrow(pred_df))
        v[!vfilter] <-
          gam::predict.Gam(m[[wm]], pred_df[!vfilter, ], type = "response")
        r[as.numeric(rownames(pred_df))] <- v
        rm(v)
      } else {
        r[as.numeric(rownames(pred_df))] <-
          gam::predict.Gam(m[[wm]], pred_df, type = "response")
      }

      model_c[[wm]] <- r
    }
  }

  #### glm models ####
  wm <- which(clss == "glm")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      r <- pred[[1]]
      values(r) <- NA

      # Test factor levels
      f <- which(sapply(m[[wm]]$data, class) == "factor")
      if (f > 0) {
        for (ii in 1:length(f)) {
          vf <- m[[wm]]$data[, f[ii]] %>% unique()
          vf2 <- pred_df[, names(f[ii])] %>% unique()
          vfilter <- list()
          if (sum(!vf2 %in% vf) > 0) {
            vfilter[[ii]] <- !pred_df[, names(f[ii])] %in% vf
          }
        }
        if (length(vfilter) > 1) {
          vfilter <- vapply(do.call("rbind", vfilter), any, logical(1))
        } else {
          vfilter <- vfilter[[1]]
        }
      }

      if (sum(vfilter) > 0) {
        v <- rep(0, nrow(pred_df))
        v[!vfilter] <-
          stats::predict(m[[wm]], pred_df[!vfilter, ], type = "response")
        r[as.numeric(rownames(pred_df))] <- v
        rm(v)
      } else {
        r[as.numeric(rownames(pred_df))] <-
          stats::predict(m[[wm]], pred_df, type = "response")
      }

      model_c[[wm]] <- r
    }
  }

  #### gbm models ####
  wm <- which(clss == "gbm")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      r <- pred[[1]]
      values(r) <- NA
      r[as.numeric(rownames(pred_df))] <-
        stats::predict(m[[wm]], pred_df, type = "response")

      model_c[[wm]] <- r
    }
  }


  #### maxnet models ####
  wm <- which(clss == "maxnet")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      model_c[[wm]] <- raster::predict(pred, m[[wm]], clamp = clamp, type = pred_type)
    }
  }

  #### nnet class ####
  wm <- which(clss == "nnet")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      r <- pred[[1]]
      values(r) <- NA

      # Test factor levels
      f <- (m[[wm]]$xlevels)

      if (length(f) > 0) {
        for (ii in 1:length(f)) {
          vf <- f[[ii]] %>%
            unique() %>%
            as.numeric()
          vf2 <- pred_df[, names(f[ii])] %>% unique()
          vfilter <- list()
          if (sum(!vf2 %in% vf) > 0) {
            vfilter[[ii]] <- !pred_df[, names(f[ii])] %in% vf
          }
        }
        if (length(vfilter) > 1) {
          vfilter <- vapply(do.call("rbind", vfilter), any, logical(1))
        } else {
          vfilter <- vfilter[[1]]
        }
      }

      if (sum(vfilter) > 0) {
        v <- rep(0, nrow(pred_df))
        v[!vfilter] <-
          stats::predict(m[[wm]], pred_df[!vfilter, ], type = "raw")
        r[as.numeric(rownames(pred_df))] <- v
        rm(v)
      } else {
        r[as.numeric(rownames(pred_df))] <-
          stats::predict(m[[wm]], pred_df, type = "raw")
      }

      model_c[[wm]] <- r
    }
  }


  #### randomforest class ####
  wm <- which(clss == "randomforest")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      r <- pred[[1]]
      values(r) <- NA

      # Test factor levels
      f <-
        m[[wm]]$forest$xlevels[sapply(m[[wm]]$forest$xlevels, function(x) {
          class(x) == "character"
        })]
      if (length(f) > 0) {
        for (ii in 1:length(f)) {
          vf <- f[[ii]] %>%
            unique() %>%
            as.numeric()
          vf2 <- pred_df[, names(f[ii])] %>% unique()
          vfilter <- list()
          if (sum(!vf2 %in% vf) > 0) {
            vfilter[[ii]] <- !pred_df[, names(f[ii])] %in% vf
          }
        }
        if (length(vfilter) > 1) {
          vfilter <- vapply(do.call("rbind", vfilter), any, logical(1))
        } else {
          vfilter <- vfilter[[1]]
        }
      }


      if (sum(vfilter) > 0) {
        v <- rep(0, nrow(pred_df))
        v[!vfilter] <-
          stats::predict(m[[wm]], pred_df[!vfilter, ] %>%
            dplyr::mutate(across(
              .cols = names(f),
              .fns = ~ droplevels(.)
            )),
          type = "prob"
          )[, 2]
        r[as.numeric(rownames(pred_df))] <- v
        rm(v)
      } else {
        r[as.numeric(rownames(pred_df))] <-
          stats::predict(m[[wm]], pred_df, type = "prob")[, 2]
      }

      model_c[[wm]] <- r
    }
  }

  #### ksvmj class ####
  wm <- which(clss == "ksvm")
  if (length(wm) > 0) {
    wm <- names(wm)
    for (i in wm) {
      r <- pred[[1]]
      values(r) <- NA

      # Test factor levels
      f_n <- which(sapply(pred_df, class) == "factor") %>% names()
      f_n2 <- m[[wm]]@xmatrix[[1]] %>% colnames()
      f <- lapply(f_n, function(x) {
        gsub(x, "", grep(x, f_n2, value = TRUE))
      })
      names(f) <- f_n

      if (length(f) > 0) {
        for (ii in 1:length(f)) {
          vf <- f[[ii]] %>%
            unique() %>%
            as.numeric()
          vf2 <- pred_df[, names(f[ii])] %>% unique()
          vfilter <- list()
          if (sum(!vf2 %in% vf) > 0) {
            vfilter[[ii]] <- !pred_df[, names(f[ii])] %in% vf
          }
        }
        if (length(vfilter) > 1) {
          vfilter <- vapply(do.call("rbind", vfilter), any, logical(1))
        } else {
          vfilter <- vfilter[[1]]
        }
      }

      if (sum(vfilter) > 0) {
        v <- rep(0, nrow(pred_df))
        v[!vfilter] <-
          kernlab::predict(m[[wm]], pred_df[!vfilter, ] %>%
            dplyr::mutate(across(
              .cols = names(f),
              .fns = ~ droplevels(.)
            )), type = "prob")[, 2]
        r[as.numeric(rownames(pred_df))] <- v
        rm(v)
      } else {
        r[as.numeric(rownames(pred_df))] <-
          kernlab::predict(m[[wm]], pred_df, type = "prob")[, 2]
      }
      model_c[[wm]] <- r
    }
  }

  df <- data.frame(alg=c('gam', 'graf', 'glm', 'gbm', 'maxnet', 'nnet', 'randomforest', 'ksvm'),
  names=c('gam', 'gau', 'glm', 'gbm', 'mx', 'nne', 'rf', 'svm'))

  names(model_c) <- left_join(data.frame(alg=clss), df, by = 'alg')[,2]
  model_c <- stack(model_c)


  #### Thresholds ####
  if(!is.null(thr)) {
    if (thr == 'selected_thr') {
      thr_df <- lapply(models, function(x)
        x[["selected_threshold"]])
    }
    if (thr == 'all_thr') {
      thr_df <- lapply(models, function(x)
        x[["threshold_table"]])
    }

    pred_thr <- list()
    for (i in 1:nlayers(model_c)) {
      print(i)
      pred_thr[[i]] <-
        lapply(thr_df[[i]]$values, function(x)
          model_c[[i]] > x) %>% stack
    }
    names(pred_thr) <- names(model_c)

  }

  return(pred)
}
