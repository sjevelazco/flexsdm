## %######################################################%##
#                                                          #
####                Auxiliary functions                 ####
#                                                          #
## %######################################################%##


pre_tr_te <- function(data, p_names, h) {
  train <- list()
  test <- list()

  if (any(c("train", "train-test", "test")
  %in%
    unique(data[, p_names[h]]))) {
    np2 <- 1

    filt <- grepl("train", data[, p_names[h]])
    train[[1]] <- data[filt, ] %>%
      dplyr::select(-p_names[!p_names == p_names[h]])

    filt <- grepl("test", data[, p_names[h]])
    test[[1]] <- data[filt, ] %>%
      dplyr::select(-p_names[!p_names == p_names[h]])
  } else {
    np2 <- max(data[p_names[h]])

    for (i in 1:np2) {
      train[[i]] <- data[data[p_names[h]] == i, ] %>%
        dplyr::select(-p_names[!p_names == p_names[h]])

      test[[i]] <- data[data[p_names[h]] != i, ] %>%
        dplyr::select(-p_names[!p_names == p_names[h]])
    }
  }
  return(list(train = train, test = test, np2 = np2))
}


# Inverse bioclim
inv_bio <- function(e, p) {
  Model <- dismo::bioclim(e, p)
  r <- dismo::predict(Model, e)
  r <- (r - raster::minValue(r)) /
    (raster::maxValue(r) - raster::minValue(r))
  r <- (1 - r) >= 0.99 # environmental constrain
  r[which(r[, ] == FALSE)] <- NA
  return(r)
}


# Inverse geo
inv_geo <- function(e, p, d) {
  r <- raster::rasterize(p, e)
  b <- raster::buffer(r, width = d)
  e[!is.na(b[])] <- NA
  return(e)
}

# Boyce
# This function calculate Boyce index performance metric. Codes were adapted from enmSdm package. Boyce have
# value between -1 and +1, with a value tending toward +1 indicating good to perfect predictions, values
# around 0 indicating predictions no different from those obtained by chance,
# and values toward -1 indicating counter-predictions.
boyce <- function(pres,
                  contrast,
                  n_bins = 101,
                  n_width = 0.1) {
  lowest <- min(c(pres, contrast), na.rm = TRUE)
  highest <- max(c(pres, contrast), na.rm = TRUE) + .Machine$double.eps
  window_width <- n_width * (highest - lowest)

  lows <- seq(lowest, highest - window_width, length.out = n_bins)
  highs <- seq(lowest + window_width + .Machine$double.eps, highest, length.out = n_bins)

  ## initiate variables to store predicted/expected (P/E) values
  freq_pres <- NA
  freq_contrast <- NA

  # tally proportion of test presences/background in each class
  for (i in 1:n_bins) {
    # number of presence predictions in a class
    freq_pres[i] <-
      sum(pres >= lows[i] & pres < highs[i], na.rm = TRUE)

    # number of background predictions in this class
    freq_contrast[i] <-
      sum(contrast >= lows[i] & contrast < highs[i], na.rm = TRUE)
  }

  # mean bin prediction
  mean_pred <- rowMeans(cbind(lows, highs))

  # add small number to each bin that has 0 background frequency but does have a presence frequency > 0
  if (any(freq_pres > 0 & freq_contrast == 0)) {
    small_value <- 0.5
    freq_contrast[freq_pres > 0 & freq_contrast == 0] <- small_value
  }

  # remove classes with 0 presence frequency
  if (any(freq_pres == 0)) {
    zeros <- which(freq_pres == 0)
    mean_pred[zeros] <- NA
    freq_pres[zeros] <- NA
    freq_contrast[zeros] <- NA
  }

  # remove classes with 0 background frequency
  if (any(0 %in% freq_contrast)) {
    zeros <- which(freqPres == 0)
    mean_pred[zeros] <- NA
    freq_pres[zeros] <- NA
    freq_contrast[zeros] <- NA
  }

  P <- freq_pres / length(pres)
  E <- freq_contrast / length(contrast)
  PE <- P / E

  # remove NAs
  rm_nas <- stats::complete.cases(data.frame(mean_pred, PE))
  mean_pred <- mean_pred[rm_nas]
  PE <- PE[rm_nas]

  # calculate Boyce index
  result <- stats::cor(x = mean_pred, y = PE, method = "spearman")
  return(result)
}
