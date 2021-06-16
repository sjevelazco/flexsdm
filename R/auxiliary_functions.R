## %######################################################%##
#                                                          #
####                Auxiliary functions                 ####
#                                                          #
## %######################################################%##


#' pre_tr_te
#'
#' @noRd
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
      train[[i]] <- data[data[p_names[h]] != i, ] %>%
        dplyr::select(-p_names[!p_names == p_names[h]])

      test[[i]] <- data[data[p_names[h]] == i, ] %>%
        dplyr::select(-p_names[!p_names == p_names[h]])
    }
  }
  return(list(train = train, test = test, np2 = np2))
}



# Inverse bioclim
inv_bio <- function(e, p) {
  e <- raster::stack(e)
  model <- dismo::bioclim(e, as.matrix(p))
  r <- terra::predict(model, e)
  r <- terra::rast(r)
  r <- (r - terra::minmax(r)[1]) /
    (terra::minmax(r)[2] - terra::minmax(r)[1])
  r <- (1 - r) >= 0.99 # environmental constrain
  r[which(r[, ] == FALSE)] <- NA
  return(r)
}


#' Inverse geo
#'
#' @noRd
#'
inv_geo <- function(e, p, d) {
  colnames(p) <- c("x", "y")
  p <- terra::vect(p, geom = c("x", "y"))
  r <- terra::rasterize(p, e)
  b <- terra::buffer(r, width = d)
  e <- mask(e, b, maskvalues = 1)
  return(e)
}

#' Boyce
#'
#' @description This function calculate Boyce index performance metric. Codes were adapted from
#' enmSdm package.
#'
#' @noRd
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
    zeros <- which(freq_pres == 0)
    mean_pred[zeros] <- NA
    freq_pres[zeros] <- NA
    freq_contrast[zeros] <- NA
  }

  P <- freq_pres / length(pres)
  E <- freq_contrast / length(contrast)
  PE <- P / E

  # remove NAs
  rm_nas <- stats::complete.cases(data.frame(mean_pred, PE))
  # mean_pred <- mean_pred[rm_nas]
  # PE <- PE[rm_nas]

  # calculate Boyce index
  result <- stats::cor(
    x = ifelse(is.na(mean_pred), 0, mean_pred),
    y = ifelse(is.na(PE), 0, PE), method = "spearman"
  )
  return(result)
}


#' Predict maxnet
#'
#' @noRd
predict_maxnet <- function(object, newdata, clamp = TRUE, type = c("link", "exponential", "cloglog", "logistic"), ...) {
  if (clamp) {
    for (v in intersect(names(object$varmax), names(newdata))) {
      newdata[, v] <- pmin(
        pmax(newdata[, v], object$varmin[v]),
        object$varmax[v]
      )
    }
  }
  terms <- sub(
    "hinge\\((.*)\\):(.*):(.*)$", "hingeval(\\1,\\2,\\3)",
    names(object$betas)
  )
  terms <- sub(
    "categorical\\((.*)\\):(.*)$", "categoricalval(\\1,\\2)",
    terms
  )
  terms <- sub(
    "thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)",
    terms
  )
  f <- formula(paste(
    "~", paste(terms, collapse = " + "),
    "-1"
  ))
  hingeval <- function(x, min, max) {
    pmin(1, pmax(0, (x - min) / (max - min)))
  }
  mm <- stats::model.matrix(f, data.frame(newdata))
  if (clamp) {
    mm <- t(pmin(
      pmax(t(mm), object$featuremins[names(object$betas)]),
      object$featuremaxs[names(object$betas)]
    ))
  }
  link <- (mm %*% object$betas) + object$alpha
  type <- match.arg(type)
  if (type == "link") {
    return(link)
  }
  if (type == "exponential") {
    return(exp(link))
  }
  if (type == "cloglog") {
    return(1 - exp(0 - exp(object$entropy + link)))
  }
  if (type == "logistic") {
    return(1 / (1 + exp(-object$entropy - link)))
  }
}

#' Outliers with Reverse Jackknife
#'
#' @noRd
#'
rev_jack <- function(v) {
  v2 <- v
  v <- unique(v)
  lgh <- length(v) - 1
  t1 <- (0.95 * sqrt(length(v))) + 0.2
  x <- sort(v)
  y <- rep(0, lgh)
  for (i in seq_len(lgh)) {
    x1 <- x[i + 1]
    if (x[i] < mean(v)) {
      y[i] <- (x1 - x[i]) * (mean(v) - x[i])
    }
    else {
      y[i] <- (x1 - x[i]) * (x1 - mean(v))
    }
  }
  my <- mean(y)
  z <- y / (sqrt(sum((y - my)^2) / lgh))
  out <- rep(0, length(v2))
  if (any(z > t1)) {
    f <- which(z > t1)
    v <- x[f]
    if (v < median(x)) {
      xa <- (v2 <= v) * 1
      out <- out + xa
    }
    if (v > median(x)) {
      xb <- (v2 >= v) * 1
      out <- out + xb
    }
  }
  else {
    out <- out
  }
  return(which(out == 1))
}
