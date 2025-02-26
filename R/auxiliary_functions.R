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
#'
#' @noRd
#'
bio <- function(data, env_layer) {
  . <- NULL
  if (class(data)[1] != "data.frame") {
    data <- data.frame(data)
  }
  if (!methods::is(env_layer, "SpatRaster")) {
    env_layer <- terra::rast(env_layer)
  }

  data <- na.omit(data)

  result <- env_layer[[1]]
  result[] <- NA

  minv <- apply(data, 2, min)
  maxv <- apply(data, 2, max)
  vnames <- names(data)

  data_2 <- data %>%
    na.omit() %>%
    apply(., 2, sort) %>%
    data.frame()

  rnk <- function(x, y) {
    b <- apply(y, 1, FUN = function(z) sum(x < z))
    t <- apply(y, 1, FUN = function(z) sum(x == z))
    r <- (b + 0.5 * t) / length(x)
    i <- which(r > 0.5)
    r[i] <- 1 - r[i]
    r * 2
  }

  var_df <- terra::as.data.frame(env_layer)
  var_df <- na.omit(var_df)

  k <- (apply(t(var_df) >= minv, 2, all) &
    apply(t(var_df) <= maxv, 2, all))

  for (j in vnames) {
    var_df[k, j] <- rnk(
      data_2[, j],
      var_df[k, j, drop = FALSE]
    )
  }
  var_df[!k, ] <- 0
  res <- apply(var_df, 1, min)
  result[as.numeric(names(res))] <- res
  return(result)
}

inv_bio <- function(e, p) {
  if (!methods::is(e, "SpatRaster")) {
    e <- terra::rast(e)
  }
  r <- bio(data = terra::extract(e, p)[-1], env_layer = e)
  r <- (r - terra::minmax(r)[1]) /
    (terra::minmax(r)[2] - terra::minmax(r)[1])
  r <- r <= 0.01 # environmental constrain
  r[which(r[, ] == FALSE)] <- NA
  return(r)
}


#' Inverse geo
#'
#' @noRd
#'
inv_geo <- function(e, p, d) {
  colnames(p) <- c("x", "y")
  p <- terra::vect(p, geom = c("x", "y"), crs = terra::crs(e))
  b <- terra::buffer(p, width = d)
  b <- terra::rasterize(b, e, background = 0)
  e <- terra::mask(e, b, maskvalues = 1)
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

# maxnet:::predict.maxnet()
#' Predict maxnet
#' @importFrom stats model.matrix
#' @noRd
predict_maxnet <- function(object, newdata, clamp = TRUE, type = c("link", "exponential", "cloglog", "logistic"), ...) {
  categoricalval <- function(x, category) {
    ifelse(x == category, 1, 0)
  }
  thresholdval <- function(x, knot) {
    ifelse(x >= knot, 1, 0)
  }
  hingeval <- function(x, min, max) {
    pmin(1, pmax(0, (x - min) / (max - min)))
  }

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
    "categorical\\((.*)\\):(.*)$", "categoricalval(\\1,\"\\2\")",
    terms
  )
  terms <- sub(
    "thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)",
    terms
  )
  f <- formula(paste("~", paste(terms, collapse = " + "), "-1"))
  mm <- model.matrix(f, data.frame(newdata))
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
    } else {
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
  } else {
    out <- out
  }
  return(which(out == 1))
}

#' Calculate amount of data for each training dataset in a given partition
#'
#' @noRd
#'
n_training <- function(data, partition) {
  . <- partt <- NULL
  if (any(c("train", "train-test", "test")
  %in%
    (data %>%
      dplyr::select(dplyr::starts_with({{ partition }})) %>%
      dplyr::pull() %>%
      unique()))) {
    nn_part <- data %>%
      dplyr::select(dplyr::starts_with({{ partition }})) %>%
      apply(., 2, table) %>%
      data.frame()
    nn_part <- nn_part %>% dplyr::mutate(partt = rownames(nn_part))
    nn_part$partt[grepl("train", nn_part$partt)] <- "train"
    nn_part <- nn_part %>%
      dplyr::filter(partt == "train") %>%
      dplyr::select(-partt)
    nn_part <- colSums(nn_part)
  } else {
    data <- data %>%
      dplyr::select(dplyr::starts_with({{ partition }}))

    nn_part <- list()
    for (ppp in 1:ncol(data)) {
      nn_part[[ppp]] <- data %>%
        dplyr::pull(ppp) %>%
        table() %>%
        c()
      sm <- nn_part[[ppp]] %>% sum()
      nn_part[[ppp]] <- sapply(nn_part[[ppp]], function(x) sum(sm - x))
    }
    nn_part <- unlist(nn_part)
  }
  return(nn_part)
}

#' Calculate number of coefficient for gam models
#'
#' @noRd
#'
n_coefficients <- function(data, predictors, predictors_f = NULL, k = 10) {
  data <- data.frame(data)
  if (k < 0) {
    k <- 10
  }
  if (!is.null(predictors_f)) {
    n_levels <- rep(NA, length(predictors_f))
    for (fff in 1:length(predictors_f)) {
      n_levels[fff] <- unique(data[, predictors_f]) %>%
        na.omit() %>%
        length()
    }
    n_levels <- sum(n_levels)
  } else {
    n_levels <- 0
  }
  n <- (k - 1) * length(predictors) + n_levels
  return(n)
}

#' Euclidean distance for extrapolation
#'
#' @noRd
#'
euc_dist <- function(x, y) {
  if (!methods::is(x, "matrix")) {
    x <- as.matrix(x)
  }
  if (!methods::is(y, "matrix")) {
    y <- as.matrix(y)
  }
  result <- matrix(0, nrow = nrow(x), ncol = nrow(y))
  for (ii in 1:nrow(y)) {
    result[, ii] <- sqrt(colSums((t(x) - y[ii, ])^2))
  }
  rownames(result) <- rownames(x)
  colnames(result) <- rownames(y)
  return(result)
}

#' Moran I, based on ape package
#'
#' @noRd
#'
morani <- function(x, weight, na.rm = FALSE, scaled = TRUE) {
  if (dim(weight)[1] != dim(weight)[2]) {
    stop("'weight' must be a square matrix")
  }
  n <- length(x)
  if (dim(weight)[1] != n) {
    stop("'weight' must have as many rows as observations in 'x'")
  }
  ei <- -1 / (n - 1)
  nas <- is.na(x)
  if (any(nas)) {
    if (na.rm) {
      x <- x[!nas]
      n <- length(x)
      weight <- weight[!nas, !nas]
    } else {
      warning("'x' has missing values: maybe you wanted to set na.rm = TRUE?")
      return(NA)
    }
  }
  rs <- rowSums(weight)
  rs[rs == 0] <- 1
  weight <- weight / rs
  s <- sum(weight)
  m <- mean(x)
  y <- x - m
  cv <- sum(weight * y %o% y)
  v <- sum(y^2)
  res <- (n / s) * (cv / v)
  if (scaled) {
    imax <- (n / s) * (sd(rowSums(weight) * y) / sqrt(v / (n - 1)))
    res <- res / imax
  }

  return(res)
}

#' Mahalanobis distance
#'
#' @noRd
#'
mah_dist <- function(x, y, cov) {
  if (!methods::is(x, "matrix")) {
    x <- as.matrix(x)
  }
  if (!methods::is(y, "matrix")) {
    y <- as.matrix(y)
  }
  result <- matrix(0, nrow = nrow(x), ncol = nrow(y))
  for (ii in 1:nrow(y)) {
    # root square of squared Mahalanobis distance
    result[, ii] <- sqrt(mahalanobis(x = x, center = y[ii, ], cov = cov))
  }
  rownames(result) <- rownames(x)
  colnames(result) <- rownames(y)
  return(result)
}

#' K-means sample
#'
#' @noRd
#'
kf <- function(df, n){
  suppressWarnings(km <- stats::kmeans(df %>% dplyr::select(-c(cell:y)),
                                       centers = n))
  result <- data.frame(cluster = km$cluster, cell = names(km$cluster)) %>%
    dplyr::as_tibble() %>%
    dplyr::arrange(cluster) %>%
    dplyr::group_by(cluster) %>%
    dplyr::sample_n(1, replace = FALSE) %>%
    dplyr::group_by() %>%
    dplyr::pull(cell) %>%
    as.numeric()
  result <- df %>%
    dplyr::filter(cell %in% result) %>%
    dplyr::select(x, y) %>%
    dplyr::as_tibble()
  result$pr_ab <- 0
  return(result)
}
