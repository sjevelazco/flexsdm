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

#' Extract maxent classes based on a formula
#'
#' @param f A formula object.
#' @return A character string concatenating the identified feature classes
#'   (e.g., "lqh"). The classes are returned in a standard order (l, q, p, h, t)
#'   for consistency.
#' @noRd
extract.maxnet.classes <- function(f) {
  # Ensure the input is a formula and get the terms as a character vector
  if (!inherits(f, "formula")) {
    stop("Input must be a formula object.")
  }

  # labels(terms(f)) is the canonical way to get the RHS terms
  formula_terms <- labels(terms(f))

  # Initialize a vector to store the classes we find
  detected_classes <- character(0)

  # --- Logic to detect each class based on term structure ---

  # 1. Check for quadratic features ('q'): terms like I(var^2)
  # The regex looks for a term starting with "I(" and containing "^2)"
  if (any(grepl("^I\\(.*\\^2\\)$", formula_terms))) {
    detected_classes <- c(detected_classes, "q")
  }

  # 2. Check for hinge features ('h'): terms like hinge(var)
  if (any(grepl("^hinge\\(", formula_terms))) {
    detected_classes <- c(detected_classes, "h")
  }

  # 3. Check for threshold features ('t'): terms like thresholds(var)
  if (any(grepl("^thresholds\\(", formula_terms))) {
    detected_classes <- c(detected_classes, "t")
  }

  # 4. Check for product features ('p'): terms containing a colon, like var1:var2
  if (any(grepl(":", formula_terms))) {
    detected_classes <- c(detected_classes, "p")
  }

  # 5. Check for linear features ('l'). This is the trickiest.
  # A linear term is a plain variable name, not wrapped in a function
  # call (like hinge()) and not part of an interaction (like var1:var2).
  # We identify these by finding terms that DO NOT match any of the other patterns.

  # Identify all terms that are *not* simple linear terms
  is_complex_or_categorical <- grepl("^I\\(|^hinge\\(|^thresholds\\(|^categorical\\(|:", formula_terms)

  # If there are any terms left over, they must be linear terms.
  if (any(!is_complex_or_categorical)) {
    detected_classes <- c(detected_classes, "l")
  }

  # --- Finalize the output ---

  # Define the canonical order to ensure consistent output (e.g., "lqph" not "hqlp")
  class_order <- c("l", "q", "p", "h", "t")

  # Filter the detected classes and sort them according to the canonical order
  ordered_classes <- intersect(class_order, detected_classes)

  # Collapse the sorted vector into a single string
  return(paste(ordered_classes, collapse = ""))
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
kf <- function(df, n) {
  suppressWarnings(km <- stats::kmeans(df %>% dplyr::select(-c(cell:y)),
    centers = n
  ))
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
