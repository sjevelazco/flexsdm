#' Perform collinearity reduction on predictors
#'
#' @param env_layer SpatRaster An object of class SpatRaster containing the predictors.
#' This function does not allow categorical variables
#'
#' @param method character. Collinearity reduction method. It is necessary to
#' provide a vector for this argument. The next methods are implemented:
#' \itemize{
#'   \item pearson: Select variables by Pearson correlation, a threshold of maximum correlation
#'   must be specified. Otherwise, a threshold of 0.7 is defined as default.
#'   Usage method = c('pearson', th='0.7').
#'   \item vif: Select variables by Variance Inflation Factor, a threshold can be specified by
#'   user. Otherwise, a threshold of 10 is defined as default.Usage method = c('vif', th = '10').
#'   \item pca: Perform a Principal Component Analysis and use the principal components as the
#'   new predictors. The selected components account for 95% of the whole variation in the system.
#'   Usage method = c('pca').
#'   \item fa: Perform a Factorial Analysis and select, from the original predictors, the ones with
#'    the highest correlation with the axis which accounts for 95% of the whole variation in the
#'    system.Usage method = c('fa').
#' }
#' @param proj character. Path to a folder that contains sub-folders for the different projection
#' scenarios. Only used for pca. Usage proj = "C:/User/Desktop/Projections"
#'
#' @return
#' If it used 'pearson' or 'vif' method, it is returned a list with the next object:
#' \itemize{
#' \item env_layer: a SpatRaster object with selected environmental variables
#' \item removed_variables: a character vector with removed environmental variables
#' \item performance: a matrix with Pearson pairwise correlation between variables
#' }
#'
#' If it used 'pca' method, it is returned a list with the next object:
#' \itemize{
#' \item env_layer: SpatRaster with scores of selected principal component (PC) that sum up 95% of the
#' whole variation or original environmental variables
#' \item coefficient: a matrix with the coefficient of principal component (PC) for predictors
#' \item cumulative_variance: a tibble with the cumulative variance explained in selected principal component (PC)
#' }
#'
#' If it used 'fa' method, it is returned a list with the next object:
#' \itemize{
#' \item env_layer: #TODO
#' \item coefficient:: #TODO
#' \item cumulative_variance:: #TODO
#' }
#'
#' @export
#' @importFrom stats cor lm prcomp factanal
#' @importFrom terra as.data.frame subset predict rast scale writeRaster
#' @importFrom dplyr tibble
#'
#' @examples
#' \dontrun{
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' # Perform pearson collinearity control
#' var <- correct_colinvar(env_layer = somevar, method = c("pearson", th = "0.8"))
#' var$env_layer
#' var$removed_variables
#' var$correlation_table
#'
#' # Perform vif collinearity control
#' var <- correct_colinvar(env_layer = somevar, method = c("vif", th = "8"))
#' var$env_layer
#' var$removed_variables
#' var$correlation_table
#'
#' # Perform pca collinearity control
#' var <- correct_colinvar(env_layer = somevar, method = c("pca"))
#' plot(var$env_layer)
#' var$env_layer
#' var$coefficient
#' var$cumulative_variance
#'
#' # Perform fa colinearity control
#' # this method only will be performed if covariance matrix is invertible.
#' # WRITE HERE A EXAMPLE THAT WORKS :)
#' # var <- correct_colinvar(env_layer = somevar, method = c("fa"))
#' }
#'
correct_colinvar <- function(env_layer,
                             method,
                             proj = NULL) {
  #TODO write documentation for fa methods and write details in methods!!!
  if (!any(c("pearson", "vif", "pca", "fa") %in% method)) {
    stop(
      "argument 'method' was misused, select one of the available methods: pearson, vif, pca, fa"
    )
  }

  if (!class(env_layer) %in% "SpatRaster") {
    stop("Raster object must be from the class SpatRaster")
  }

  if (any(method %in% "pearson")) {
    if (is.null(method["th"])) {
      th <- 0.7
    } else {
      th <- as.numeric(method["th"])
    }

    h <- terra::as.data.frame(env_layer)
    h <- abs(stats::cor(h, method = "pearson"))
    diag(h) <- 0

    res <- as.list(1:10000)
    for (i in 1:10000) {
      ord <- sample(1:ncol(h))
      h2 <- h[ord, ord]
      h2[upper.tri(h2)] <- 0
      res[[i]] <-
        colnames(h2)[!apply(h2, 2, function(x) {
          any(x > th)
        })]
    }

    len <- sapply(res, function(x) {
      length(x)
    })
    sel <- res[[sample(which(len == max(len)), 1)]]
    rem <- names(env_layer)[!names(env_layer) %in% sel]
    env_layer <- terra::subset(env_layer, subset = sel)

    result <- list(
      env_layer = env_layer,
      removed_variables = rem,
      correlation_table = h
    )
  }

  if (any(method %in% "vif")) {
    if (is.null(method["th"])) {
      th <- 10
    } else {
      th <- as.numeric(method["th"])
    }

    x <- terra::as.data.frame(env_layer)
    LOOP <- TRUE
    if (nrow(x) > 10000) {
      x <- x[sample(1:nrow(x), 10000), ]
    }
    n <- list()
    n$variables <- colnames(x)
    exc <- c()

    while (LOOP) {
      v <- rep(NA, ncol(x))
      names(v) <- colnames(x)
      for (i in 1:ncol(x)) {
        v[i] <- 1 / (1 - summary(lm(x[, i] ~ ., data = x[-i]))$r.squared)
      }
      if (v[which.max(v)] >= th) {
        ex <- names(v[which.max(v)])
        exc <- c(exc, ex)
        x <- x[, -which(colnames(x) == ex)]
      } else {
        LOOP <- FALSE
      }
    }
    if (length(exc) > 0) {
      n$excluded <- exc
    }

    v <- rep(NA, ncol(x))
    names(v) <- colnames(x)
    for (i in 1:ncol(x)) {
      v[i] <- 1 / (1 - summary(stats::lm(x[, i] ~ ., data = x[-i]))$r.squared)
    }

    n$corMatrix <- stats::cor(x, method = "pearson")
    n$results <- data.frame(Variables = names(v), VIF = as.vector(v))

    diag(n$corMatrix) <- 0
    env_layer <-
      terra::subset(env_layer, subset = n$results$Variables)

    result <- list(
      env_layer = env_layer,
      removed_variables = n$excluded,
      correlation_table = n$corMatrix
    )
  }

  if (any(method %in% "pca")) {
    p <- terra::as.data.frame(env_layer, xy = FALSE, na.rm = TRUE)

    p <- stats::prcomp(p,
      retx = TRUE,
      scale. = TRUE,
      center = TRUE
    )

    means <- p$center
    stds <- p$scale
    cof <- p$rotation

    cvar <- summary(p)$importance["Cumulative Proportion", ]
    naxis <- Position(function(x) {
      x >= 0.95
    }, cvar)
    cvar <- data.frame(cvar)
    env_layer <- terra::predict(env_layer, p, index = 1:naxis)

    result <- list(
      env_layer = env_layer,
      coefficients = cof,
      cumulative_variance = dplyr::tibble(cvar)
    )

    if (!is.null(proj)) {
      dpca <- file.path(dirname(proj), "Projection_PCA")
      dir.create(dpca)
      subfold <- list.files(proj)
      subfold <- as.list(file.path(dpca, subfold))
      sapply(subfold, function(x) {
        dir.create(x)
      })

      proj <- list.files(proj, full.names = TRUE)
      for (i in 1:length(proj)) {
        scen <- terra::rast(list.files(proj[i], full.names = TRUE))
        scen <- terra::scale(scen, center = means, scale = stds)
        scen <- predict(scen, p, index = 1:naxis)
        terra::writeRaster(
          scen,
          file.path(subfold[i], paste0(names(scen), ".tif")),
          filetype = "GTiff",
          NAflag = -9999
        )
      }

      result <- list(result,
        proj = dpca
      )
    }
  }

  if (any(method %in% "fa")) {
    p <- terra::scale(env_layer, center = TRUE, scale = TRUE)
    p <- terra::as.data.frame(p, xy = FALSE, na.rm = TRUE)

    if (nrow(p) > 10000) {
      p <- p[sample(1:nrow(p), 10000), ]
    }

    e <- eigen(terra::predict(p))
    len <- length(e$values)
    a <- NULL
    r <- NULL

    for (j in 1:len) {
      a[j] <- 1 / len * sum(1 / (j:len))
      r[j] <- e$values[j] / (sum(e$values))
    }

    ns <- length(which(r > stats::cor))

    fit <-
      tryCatch(
        stats::factanal(
          x = p,
          factors = ns,
          rotation = "varimax",
          lower = 0.01
        ),
        error = function(e) {
          stop(
            "Covariance matrix is not invertible. Consider choosing another method to control collinearity.",
            call. = F
          )
        }
      )

    sel <-
      row.names(fit$loadings)[apply(fit$loadings, 2, which.max)]
    rem <-
      row.names(fit$loadings)[!row.names(fit$loadings) %in% sel]

    env_layer <- terra::subset(env_layer, sel)

    result <- list(
      env_layer = env_layer,
      removed_variables = rem,
      correlation_table = fit$loadings
    )
  }

  return(result)
}
