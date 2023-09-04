#' Collinearity reduction of predictor variables
#'
#' @param env_layer SpatRaster An object of class SpatRaster containing the predictors.
#' This function does not allow categorical variables
#'
#' @param method character. Collinearity reduction method. It is necessary to
#' provide a vector for this argument. The next methods are implemented:
#' \itemize{
#'   \item pearson: Highlights correlated variables according to Pearson correlation. A threshold of maximum correlation
#'   must be specified. Otherwise, a threshold of 0.7 is defined as default.
#'   Usage method = c('pearson', th='0.7').
#'   \item vif: Select variables by Variance Inflation Factor, a threshold can be specified by
#'   user. Otherwise, a threshold of 10 is defined as default.Usage method = c('vif', th = '10').
#'   \item pca: Perform a Principal Component Analysis and use the principal components as the
#'   new predictors. The selected components account for 95\% of the whole variation in the system.
#'   Usage method = c('pca').
#'   \item fa: Perform a Factorial Analysis and select, from the original predictors, the number of factors is defined by Broken-Stick and variables with the highest correlation to the factors are selected.  Usage method = c('fa').
#' }
#' @param proj character. Only used for pca method. Path to a folder that contains sub-folders for the different projection
#' scenarios. Variables names must have the same names as in the raster used in env_layer argument. Usage proj = "C:/User/Desktop/Projections" (see in Details more about the use of this argument)
#' @param maxcell numeric. Number of raster cells to be randomly sampled. Taking a sample could be
#' useful to reduce memory usage for large rasters. If NULL, the function will use all
#' raster cells. Default NULL. Usage maxcell = 50000.
#'
#' @return
#' #' If 'pearson', returns a list with the following elements:
#' \itemize{
#' \item cor_table: a matrix object with pairwise correlation values of the environmental variables
#' \item cor_variables: a list object with the same length of the number of environmental values containing the pairwise relations that exceeded the correlation threshold for each one of the environmental variables
#' }
#'
#' If 'vif' method, returns a list with the following elements:
#' \itemize{
#' \item env_layer: a SpatRaster object with selected environmental variables
#' \item removed_variables: a character vector with removed environmental variables
#' \item vif_table: a data frame with VIF values for all environmental variables
#' }
#'
#' If 'pca' method, returns a list with the following elements:
#' \itemize{
#' \item env_layer: SpatRaster with scores of selected principal component (PC) that sum up 95\% of the
#' whole variation or original environmental variables
#' \item coefficients: a matrix with the coefficient of principal component (PC) for predictors
#' \item cumulative_variance: a tibble with the cumulative variance explained in selected principal component (PC)
#' }
#'
#' If 'fa' method, returns a list with the following elements:
#' \itemize{
#' \item env_layer: SpatRaster with scores of selected variables due to correlation to factors.
#' \item number_factors: number of factors selected according to the Broken-Stick criteria,
#' \item removed_variables: removed variables,
#' \item uniqueness: uniqueness of each environmental variable according to the factorial analysis,
#' \item loadings: environmental variables loadings in each of the chosen factors
#' }
#'
#' @details In the case of having environmental variables for the current conditions and other time
#' periods (future or present), it is recommended to perform the PCA analysis with the current
#' environmental condition and project the PCA for the other time periods. To do so, it is necessary
#' to use “proj” argument. Path to a folder (e.g., projections) that contains sub-folders for the
#' different projection scenarios (e.g., years and emissions). Within each sub-folder must be stored
#' single or multiband rasters with the environmental variables.
#'
#' For example:
#'
#' C:/Users/my_pc/projections/ \cr
#'     ├── MRIESM_2050_ssp126 \cr
#'     │   └── var1.tif\cr
#'     │   └── var2.tif\cr
#'     │   └── var3.tif\cr
#'     ├── MRIESM_2080_ssp585\cr
#'     │   └── var1.tif\cr
#'     │   └── var2.tif\cr
#'     │   └── var3.tif\cr
#'     ├── UKESM_2050_ssp370\cr
#'     │   └── var1.tif\cr
#'     │   └── var2.tif\cr
#'     │   └── var3.tif
#'
#' If pca method is run with time projections, correct_colinvar function will create the
#' Projection_PCA (the exact path is in the path object returned by the function) with the same
#' system of sub-folders and multiband raster with the principal components (pcs.tif)
#'
#' C:/Users/my_pc/Projection_PCA/\cr
#'       ├── MRIESM_2050_ssp126\cr
#'       │   └── pcs.tif           # a multiband tif with principal components\cr
#'       ├── MRIESM_2080_ssp585\cr
#'       │   └── pcs.tif\cr
#'       ├── UKESM_2050_ssp370\cr
#'       │   └── pcs.tif
#'
#' @export
#' @importFrom dplyr tibble
#' @importFrom stats na.omit cor lm prcomp factanal
#' @importFrom terra rast as.data.frame subset predict scale writeRaster global
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#'
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' # Perform pearson collinearity control
#' var <- correct_colinvar(env_layer = somevar, method = c("pearson", th = "0.7"))
#' var$cor_table
#' var$cor_variables
#'
#' # For all correct_colinvar methods it is possible to take a sample or raster to reduce memory
#' var <- correct_colinvar(env_layer = somevar, method = c("pearson", th = "0.7"), maxcell = 10000)
#' var$cor_table
#' var$cor_variables
#'
#' # Perform vif collinearity control
#' var <- correct_colinvar(env_layer = somevar, method = c("vif", th = "8"))
#' var$env_layer
#' var$removed_variables
#' var$vif_table
#'
#' # Perform pca collinearity control
#' var <- correct_colinvar(env_layer = somevar, method = c("pca"))
#' plot(var$env_layer)
#' var$env_layer
#' var$coefficients
#' var$cumulative_variance
#'
#'
#' # Perform pca collinearity control with different projections
#' ## Below will be created a set of folders to simulate the structure of the directory where
#' ## environmental variables are stored for different scenarios
#' dir_sc <- file.path(tempdir(), "projections")
#' dir.create(dir_sc)
#' dir_sc <- file.path(dir_sc, c('scenario_1', 'scenario_2'))
#' sapply(dir_sc, dir.create)
#'
#' somevar <-
#'   system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' terra::writeRaster(somevar, file.path(dir_sc[1], "somevar.tif"), overwrite=TRUE)
#' terra::writeRaster(somevar, file.path(dir_sc[2], "somevar.tif"), overwrite=TRUE)
#'
#' ## Perform pca with projections
#' dir_w_proj <- dirname(dir_sc[1])
#' dir_w_proj
#' var <- correct_colinvar(env_layer = somevar, method = "pca", proj = dir_w_proj)
#' var$env_layer
#' var$coefficients
#' var$cumulative_variance
#' var$proj
#'
#'
#' # Perform fa colinearity control
#' var <- correct_colinvar(env_layer = somevar, method = c("fa"))
#' var$env_layer
#' var$number_factors
#' var$removed_variables
#' var$uniqueness
#' var$loadings
#' }
#'
correct_colinvar <- function(env_layer,
                             method,
                             proj = NULL,
                             maxcell = NULL) {
  . <- NULL
  if (!any(c("pearson", "vif", "pca", "fa") %in% method)) {
    stop(
      "argument 'method' was misused, select one of the available methods: pearson, vif, pca, fa"
    )
  }

  if (class(env_layer)[1] != "SpatRaster") {
    env_layer <- terra::rast(env_layer)
  }

  if (any(method %in% "pearson")) {
    if (is.na(method["th"])) {
      th <- 0.7
    } else {
      th <- as.numeric(method["th"])
    }

    if(is.null(maxcell)){
      h <- terra::as.data.frame(env_layer) %>% stats::na.omit()
    } else {
      # Raster random sample
      set.seed(10)
      h <- terra::as.data.frame(env_layer[[1]], cells=TRUE)[,1] %>%
        sample(., size = maxcell, replace = FALSE) %>%
        sort()
      h <- env_layer[h] %>%
        stats::na.omit()
    }
    h <- abs(stats::cor(h, method = "pearson"))
    diag(h) <- 0

    cor_var <- h>th
    cor_var <- apply(cor_var,2, function(x) colnames(h)[x])
    if(length(cor_var)==0){
      cor_var <- 'No pair of variables reached the specified correlation threshold.'
    }

    result <- list(
      cor_table = h,
      cor_variables = cor_var
    )
  }

  if (any(method %in% "vif")) {
    if (is.null(method["th"])) {
      th <- 10
    } else {
      th <- as.numeric(method["th"])
    }

    if(is.null(maxcell)){
      x <- terra::as.data.frame(env_layer)
    } else {
      # Raster random sample
      set.seed(10)
      x <- terra::as.data.frame(env_layer[[1]], cells=TRUE)[,1] %>%
        sample(., size = maxcell, replace = FALSE) %>%
        sort()
      x <- env_layer[x] %>%
        stats::na.omit()
    }

    LOOP <- TRUE
    if (nrow(x) > 200000) {
      x <- x[sample(1:nrow(x), 200000), ]
    }
    n <- list()
    n$variables <- colnames(x)
    exc <- c()

    while (LOOP) {
      v <- rep(NA, ncol(x))
      names(v) <- colnames(x)
      for (i in 1:ncol(x)) {
        suppressWarnings(v[i] <- 1 / (1 - summary(lm(x[, i] ~ ., data = x[-i]))$r.squared))
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

    # n$corMatrix <- stats::cor(x, method = "pearson")
    n$results <- data.frame(Variables = names(v), VIF = as.vector(v))

    # diag(n$corMatrix) <- 0
    env_layer <-
      terra::subset(env_layer, subset = n$results$Variables)

    result <- list(
      env_layer = env_layer,
      removed_variables = n$excluded,
      vif_table = dplyr::tibble(n$results)
    )
  }

  if (any(method %in% "pca")) {

    # mean
    means <- t(terra::global(env_layer, 'mean', na.rm=T)) %>% c()
    names(means) <- names(env_layer)
    # SD
    stds <- t(terra::global(env_layer, 'sd', na.rm=T)) %>% c()
    names(stds) <- names(env_layer)

    # Standardize raster values
    env_layer <- terra::scale(env_layer, center = means, scale = stds)
    vnmes <- names(means)


    if(is.null(maxcell)){
      p0 <- terra::as.data.frame(env_layer, xy = FALSE, na.rm = TRUE)
    } else {
      # Raster random sample
      set.seed(10)
      p0 <- terra::as.data.frame(env_layer[[1]], cells=TRUE)[,1] %>%
        sample(., size = maxcell, replace = FALSE) %>%
        sort()
      p0 <- env_layer[p0] %>%
        stats::na.omit()
    }

    p <- stats::prcomp(p0,
                         retx = TRUE,
                         scale. = FALSE,
                         center = FALSE
                       )
    cof <- p$rotation

    cvar <- summary(p)$importance["Cumulative Proportion", ]
    naxis <- Position(function(x) {
      x >= 0.95
    }, cvar)
    cvar <- data.frame(cvar)


    # p <- terra::as.data.frame(env_layer, xy = FALSE, na.rm = TRUE)
    p <- stats::prcomp(p0, retx = TRUE, scale. = FALSE, center = FALSE, rank. = naxis)
    env_layer <- terra::predict(env_layer, p)

    rm(p0)

    result <- list(
      env_layer = env_layer,
      coefficients = data.frame(cof) %>% dplyr::tibble(variable = rownames(.), .),
      cumulative_variance = dplyr::tibble(PC = 1:nrow(cvar), cvar)
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
        scen <- scen[[names(means)]]
        scen <- terra::scale(scen, center = means, scale = stds)
        scen <- terra::predict(scen, p)
        terra::writeRaster(
          scen,
          file.path(subfold[[i]], "pcs.tif"),
          overwrite=TRUE
        )
      }

      result$proj <- dpca
    }
  }

  if (any(method %in% "fa")) {
    p <- terra::scale(env_layer, center = TRUE, scale = TRUE)

    if(is.null(maxcell)){
      p <- terra::as.data.frame(p, xy = FALSE, na.rm = TRUE)
    } else {
      # Raster random sample
      set.seed(10)
      p <- terra::as.data.frame(env_layer[[1]], cells=TRUE)[,1] %>%
        sample(., size = maxcell, replace = FALSE) %>%
        sort()
      p <- env_layer[p] %>%
        stats::na.omit()
    }

    if (nrow(p) > 200000) {
      p <- p[sample(1:nrow(p), 200000), ]
    }

    e <- eigen(stats::cor(p))
    len <- length(e$values)
    a <- NULL
    r <- NULL

    for (j in 1:len) {
      a[j] <- 1 / len * sum(1 / (j:len))
      r[j] <- e$values[j] / (sum(e$values))
    }

    ns <- length(which(r > a))

    fit <-
      tryCatch(
        stats::factanal(
          x = p,
          factors = ns,
          rotation = "varimax",
          lower = 0.001
        ),
        error = function(e) {
          stats::factanal(
            x = p,
            factors = ns - 1,
            rotation = "varimax",
            lower = 0.001
          )
        }
      )

    sel <-
      row.names(fit$loadings)[apply(fit$loadings, 2, which.max)]
    rem <-
      row.names(fit$loadings)[!row.names(fit$loadings) %in% sel]

    env_layer <- terra::subset(env_layer, sel)

    h <- fit$loadings %>%
      matrix() %>%
      data.frame()
    colnames(h) <- paste("Factor", 1:ncol(h), sep = "_")

    result <- list(
      env_layer = env_layer,
      number_factors = fit$factors,
      removed_variables = rem,
      uniqueness = fit$uniquenesses,
      loadings = dplyr::tibble(Variable = rownames(fit$loadings), h)
    )
  }

  return(result)
}
