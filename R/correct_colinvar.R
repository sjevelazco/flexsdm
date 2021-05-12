#' Perform collinearity reduction on predictors
#'
#' @param rstack SpatRaster An object of class SpatRaster containing the predictors.
#' @param method character. Collinearity reduction method. It is necessary to
#' provide a vector for this argument. The next methods are implemented:
#' \itemize{
#'   \item pearson: Select variables by Pearson correlation, a threshold of maximum correlation must be specified by user. Otherwise a threshold of 0.7 is defined as default. Usage method = c('pearson', th='0.7').
#'   \item vif: Select variables by Variance Inflation Factor, a threshold can be specified by user. Otherwise a threshold of 10 is defined as default.Usage method = c('vif', th = '10').
#'   \item pca: Perform a Principal Component Analysis and use the principal components as the new predictors. The selected components account for 95% of the whole variation in the system.Usage method = c('pca').
#'   \item fa: Perform a Factorial Analysis and select, from the original predictors, the ones with the highest correlation with the axis which account for 95% of the whole variation in the system.Usage method = c('fa').
#' }
#' @param proj character. Path to a folder which contains sub-folders for the different projection scenarios. Only used for pca. Usage "C:/User/Desktop/Projections"
#'
#' @return
#' @export
#'
#' @importFrom stats cor prcomp factanal
#' @importFrom terra as.data.frame subset predict rast scale writeRaster
#'
#' @examples
#' \dontrun{
#' somevar <- system.file("external/somevar.tif", package = "flexsdm")
#' somevar <- terra::rast(somevar)
#'
#' #Perform pearson collinearity control
#' somevar <- correct_colinvar(rstack = somevar, method = c('pearson', th='0.8'))
#' somevar$rstack
#' somevar$removed_variables
#' somevar$correlation_table
#'
#' #Perform vif collinearity control
#' somevar <- correct_colinvar(rstack = somevar, method = c('vif', th='8'))
#' somevar$rstack
#' somevar$removed_variables
#' somevar$correlation_table
#'
#' #Perform pca collinearity control
#' somevar <- correct_colinvar(rstack = somevar, method = c('pca'))
#' somevar$rstack
#' somevar$coeficients
#' somevar$cumulative_variance
#'
#' #Perform pca collinearity control for projections
#' somevar <- correct_colinvar(rstack = somevar, method = c('pca', proj = ))
#' somevar$rstack
#' somevar$coeficients
#' somevar$cumulative_variance
#' somevar$proj
#'
#' #Perform fa collinearity control
#' somevar <- correct_colinvar(rstack = somevar, method = c('fa'))
#' somevar$rstack
#' somevar$removed_variables
#' somevar$correlation_table
#' }
#'
correct_colinvar <- function(rstack,
                             method,
                             proj = NULL) {
  if (!any(c("pearson", "vif", "pca", "fa") %in% method)) {
    stop(
      "argument 'method' was misused, select one of the available methods: pearson, vif, pca, fa"
    )
  }

  if (!class(rstack) %in% 'SpatRaster') {
    stop("Raster object must be from the class SpatRaster")
  }

  if (any(method %in% "pearson")) {
    if (is.null(method['th'])) {
      th <- 0.7
    } else{
      th <- as.numeric(method['th'])
    }

    h <- terra::as.data.frame(rstack)
    h <- base::abs(stats::cor(h, method = 'pearson'))
    diag(h) <- 0

    res <- as.list(1:10000)
    for (i in 1:10000) {
      ord <- sample(1:ncol(h))
      h2 <- h[ord, ord]
      h2[upper.tri(h2)] <- 0
      res[[i]] <-
        colnames(h2)[!apply(h2, 2, function(x)
          any(x > th))]
    }

    len <- sapply(res, function(x)
      length(x))
    sel <- res[[sample(which(len == max(len)), 1)]]
    rem <- names(rstack)[!names(rstack) %in% sel]
    rstack <- terra::subset(rstack, subset = sel)

    result <- list(
      rstack = rstack,
      removed_variables = rem,
      correlation_table = h
    )
  }

  if (any(method %in% "vif")) {
    if (is.null(method['th'])) {
      th <- 10
    } else{
      th <- as.numeric(method['th'])
    }

    x <- terra::as.data.frame(rstack)
    LOOP <- TRUE
    if(nrow(x) > 10000){
      x <- x[sample(1:nrow(x),10000),]
    }
    n <- list()
    n$variables <- colnames(x)
    exc <- c()

    while (LOOP) {
      v<-rep(NA,ncol(x))
      names(v) <- colnames(x)
      for (i in 1:ncol(x)) {
        v[i] <-  1/(1-summary(lm(x[,i]~.,data=x[-i]))$r.squared)
      }
      if (v[which.max(v)] >= th) {
        ex <- names(v[which.max(v)])
        exc <- c(exc,ex)
        x <- x[,-which(colnames(x) == ex)]
      } else {
        LOOP=FALSE
      }
    }
    if (length(exc) > 0){
      n$excluded <- exc
    }

    v<-rep(NA,ncol(x))
    names(v) <- colnames(x)
    for (i in 1:ncol(x)) {
      v[i] <-  1/(1-summary(lm(x[,i]~.,data=x[-i]))$r.squared)
    }

    n$corMatrix <- cor(x, method="pearson")
    n$results <- data.frame(Variables=names(v),VIF=as.vector(v))

    diag(n$corMatrix) <- 0
    rstack <-
      terra::subset(rstack, subset = n$results$Variables)

    result <- list(
      rstack = rstack,
      removed_variables = n$excluded,
      correlation_table = n$corMatrix
    )
  }

  if (any(method %in% "pca")) {
    p <- terra::as.data.frame(rstack, xy = FALSE, na.rm = TRUE)

    p <- stats::prcomp(p,
                       retx = TRUE,
                       scale. = TRUE,
                       center = TRUE)

    means <- p$center
    stds <- p$scale
    cof <- p$rotation

    cvar <- summary(p)$importance["Cumulative Proportion",]
    naxis <- Position(function(x)
      x >= 0.95, cvar)
    cvar <- data.frame(cvar)
    rstack <- terra::predict(rstack, p, index = 1:naxis)

    result <- list(
      rstack = rstack,
      coeficients = cof,
      cumulative_variance = cvar
    )

    if (!is.null(proj)) {
      dpca <- file.path(dirname(proj), 'Projection_PCA')
      dir.create(dpca)
      subfold <- list.files(proj)
      subfold <- as.list(file.path(dpca, subfold))
      sapply(subfold, function(x)
        dir.create(x))

      proj <- base::list.files(proj, full.names = TRUE)
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
                     proj = dpca)
    }
  }

  if (any(method %in% "fa")) {
    p <- terra::scale(rstack, center = TRUE, scale = TRUE)
    p <- terra::as.data.frame(p, xy = FALSE, na.rm = TRUE)

    if(nrow(p) > 10000){
      p <- p[sample(1:nrow(p),10000),]
    }

    e <- eigen(cor(p))
    len <- length(e$values)
    a <- NULL
    r <- NULL

    for(j in 1:len){
      a[j] <- 1/len * sum(1/(j:len) )
      r[j] <- e$values[j]/(sum(e$values))
    }

    ns <- length(which(r > a))

    fit <-
      tryCatch(
        stats::factanal(
          x = p,
          factors = nS$Components$noc,
          rotation = "varimax",
          lower = 0.01
        ),
        error = function(e)
          stop(
            "Covariance matrix is not invertible. Consider choosing another method to control collinearity.",
            call. =F)
      )

    sel <-
      row.names(fit$loadings)[apply(fit$loadings, 2, which.max)]
    rem <-
      row.names(fit$loadings)[!row.names(fit$loadings) %in% sel]

    rstack <- terra::subset(rstack, sel)

    result <- list(
      rstack = rstack,
      removed_variables = rem,
      correlation_table = fit$loadings
    )
  }

  return(result)
}



