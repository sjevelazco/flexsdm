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
#' @importFrom usdm vifstep exclude
#' @importFrom terra as.data.frame subset predict rast scale names
#' @importFrom nFactors parallel nScree
#'
#' @examples
#' \donotrun{
#' data('somevar')
#'
#' #Perform pearson collinearity control
#' h <- terra::as.data.frame(somevar)
#' h <- base::abs(stats::cor(h, method = 'pearson'))
#' diag(h) <- 0
#'
#' res <- as.list(1:10000)
#' for(i in 1:10000){
#'   ord <- base::sample(1:ncol(h))
#'   h2 <- h[ord,ord]
#'   h2[upper.tri(h2)] <- 0
#'   res[[i]] <- terra::names(somevar)[!apply(h2,2,function(x) any(x > 0.7))]
#' }
#'
#' len <- sapply(res, function(x) length(x))
#' sel <- res[[sample(which(len==max(len)),1)]]
#' rem <- terra::names(rstack)[!names(rstack)%in%sel]
#' rstack <- terra::subset(rstack,subset=sel)
#' result <- list(
#' rstack = rstack,
#' removed_variables = rem,
#' correlation_table = h
#' )
#'
#' #Perform vif collinearity control
#' VF <- usdm::vifstep(terra::as.data.frame(somevar), th = 10)
#' rem <- VF@excluded
#' h <- VF@corMatrix
#' diag(h) <- 0
#' rstack <- terra::subset(rstack,subset=VF@variables[!VF@variables%in%VF@excluded])
#' result <- list(
#' rstack = rstack,
#' removed_variables = rem,
#' correlation_table = h
#' )
#'
#' #Perform pca collinearity control
#' p <- terra::as.data.frame(somevar,xy=F,na.rm=TRUE)
#' p <- stats::prcomp(p, retx = T, scale. = T, center = T)
#' means <- p$center
#' stds <- p$scale
#' cof <- p$rotation
#' cvar <- summary(p)$importance["Cumulative Proportion", ]
#' naxis <- Position(function(x) x>=0.95,cvar)
#' cvar <- data.frame(cvar)
#' rstack <- terra::predict(rstack,p,index=1:naxis)
#' result <- list(
#' rstack = rstack,
#' coeficients = cof,
#' cumulative_variance = cvar
#' )
#'
#' #Perform fa collinearity control
#' p <- terra::scale(rstack,center = T, scale = T)
#' p <- terra::as.data.frame(p,xy=F,na.rm=TRUE)
#' e <- eigen(cor(p))
#' ap <- nFactors::parallel(subject=nrow(p),var=ncol(p),rep=100,cent=.05)
#' nS <- nFactors::nScree(x=e$values, aparallel=ap$eigen$qevpea)
#' fit <- tryCatch(stats::factanal(x = p, factors = nS$Components$noc, rotation = "varimax", lower = 0.01),
#' error = function(e) message("Warning: covariance matrix is not invertible. Consider removing highly correlated variables
#' and trying again or choosing another method to control collinearity."))
#' sel <- row.names(fit$loadings)[apply(fit$loadings,2,which.max)]
#' rem <- row.names(fit$loadings)[!row.names(fit$loadings) %in% sel]
#' rstack <- terra::subset(rstack,sel)
#' result <- list(
#' rstack = rstack,
#' removed_variables = rem,
#' correlation_table = fit$loadings
#' )
#'
#' }
#'
correct_colinvar <- function(rstack,
                             method,
                             proj = NULL){
  if (!any(c("pearson", "vif", "pca", "fa") %in% method)) {
    stop("argument 'method' was misused, select one of the available methods: pearson, vif, pca, fa")
  }

  if (!class(rstack)%in%'SpatRaster'){
    stop("Raster object must be from the class SpatRaster")
  }

  if (any(method %in% "pearson")) {
    if(is.null(method['th'])){
      th <- 0.7
    } else{
      th <- as.numeric(method['th'])
    }

    h <- terra::as.data.frame(rstack)
    h <- base::abs(stats::cor(h, method = 'pearson'))
    diag(h) <- 0

    res <- as.list(1:10000)
    for(i in 1:10000){
      ord <- sample(1:ncol(h))
      h2 <- h[ord,ord]
      h2[upper.tri(h2)] <- 0
      res[[i]] <- names(rstack)[!apply(h2,2,function(x) any(x > th))]
    }

    len <- sapply(res, function(x) length(x))
    sel <- res[[sample(which(len==max(len)),1)]]
    rem <- names(rstack)[!names(rstack)%in%sel]
    rstack <- terra::subset(rstack,subset=sel)

    result <- list(
      rstack = rstack,
      removed_variables = rem,
      correlation_table = h
    )
  }

  if (any(method %in% "vif")) {
    if(is.null(method['th'])){
      th <- 10
    } else{
      th <- as.numeric(method['th'])
    }

    VF <- usdm::vifstep(terra::as.data.frame(rstack), th = th)
    rem <- VF@excluded
    h <- VF@corMatrix
    diag(h) <- 0
    rstack <- terra::subset(rstack,subset=VF@variables[!VF@variables%in%VF@excluded])

    result <- list(
      rstack = rstack,
      removed_variables = rem,
      correlation_table = h
    )
  }

  if (any(method %in% "pca")) {
    p <- terra::as.data.frame(rstack,xy=F,na.rm=TRUE)

    p <- stats::prcomp(p, retx = T, scale. = T, center = T)

    means <- p$center
    stds <- p$scale
    cof <- p$rotation

    cvar <- summary(p)$importance["Cumulative Proportion", ]
    naxis <- Position(function(x) x>=0.95,cvar)
    cvar <- data.frame(cvar)
    rstack <- terra::predict(rstack,p,index=1:naxis)

    result <- list(
      rstack = rstack,
      coeficients = cof,
      cumulative_variance = cvar
    )

    if(!is.null(proj)){

      dpca <- file.path(dirname(proj), 'Projection_PCA')
      dir.create(dpca)
      subfold <- list.files(proj)
      subfold <- as.list(file.path(dpca, subfold))
      sapply(subfold, function(x) dir.create(x))

      proj <- base::list.files(proj,full.names = T)
      for(i in 1:length(proj)){
        scen <- terra::rast(list.files(proj[i],full.names = T))
        scen <- terra::scale(scen, center = means, scale = stds)
        scen <- predict(scen, p, index = 1:naxis)
        terra::writeRaster(scen,file.path(subfold[i],paste0(names(scen),".tif")),filetype = "GTiff", NAflag = -9999)
      }

      result <- list(
        result,
        proj = dpca
      )
    }
  }

  if (method%in%"fa"){
    p <- terra::scale(rstack,center = T, scale = T)
    p <- terra::as.data.frame(p,xy=F,na.rm=TRUE)

    e <- eigen(cor(p))
    ap <- nFactors::parallel(subject=nrow(p),var=ncol(p),
                   rep=100,cent=.05)
    nS <- nFactors::nScree(x=e$values, aparallel=ap$eigen$qevpea)
    fit <- tryCatch(stats::factanal(x = p, factors = nS$Components$noc, rotation = "varimax", lower = 0.01),
             error = function(e) message("Warning: covariance matrix is not invertible. Consider removing highly correlated variables and trying again or choosing another method to control collinearity."))

    sel <- row.names(fit$loadings)[apply(fit$loadings,2,which.max)]
    rem <- row.names(fit$loadings)[!row.names(fit$loadings) %in% sel]

    rstack <- terra::subset(rstack,sel)

    result <- list(
      rstack = rstack,
      removed_variables = rem,
      correlation_table = fit$loadings
    )
  }

  return(result)
}



