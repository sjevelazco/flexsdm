
evaluate <- function(){
  
}


ecospat.boyce(RastPart[["BIO"]][[i]],PredPoint[PredPoint$PresAbse==1,2],PEplot=F)$Spearman.cor
ecospat::ecospat.boyce
require(ecospat)

obs <- (ecospat.testData$glm_Saxifraga_oppositifolia
        [which(ecospat.testData$Saxifraga_oppositifolia==1)])

ecospat.boyce(fit = ecospat.testData$glm_Saxifraga_oppositifolia , obs,
               window.w="default")


e <- ecospat.boyce(fit = c(p, a), obs=rep(c(1,0), each=50))

cor(ecospat.testData$glm_Saxifraga_oppositifolia , obs)
require(dismo)

p <- rnorm(50, mean=0.7, sd=0.3)
# a has the predicted values for 50 background locations (or absence)
a <- rnorm(50, mean=0.4, sd=0.4)
e <- dismo::evaluate(p=p, a=a)
e
plot(e, "ROC")
plot(e, "TPR")
plot(e, "TNR")
e@cor
TSS <- ((e@TPR + e@TNR) - 1)
dismo::threshold(e) # add LPT
?threshold
e <- evaluate(p=p, a=a, tr = 0.5)
e <- evaluate(p=p, a=a)


dismo::evaluate(PredPoint[PredPoint$PresAbse == 1, 2],
                PredPoint[PredPoint$PresAbse == 0, 2],tr=Thr)

Eval_Jac_Sor_TMLA <- function(p,
                              a,
                              thr=NULL) {
  #Evaluate ENMs using Jaccard and Sorensen Indexes
  
  #Parameters:
  #p:presence points suitability
  #a:absence points suitability
  #tr:numeric vector with threshold values
  #Initialisation
  np <- length(p)
  na <- length(a)
  if (na == 0 | np == 0) {
    stop('cannot evaluate a model without absence and presence data that are not NA')
  }
  
  #Threshold breaks:
  if(is.null(thr)){
    if (length(p) > 1000) {
      tr <- as.vector(quantile(p, 0:1000 / 1000))
    } else {
      tr <- p
    }
    if (length(a) > 1000) {
      tr <- c(tr, as.vector(quantile(a, 0:1000 / 1000)))
    } else {
      tr <- c(tr, a)
    }
    tr <- sort(unique(round(tr, 8)))
    tr <- c(tr - 0.0001, tr[length(tr)] + c(0, 0.0001))
  }else{
    tr <- thr
  }
  
  res <- matrix(ncol = 4, nrow = length(tr))
  colnames(res) <- c('tp', 'fp', 'fn', 'tn')
  #Confusion Matrix
  for (i in 1:length(tr)) {
    res[i, 1] <- length(p[p >= tr[i]])  # a  true positives
    res[i, 2] <- length(a[a >= tr[i]])  # b  false positives
    res[i, 3] <- length(p[p < tr[i]])    # c  false negatives
    res[i, 4] <- length(a[a < tr[i]])    # d  true negatives
  }
  a = res[, 1]
  b = res[, 2]
  c = res[, 3]
  d = res[, 4]
  
  #Sorensen Index
  SOR <- 2 * a / (c + 2 * a + b)
  if(is.null(thr)){
    SorTHR <- tr[which(SOR == max(SOR))][1]  
  }else{
    SorTHR <- tr  
  }
  
  #Jaccard Index
  JAC <- a / (c + a + b)
  if(is.null(thr)){
    JacTHR <- tr[which(JAC == max(JAC))][1]
  }else{
    JacTHR <- tr  
  }
  
  #Fpb
  Fpb <- 2 * JAC
  if(is.null(thr)){
    FpbTHR <- tr[which(Fpb == max(Fpb))][1]
  }else{
    FpbTHR <- tr  
  }
  
  #Final Result Object
  xc <- list()
  xc$presences <- np
  xc$absences <- na
  xc$Sorensen <- SOR
  xc$SorensenTHR <- SorTHR
  xc$Jaccard <- JAC
  xc$JaccardTHR <- JacTHR
  xc$Fpb <- Fpb
  xc$t <- tr
  xc$TPR <- a / (a + c)
  xc$TNR <- d / (b + d)
  
  #Return final result
  return(xc)
}

Validation<-Validation_Table_TMLA(Eval=Eval,Eval_JS=Eval_JS,N=N,Thr=Threshold)