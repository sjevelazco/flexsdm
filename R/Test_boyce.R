source("./R/boyce.R")
source("https://raw.githubusercontent.com/adamlilith/enmSdm/master/R/contBoyce.r")
require(ggplot2)
require(ecospat)

# compare stability of CBI calculated with ecospat.boyce() in ecospat package
library(ecospat)
set.seed(123)
results <- data.frame()
for (perform in c(1, 1.5, 2)) {
  for (i in 1:100) {
    print(i)
    set.seed(i)
    pres <- runif(100)^(1 / perform)
    set.seed(i)
    contrast <- runif(1000)
    
    cbi_enmSdm <- contBoyce(pres, contrast)
    cbi_adapted <- boyce(pres = pres, contrast = contrast, n_bins = 101, n_width = 0.1)
    cbi_ecospat <- ecospat.boyce(contrast, pres, PEplot=FALSE)$Spearman.cor
    
    results <- rbind(
      results,
      data.frame(
        performance = rep(perform, 3),
        method = c('enmSdm', 'cbi_adapted', 'ecospat'),
        cbi = c(cbi_enmSdm, cbi_adapted, cbi_ecospat)
      )
    )
  }
}

results$performance[results$performance == 1] <- 'poor'
results$performance[results$performance == 1.5] <- 'OK'
results$performance[results$performance == 2] <- 'good'



results$category <- paste0(results$method, '\n', results$performance)
require(ggplot2)
ggplot(results, aes(method, cbi, fill=performance))+geom_boxplot()

results$cbi %>% hist
freq()
results$cbi
# Difference between adapted function and enmSdm
results[results$method=="cbi_adapted",'cbi']%in%
  results[results$method=="enmSdm" ,'cbi'] # No difference between to function 
plot(results[results$method=="cbi_adapted",'cbi']-
       results[results$method=="enmSdm" ,'cbi'])

# Difference between adapted function and ecospat
results[results$method=="cbi_adapted",'cbi']-
  results[results$method=="ecospat" ,'cbi']
hist(results[results$method=="cbi_adapted",'cbi']-
        results[results$method=="ecospat" ,'cbi']
)
summary(results[results$method=="cbi_adapted",'cbi']-
       results[results$method=="ecospat" ,'cbi']
)

plot(results[results$method=="cbi_adapted",'cbi']-
  results[results$method=="ecospat" ,'cbi'])

par(mfrow=c(1,1))
plot(results[results$method=="cbi_adapted",'cbi'], results[results$method=="ecospat" ,'cbi'])

