require('microbenchmark')
require(ggplot2)

col <- names(env_layer)
rown <- 1:100

# cl <- parallel::makeCluster(cores, outfile = "")
# doParallel::registerDoParallel(cl)

im <- sapply(species2[rown, col],
       function(x)
         ape::Moran.I(x,
                      dist2[rown, rown],
                      na.rm = T,
                      scaled = T)$observed)

results <- microbenchmark(
  apply_ap = apply(species2[rown, col],2,
                     function(x)
                       ape::Moran.I(x,
                                    dist2[rown,rown],
                                    na.rm = T,
                                    scaled = T)$observed),
  sapply_ap = sapply(species2[rown, col],
                          function(x)
                            ape::Moran.I(x,
                                         dist2[rown,rown],
                                         na.rm = T,
                                         scaled = T)$observed),
  
  # mclapplya_ap=unlist(parallel::mclapply(species2[rown, col],
  #                    function(x){
  #                      ape::Moran.I(x,
  #                                   dist2[rown,rown],
  #                                   na.rm = T,
  #                                   scaled = T)$observed
  #                    },  mc.cores = 1)),
                     
  # foreach_ap = foreach(
  #   ci = 1:length(col),
  #   .packages = c("ape"),
  #   .final = unlist
  # ) %dopar% {
  #   ape::Moran.I(species2[rown, col[ci]],
  #                dist2[rown, rown],
  #                na.rm = T,
  #                scaled = T)$observed
  # },
  # 
  # do_ap = foreach(
  #   ci = 1:length(col),
  #   .packages = c("ape"),
  #   .final = unlist
  # ) %do% {
  #   ape::Moran.I(species2[rown, col[ci]],
  #                dist2[rown, rown],
  #                na.rm = T,
  #                scaled = T)$observed
  # }
  
  times = 25)


require(ggplot2)

results <- microbenchmark(
  for_m = {eval <- list()
  for(i in 1:length(pred_test)) {
    eval[[i]] <-
      enmtml_evaluate(p = pred_test[[i]]$pred[pred_test[[i]]$pr_ab == 1],
                      a = pred_test[[i]]$pred[pred_test[[i]]$pr_ab == 0],
                      thr = thr)
  }},
  
  sappli_m = lapply(pred_test, function(x) {
    enmtml_evaluate(p = x$pred[x$pr_ab == 1],
                    a = x$pred[x$pr_ab == 0],
                    thr = thr)
  }),
  times = 10)
autoplot(results)
results %>% class
results[[1]]
