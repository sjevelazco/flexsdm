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
  braquets = variables[[i]][] %>% na.omit() %>% c(),
  getvalues = raster::getValues(variables[[i]]) %>% na.omit() %>% c(),
  times = 10)
autoplot(results)
results %>% class
results[[1]]
