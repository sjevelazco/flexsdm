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



for (p in 1:length(grid)) {
  ncell3 <- ncell[,p]
  part3 <- c(part[,p])
  filt <- data.frame(
    nrow = 1:length(ncell3),
    ncell = ncell3,
    group = part3,
    pr_ab = presences2@data[c('pr_ab')]
  ) %>%
    dplyr::group_by(ncell, group, pr_ab) %>% 
    dplyr::slice_sample(n = 1) %>%
    dplyr::pull(nrow) %>% 
    sort()
  print(paste(p, length(filt)))
}
