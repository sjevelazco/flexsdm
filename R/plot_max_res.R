plot_max_res <- function(r, max_res_mult){
  r0 <- r
  r0[!is.na(values(r0))] <- 1
  
  r[] <- 0
  res(r) <- res(r)*max_res_mult 
  plot(r0)
  plot(rasterToPolygons(r), add=TRUE)
}
