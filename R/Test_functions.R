##%######################################################%##
#                                                          #
####                 Test some function                 ####
#                                                          #
##%######################################################%##


require(ENMTML)
require(raster)
require(dplyr)
require(sf)

# Basic layer database
# r_base <- "C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/StudyAreaRaster.tif" %>%
#   raster::raster()
# dir.create("./Data")
# save(r_base, file=file.path(here::here('Data', 'basic_layer.RData')), compress = 'bzip2')
load(here::here('Data', 'basic_layer.RData'))

# Create a presences absences database 
set.seed(10)
sp_db <- dismo::randomPoints(r_base, n = 1000) %>%
  data.frame()
plot(r_base)
points(sp_db)

filt <- raster::extract(r_base, sp_db)
sp_db$pr_ab
sp_db <- sp_db %>% 
  mutate(pr_ab = 0) %>% 
  mutate(pr_ab = ifelse(filt > quantile(filt)[4], 1, 0))

sp_db <- data.frame(species='sp1', sp_db)
plot(r_base, col=gray.colors(10))
points(sp_db[,2:3], col=c('red','blue')[sp_db$pr_ab+1], pch=19)


ggplot()+
  geom_sf(data = grid[[i]])


shapefile(grid[[i]], 'test.shp')
