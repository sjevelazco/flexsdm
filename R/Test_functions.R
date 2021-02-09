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


# somevar <- "C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/CFP_environmental_stack.tif" %>%
  # raster::stack()
# somevar <- somevar[[c(1:4)]]
# save(somevar, file=file.path(here::here('Data', 'somevar.RData')), compress = 'bzip2')
# load(here::here('Data', 'somevar.RData'))

# Create a presences absences database 
load(here::here('Data', 'somevar.RData'))


# plot(somevar)

# Species 1
set.seed(10)
sp_db <- dismo::randomPoints(somevar, n = 1000) %>%
  data.frame()
filt <- raster::extract(somevar[[1]], sp_db)
sp_db$pr_ab
sp_db <- sp_db %>% 
  mutate(pr_ab = 0) %>% 
  mutate(pr_ab = ifelse(filt > quantile(filt)[4], 1, 0))
sp_db <- data.frame(species='sp1', sp_db)

# Species 2
set.seed(15)
sp_db2 <- dismo::randomPoints(somevar, n = 100) %>%
  data.frame()
filt <- raster::extract(somevar[[3]], sp_db2)

sp_db2 <- sp_db2 %>% 
  mutate(pr_ab = 0) %>% 
  mutate(pr_ab = ifelse(filt > quantile(filt)[4], 1, 0))
sp_db2 <- data.frame(species='sp2', sp_db2)

# Species 3
set.seed(13)
sp_db3 <- dismo::randomPoints(somevar, n = 50, ext = raster::extent(somevar)-500000) %>%
  data.frame()
filt <- raster::extract(somevar[[2]], sp_db3)
sp_db3 <- sp_db3 %>% 
  mutate(pr_ab = 0) %>% 
  mutate(pr_ab = ifelse(filt < quantile(filt)[2], 1, 0))
sp_db3 <- data.frame(species='sp3', sp_db3)

spp <- dplyr::bind_rows(sp_db, sp_db2, sp_db3)

head(spp)
tail(spp)
spp %>% group_by(species, pr_ab) %>% count()

# plot(somevar[[1]], col=gray.colors(10))
# points(sp_db3[,2:3], col=c('red','blue')[sp_db3$pr_ab+1], pch=19)
rm(sp_db); rm(sp_db2); rm(sp_db3)

# save(spp,
     # file = file.path(here::here('Data', 'spp.RData')))

source('./R/block_partition.R')
load('./Data/spp.RData')
load('./Data/somevar.RData')

res(somevar)*1000
ini <- Sys.time()
plot(somevar)
part <- block_partition_pa(
  env_layer = somevar,
  occ_data = spp,
  sp = 'species',
  x = 'x',
  y = 'y',
  pr_ab = 'pr_ab',
  max_res_mult = 500,
  num_grids = 30,
  n_part = 4,
  cores = 3,
  save_part_raster = TRUE,
  dir_save = getwd() #Write the directory path to save results
)
Sys.time()-ini

dim(part)

part %>% group_by(sp, pr_ab, partition) %>% count



##%######################################################%##
#                                                          #
#### Second version to test block_partition_pa function ####
#                                                          #
##%######################################################%##

# Packages
require(raster)
require(dplyr)

# Import function
source('https://raw.githubusercontent.com/sjevelazco/spatial_sp_traits/main/R/block_partition.R')


### Two occurrences databases
# 1-Virtual species 
URL <- "https://github.com/sjevelazco/spatial_sp_traits/blob/main/Data/spp.RData?raw=true"
load(url(URL))
spp

# 2-Real species 
spp2 <- vroom::vroom("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/1-PresencesOnlySDMs/BlockCrossValidation//test_sp_pres_abs.gz")
spp2 %>% dplyr::select(ends_with('coords')) %>% plot
spp2 %>% dplyr::select(CADE27) %>% table # "I assume that the column with presences-absences data is CADE27"
spp2 %>% tail
spp2 <- tibble(spp='Abies magnifica', spp2) # Function require a columns with species names I will assume that it is speceis database for Abies magnifica :)


### Load environmental variables
# URL <- "https://github.com/sjevelazco/spatial_sp_traits/raw/main/Data/somevar.RData"
# td <- tempdir()
# download.file(URL, file.path(td, 'somevar.RData'))
# load(file.path(td, 'somevar.RData')) # some env. variables. created by Brooke for CFP

env.stack <- brick('C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/Predictors/BCM1981_2010_CFP_Stack.grd') # this database is 
# in Brook's folder https://drive.google.com/drive/u/1/folders/1WDW4ryl1N29aK6Jmf3VK8XDi3NeYw_mG


# Previews plot to check if raster data and occurrence overlap
plot(env.stack[[1]])
points(spp2 %>% dplyr::select(ends_with('coords')), cex=0.5, pch=19)

env.stack <- homogenize_na(env.stack)
plot(env.stack[[1]])

part2 <- block_partition_pa(
  env_layer = env.stack,
  occ_data = spp, 
  sp = 'species', 
  x = 'x', 
  y = 'y', 
  pr_ab = 'pr_ab', 
  max_res_mult = 500, 
  num_grids = 30, 
  n_part = 2, 
  cores = 3, 
  save_part_raster = TRUE, 
  dir_save = getwd()
)

part2 %>% 
  dplyr::group_by(sp, partition, pr_ab) %>% 
  count 
