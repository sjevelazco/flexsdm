##%######################################################%##
#                                                          #
####                 Test some function                 ####
#                                                          #
##%######################################################%##

require(raster)
require(dplyr)
require(sf)

##%######################################################%##
#                                                          #
####                Basic layer database                ####
#                                                          #
##%######################################################%##

# r_base <- "C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/StudyAreaRaster.tif" %>%
#   raster::raster()
# dir.create("./Data")
# save(r_base, file=file.path(here::here('Data', 'basic_layer.RData')), compress = 'bzip2')


# somevar <- "C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/CFP_environmental_stack.tif" %>%
  # raster::stack()
# somevar <- somevar[[c(1:4)]]
# save(somevar, file=file.path(here::here('Data', 'somevar.RData')), compress = 'bzip2')
# load(here::here('Data', 'somevar.RData'))

##%######################################################%##
#                                                          #
####       Create a presences absences database         ####
#                                                          #
##%######################################################%##

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



# Save Rdata a real occurrence database
spp_real <- vroom::vroom("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/1-PresencesOnlySDMs/BlockCrossValidation/study_plots.gz")
spp_real <- rename(spp_real, pr_ab=ABMA, x=x_tran, y=y_tran)
spp_real <- spp_real %>% dplyr::select(-c(new_id:geometry))

# save(spp_real, file = file.path(here::here('Data', 'spp_real.RData')))
rm(spp_real)



##%######################################################%##
#                                                          #
####                  Teste function 1                  ####
#                                                          #
##%######################################################%##
source('./R/block_partition_pa.R')
load('./Data/spp.RData')
load('./Data/somevar.RData')
env.stack <- brick('C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/Predictors/BCM1981_2010_CFP_Stack_2v.grd')

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
  n_part = 2,
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

# Import function and occurrence database 
source('./R/block_partition_pa.R')
load('./Data/spp_real.RData')

### Load environmental variables
env_stack <- brick('C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/Predictors/BCM1981_2010_CFP_Stack_2v.grd')

# Filtering record placed in the same cell 
ncell <- raster::cellFromXY(env_stack, spp_real %>% dplyr::select(x, y) %>% data.frame)
spp_real <- tibble(ncell, spp_real) %>% arrange(ncell, desc(pr_ab))
spp_real <- spp_real %>% dplyr::filter(!duplicated(ncell))  
spp_real %>% count(pr_ab)

ini <- Sys.time()
part2 <- block_partition_pa(
  env_layer = env_stack,
  occ_data = spp_real, 
  sp = 'species', 
  x = 'x', 
  y = 'y', 
  pr_ab = 'pr_ab',
  min_res_mult = 10, 
  max_res_mult = 500, 
  num_grids = 30, 
  n_part = 2, 
  cores = 3, 
  save_part_raster = TRUE, 
  dir_save = getwd()
)
Sys.time()-ini

part2 %>% 
  dplyr::group_by(sp, partition, pr_ab) %>% 
  count 
part2 %>% 
  dplyr::group_by(sp, partition) %>% 
  count 


##%######################################################%##
#                                                          #
####          Test function second Brooke bug          ####
#                                                          #
##%######################################################%##
tibble(pa_data[[1]])

source('https://raw.githubusercontent.com/sjevelazco/spatial_sp_traits/main/R/block_partition_pa.R')
spat_part <- block_partition_pa(
  env_layer = env_stack,
  occ_data = pa_data[[1]], #
  sp = 'species',
  x = 'x_tran',
  y = 'y_tran',
  pr_ab = 'sp_pa',
  min_res_mult=2, # I recommend set this parameter with higher values maybe 3 or 4 
  max_res_mult = 500, # you can use a higher value here
  num_grids = 5, # when really use this function for your data set higher number here maybe > 30 
  n_part = 2,
  cores = 3, # don't need multiple cores for single species
  save_part_raster = TRUE,
  dir_save = getwd()
)


# A way to visualize how looks like a raster setting a given value in the max_res_mult parameter 
source('https://raw.githubusercontent.com/sjevelazco/spatial_sp_traits/main/R/plot_max_res.R')
plot_max_res(env_stack[[1]], max_res_mult = 800)


##%######################################################%##
#                                                          #
####           Detection of outliers records            ####
####       based on environmental characteristics       ####
#                                                          #
##%######################################################%##

# Pres pseudo-absences database
URL <- "https://raw.githubusercontent.com/sjevelazco/spatial_sp_traits/main/Data/spp_pres_psabs.RData?raw=true"
load(url(URL))
spp_pres_psabs

# environmental variables
URL <- "https://github.com/sjevelazco/spatial_sp_traits/raw/main/Data/somevar.RData"
td <- tempdir()
download.file(URL, file.path(td, 'somevar.RData'))
load(file.path(td, 'somevar.RData')) # some env. vari
somevar


# Function
source("https://raw.githubusercontent.com/sjevelazco/spatial_sp_traits/main/R/env_outliers.R")


somevar # raster stack

spp <- spp_pres_psabs$search_name %>% unique
table(spp_pres_psabs$search_name)
spp_pres_psabs2 <- spp_pres_psabs %>% dplyr::filter(search_name%in%spp[1:3]) # filter database to get smaller database

out <- env_outliers(da = spp_pres_psabs2,
                    species = 'search_name',
                    x = 'longitude_m',
                    y = 'latitude_m',
                    pr_ab = 'pr_ab',
                    envr = somevar,
                    id = 'IDr')

head(out)

nrow(spp_pres_psabs2)
spp_pres_psabs3 <- left_join(spp_pres_psabs2, out, by=c('IDr' = 'id')) 


##%######################################################%##
#                                                          #
####            Test  env_filtering function            ####
#                                                          #
##%######################################################%##
require(dplyr)
require(raster)
source('./R/env_filtering.R')
load('./Data/somevar.RData')
load('./Data/spp.RData')
# somevar %>% plot
spp <- spp %>% dplyr::mutate(ID=as.character(1:nrow(spp))) %>% 
  dplyr::rename(lon=x, lat=y) %>% tibble()
somevar <- raster::brick(somevar)

# env_variables <- raster::brick('C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/Predictors/BCM1981_2010_CA_CFP.grd') # if you read or transform a raster stack to brick, the raster::extract() function will work considerably faster.
table(spp$species)
occ_filtered <-
  env_filtering(
    da = spp %>% dplyr::filter(species == "sp3") %>% data.frame,
    x = 'lon',
    y = 'lat',
    id = 'ID',
    variables = somevar,
    nbins = 20,
    cores = 3,
    plot = T
  )
dim(occ_filtered)
dim(spp)



