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


##%######################################################%##
#                                                          #
####            Test tuning Brooke function             ####
#                                                          #
##%######################################################%##
load('./Data/spp_pres_psabs.RData')
env <- raster::brick("C:/Users/SVelazco/Documents/Santiago/1-PresencesOnlySDMs/3-Variables/Predictors/BCM1981_2010_CA_CFP.grd")

require(raster)
require(dplyr)
require(data.table)

envext <- raster::extract(env, spp_pres_psabs %>% dplyr::select(longitude_m, latitude_m))
env_preds <- colnames(envext)
df <- data.frame(spp_pres_psabs, envext)
head(df)
df <- split(df,df$pr_ab)
set.seed(123)
train0 <- dplyr::slice_sample(df$`0`, prop = 0.7)
set.seed(123)
train1 <- dplyr::slice_sample(df$`1`, prop = 0.7)
test0 <- anti_join(df$`0`, train0)
test1 <- anti_join(df$`1`, train1)

calib <- bind_rows(train0, train1)
eval <- bind_rows(test0, test1)
df <- bind_rows(df)
df_clean <- df %>% dplyr::select(pr_ab, all_of(env_preds))

calib <- calib[c('pr_ab', env_preds)]
eval <- eval[c('pr_ab', env_preds)]
dim(calib)
dim(eval)

source("./R/sdm_function_BRv2.R")

sdms(
  df = split_data[[i]],
  eval = eval_data[[i]],
  calib = calib_data[[i]],
  pr_ab = "pr_ab",
  env_preds = env_preds,
  sp_area = sp_areas[[i]],
  pred_raster = env_stack,
  species_name = test_species[[i]],
  dir_save = dir_save[[i]], cores = 8
)

require(parallel)
detectCores()


##%######################################################%##
#                                                          #
####              Test  evaluate function               ####
#                                                          #
##%######################################################%##
require(dismo)
require(dplyr)
source("./R/boyce.R")
source("./R/enm_eval.R")

set.seed(0)
p <- rnorm(50, mean=0.7, sd=0.3) %>% abs()
p[p>1] <- 1
p[p<0] <- 0

set.seed(0)
a <- rnorm(50, mean=0.4, sd=0.4) %>% abs()
a[a>1] <- 1
a[a<0] <- 0

set.seed(0)
backg <- rnorm(1000, mean=0.4, sd=0.4) %>% abs()
backg[backg>1] <- 1
backg[backg<0] <- 0


# Use function without threshold specification 
e <- enm_eval(p, a)
e$performance
e$threshold
e$threshold_table

# Different ways to use thr argument
enm_eval(p, a, thr=c(type=c('MAX_KAPPA')))
enm_eval(p, a, thr=c(type=c('MAX_KAPPA')), bg=backg)
enm_eval(p, a, thr=c(type=c('LPT', 'MAX_TSS', 'MAX_JACCARD')))
enm_eval(p, a, thr=c(type=c('LPT', 'MAX_TSS', 'SENSITIVITY'))) # wrong way to SENSITIVITY threshold
enm_eval(p, a, thr=c(type=c('LPT', 'MAX_TSS', 'SENSITIVITY'), sens='0.8')) # correct way to use SENSITIVITY threshold

# Use of bg argument (it will only be used for calculating BOYCE index)
enm_eval(p, a, thr=c(type=c('MAX_TSS')))[[1]]
enm_eval(p, a, thr=c(type=c('MAX_TSS')), bg=backg)[[1]]
# I the case it is needed use background for calculate all other metric background values can be used in "a" argument 
enm_eval(p, backg, thr=c(type=c('MAX_TSS')))[[1]]


##%######################################################%##
#                                                          #
####                    Test mx_tune                    ####
#                                                          #
##%######################################################%##
# require(prefixer)
# prefixer::prefixer()
# prefixer::import_from("tune_mx")

# mammals_data <- read.csv("https://raw.githubusercontent.com/vdicolab/hsdm/master/data/tabular/species/mammals_and_bioclim_table.csv", row.names=1)
# ## Create the RF model
# mammals_data$VulpesVulpes %>% table
# mammals_data2 <- mammals_data %>% rename(pr_ab=VulpesVulpes, x=X_WGS84, y=Y_WGS84)
# mammals_data2 <- mammals_data2 %>% dplyr::select(pr_ab, starts_with('bio'))
# readr:::write_tsv(mammals_data2, './Data/mammals_data2.txt')
ex_data <- readr::read_tsv('./Data/mammals_data2.txt')

# Background points (only for this example, bg will be assumed as the union of presences and absences records)
backgd <- ex_data

# Absences
na <- ex_data %>% dplyr::filter(pr_ab == 0) %>% nrow
# Presences
np <- ex_data %>% dplyr::filter(pr_ab == 1) %>% nrow
# Background
bg <- backgd %>% nrow

N <- 2
na <- sample(rep(1:N, length.out=na))
np <- sample(rep(1:N, length.out=np))
bg <- sample(rep(1:N, length.out=bg))
ex_data <- ex_data %>% dplyr::mutate(Partition=c(na, np))
backgd <- backgd %>% dplyr::mutate(Partition=bg)

gridtest <-
  expand.grid(regmult = seq(0.1, 3, 0.5),
              classes = c("l", "lq", "lqh", "lqhp", "lqhpt"))

  # if (np < 10) {
  #   classes <- "l"
  # } else if (np < 15){
  #   classes <- "lq"
  # } else if (np < 80) {
  #   classes <- "lqh"
  # }
  #   

source('./R/tune_mx.R')
source("./R/boyce.R")
source("./R/enm_eval.R")

r <-
  tune_mx(
    data = ex_data,
    background = backgd,
    response = "pr_ab",
    predictors = c("bio3", 'bio4', 'bio7', 'bio11', 'bio12'),
    predictors_f = NULL,
    partition = "Partition",
    grid = gridtest,
    thr = "MAX_TSS",
    metric = 'TSS',
    clamp = FALSE,
    pred_type = "cloglog" 
  )

r$model %>% plot(type="cloglog")
r$tune_performance
r$tune_performance$BOYCE_mean

r$best_hyperparameter
r$threshold

require(ggplot2)
ggplot(r$tune_performance, aes(regmult, TSS_mean, col=classes)) + 
  geom_errorbar(aes(ymin=TSS_mean-TSS_sd, ymax=TSS_mean+TSS_sd), width=0.1)+
  geom_point() + 
  geom_line()




##%######################################################%##
#                                                          #
####                   test tune_svm                    ####
#                                                          #
##%######################################################%##
require(maxnet)
data("bradypus")
names(bradypus)
# Absences
na <- bradypus %>% dplyr::filter(presence == 0) %>% nrow
# Presences
np <- bradypus %>% dplyr::filter(presence == 1) %>% nrow

N <- 10
na <- sample(rep(1:N, length.out=na))
np <- sample(rep(1:N, length.out=np))
bradypus <- bradypus %>% dplyr::mutate(Partition=c(na, np))
colnames(bradypus)


source("./R/tune_svm.R")
source("./R/boyce.R")
source("./R/enm_eval.R")
require(ggplot2)

tune_grid <-
  expand.grid(C = c(1, 2, 4, 8, 16, 20),
              sigma = c(0.001, 0.01, 0.1, 0.2, 0.3, 0.4))
svmt <-
  tune_svm(
    data = bradypus,
    response = "presence",
    predictors = c("cld6190_ann", "dtr6190_ann", "h_dem", "pre6190_ann", "tmn6190_ann", "vap6190_ann"),
    predictors_f = c("ecoreg"),
    partition = "Partition",
    grid = tune_grid,
    thr = "MAX_TSS",
    metric = 'TSS',
  )

svmt$model
svmt$tune_performance
svmt$best_hyper_performance
svmt$best_hyper
svmt$selected_threshold
svmt$threshold_table

ggplot(svmt$tune_performance, aes(factor(sigma), 
                                  TSS_mean, col=factor(C))) + 
  geom_errorbar(aes(ymin=TSS_mean-TSS_sd, ymax=TSS_mean+TSS_sd), width=0.1)+
  geom_point() + 
  geom_line(data = svmt$tune_performance, aes(as.numeric(factor(sigma)), 
                                              TSS_mean, col=factor(C)))

ggplot(svmt$tune_performance, aes(factor(sigma), 
                                  AUC_mean, col=factor(C))) + 
  # geom_errorbar(aes(ymin=AUC_mean-AUC_sd, ymax=AUC_mean+AUC_sd), width=0.1)+
  geom_point() + 
  geom_line(data = svmt$tune_performance, aes(as.numeric(factor(sigma)), 
                                              AUC_mean, col=factor(C)))

ggplot(svmt$tune_performance, aes(factor(sigma), 
                                  TSS_mean, col=factor(C))) + 
  geom_point() + 
  geom_line(data = svmt$tune_performance, aes(as.numeric(factor(sigma)), 
                                              TSS_mean, col=factor(C)))


##%######################################################%##
#                                                          #
####                   test tune_rf                    ####
#                                                          #
##%######################################################%##
source("./R/tune_rf.R")
source("./R/boyce.R")
source("./R/enm_eval.R")
require(maxnet)

tune_grid <-
  expand.grid(mtry = seq(1, 7, 1)) # The maximum mtry must be equal to total number of predictors

rf_t <-
  tune_rf(
    data = bradypus,
    response = "presence",
    predictors = c("cld6190_ann", "dtr6190_ann", "h_dem", "pre6190_ann", "tmn6190_ann", "vap6190_ann"),
    predictors_f = c("ecoreg"),
    partition = "Partition",
    grid = tune_grid,
    thr = "MAX_TSS",
    metric = 'TSS',
  )

rf_t$model %>% plot
rf_t$tune_performance %>% arrange(desc(TSS_mean)) 
rf_t$best_hyper_performance
rf_t$best_hyper
rf_t$selected_threshold
rf_t$threshold_table

require(ggplot2)
ggplot(rf_t$tune_performance, aes(factor(mtry), 
                                  TSS_mean)) + 
  geom_errorbar(aes(ymin=TSS_mean-TSS_sd, ymax=TSS_mean+TSS_sd), width=0.1)+
  geom_point() + 
  geom_line(data = rf_t$tune_performance, aes(as.numeric(factor(mtry)), 
                                              TSS_mean))

##%######################################################%##
#                                                          #
####                   test tune_nnet                   ####
#                                                          #
##%######################################################%##
source("./R/tune_nnet.R")
source("./R/boyce.R")
source("./R/enm_eval.R")

tune_grid <-
  expand.grid(size = c(2, 4, 6, 8, 10), 
              decay = c(0.001, 0.05, 0.1,1, 3, 4, 5, 10))

nnet_t <-
  tune_nnet(
    data = bradypus,
    response = "presence",
    predictors = c("cld6190_ann", "dtr6190_ann", "h_dem", "pre6190_ann", "tmn6190_ann", "vap6190_ann"),
    predictors_f = c("ecoreg"),
    partition = "Partition",
    grid = tune_grid,
    thr = "MAX_TSS",
    metric = 'TSS',
  )

nnet_t$model
nnet_t$tune_performance %>% arrange(desc(TSS_mean)) 
nnet_t$best_hyper_performance
nnet_t$best_hyper
nnet_t$selected_threshold
nnet_t$threshold_table

require(ggplot2)

ggplot(nnet_t$tune_performance, aes(factor(decay), 
                                  TSS_mean, col=factor(size))) + 
  # geom_errorbar(aes(ymin=TSS_mean-TSS_sd, ymax=TSS_mean+TSS_sd), width=0.1)+
  geom_point() + 
  geom_line(data = nnet_t$tune_performance, aes(as.numeric(factor(decay)), 
                                              TSS_mean, col=factor(size)))

ggplot(nnet_t$tune_performance, aes(factor(decay), 
                                    JACCARD_mean, col=factor(size))) + 
  # geom_errorbar(aes(ymin=TSS_mean-TSS_sd, ymax=TSS_mean+TSS_sd), width=0.1)+
  geom_point() + 
  geom_line(data = nnet_t$tune_performance, aes(as.numeric(factor(decay)), 
                                                JACCARD_mean, col=factor(size)))

ggplot(nnet_t$tune_performance, aes(factor(decay), 
                                    AUC_mean, col=factor(size))) + 
  # geom_errorbar(aes(ymin=TSS_mean-TSS_sd, ymax=TSS_mean+TSS_sd), width=0.1)+
  geom_point() + 
  geom_line(data = nnet_t$tune_performance, aes(as.numeric(factor(decay)), 
                                                AUC_mean, col=factor(size)))
