##%######################################################%##
#                                                          #
####                Plot Data Cleaning                  ####
#                                                          #
##%######################################################%##

# Written by Brooke Rose (with parts borrowed from Santiago)

wd <- list()
# commonly used paths in my working directory
wd$data   <- "/Users/Brooke Rose/Google Drive/Franklin_grant/project/data/"
wd$output <- "/Users/Brooke Rose/Google Drive/Franklin_grant/project/output/"
wd$scripts <- "/Users/Brooke Rose/Google Drive/Franklin_grant/project/scripts/"
wd$software <- "/Users/Brooke Rose/Google Drive/Franklin_grant/project/software/"

# packages and functions
require(here)
require(sf)
require(readxl)
require(dplyr)
require(janitor)
require(readr)
require(stringr)
require(textclean)
require(raster)
require(ggplot)
require(tidyverse)
require(fuzzySim)

# Study species list
study_sp <- data.table::fread(paste0(wd$data, 'plots/study_species.csv')) %>% tibble()
study_sp <- study_sp %>% janitor::clean_names() %>% na.omit()
study_sp_names <- study_sp$species_name
study_sp_names <- gsub(x = study_sp_names, pattern = "\\.", replacement = " ") 

# California Floristic province shapefile
cfp <- st_read(paste0(wd$data,'shapefiles/CFP/CFP_GIS_California.shp'))
cfp_trans <- st_transform(cfp, crs="+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

# Environmental rasters
env_stack <- raster::brick(paste0(wd$output, 'predictors/env_stack/BCM1981_2010_CA_CFP.grd'))

##%######################################################%##
#                                                          #
####                Prepare plant codes                 ####
#                                                          #
##%######################################################%##
plant_codes <- data.table::fread(paste0(wd$data, 'Thorne_plots/California USDA Plants codes.csv'))
plant_codes <- plant_codes %>% janitor::clean_names()
plant_codes <- tibble(plant_codes, scientific_name=NA)
for( i in 1:nrow(plant_codes)){
  print(i)
  plant_codes$scientific_name[i] <- 
    flora::remove.authors(plant_codes$scientific_name_with_author[i])
}

plant_codes
tail(plant_codes)

data.table::fwrite(plant_codes, paste0(wd$data, 'plots/California USDA Plants codes cleaned.gz'))


##%######################################################%##
#                                                          #
####                 Thorne Relevees                    ####
#                                                          #
##%######################################################%##

all_releves <-
  list.files(paste0(wd$data, 'Thorne_plots'),
             pattern = 'Relevees_',
             full.names = T)  %>%
  sapply(., vroom::vroom)

all_releves <- dplyr::bind_rows(all_releves) 

colnames(all_releves)[1:7] <- 
  all_releves[1:7] %>% janitor::clean_names() %>% names

### plot of raw data (clearly issues with coordinate systems)
all_releves %>% ggplot(aes(lon, lat)) + 
  geom_hex(bins=200) + 
  coord_equal() + theme_minimal()

### replacing NA's with 0's in sp columns
filt <- names(all_releves)[1:7]
pres_col <- names(all_releves)[-c(1:7)]

###  removing plots with missing coordinates
all_releves <- all_releves %>%
  dplyr::filter(complete.cases(all_releves %>% dplyr::select(lon, lat)))

all_releves <-
  all_releves %>% 
  dplyr::mutate(dplyr::across(pres_col, ~ tidyr::replace_na(.x, 0)))

all_releves <- all_releves %>% arrange(new_id) #sort rows based on new_id

### Remove  identical rows coordinates and summarize information at new_id level 
all_releves %>% dplyr::count(new_id) %>% arrange(desc(n)) #almost all new_id is repeated

df2 <- # unique data for each plot
  all_releves %>% dplyr::select(-c('total_of_old_plot_id', pres_col)) %>% unique

df3 <- all_releves %>% # summarizing species records by plot 
  dplyr::select(c('new_id', c('total_of_old_plot_id', pres_col))) %>%
  group_by(new_id) %>%
  summarise(across(c('total_of_old_plot_id', pres_col), ~ sum(.x)))

all_releves2 <- dplyr::left_join(df2, df3, by = "new_id")
#dim(all_releves)
cat("Dimensions of clean data frame: ", dim(all_releves2))
all_releves <- all_releves2
rm(all_releves2)

# nrow(unique(all_releves)) == nrow(all_releves) # all rows are unique. however there are plot with same coordinates 
# which not does not necessarily represent duplicate coordinates with different species composition

### Exploring plot without any species
rowfilt <- all_releves[pres_col] %>% colSums()

#### Plant codes
# Database with accepted names
filt <- colnames(all_releves)[-c(1:7)] %>%  sort
filt#vector with species symbol

plant_codes <- vroom::vroom(paste0(wd$data, 'plots/California USDA Plants codes cleaned.gz'))

plant_codes_accepted <-
  plant_codes %>% 
  dplyr::filter(symbol %in% filt) %>% 
  dplyr::filter(is.na(synonym_symbol))

# Database with accepted and synonyms names for those columns names with synonyms
plant_codes_synonym <-
  plant_codes %>% dplyr::filter(synonym_symbol %in% filt)
plant_codes_synonym <-
  plant_codes %>% dplyr::filter(symbol %in% plant_codes_synonym$symbol)
plant_codes_synonym <-
  plant_codes_synonym %>% filter(synonym_symbol %in% filt |
                                   is.na(synonym_symbol))

# Update columns names
filt <- names(all_releves)[1:7]
pres_col <- names(all_releves)[-c(1:7)]

for(i in 1:length(pres_col)){
  message('processing species name ', i)
  sp_symbol <- pres_col[i]
  if(sp_symbol %in% plant_codes_accepted$symbol) {
    sp_names <- plant_codes_accepted %>%
      dplyr::filter(symbol %in% sp_symbol) %>%
      pull(scientific_name)
    pres_col2 <- pres_col[i]
    names(pres_col2) <- sp_names
    colnames(all_releves)[colnames(all_releves)==pres_col2] <- names(pres_col2)
    # all_releves <- dplyr::rename(all_releves, pres_col2)
    
  } else if (sp_symbol %in% plant_codes_synonym$synonym_symbol) {
    message('Updating species names')
    sp_symbol2 <- plant_codes_synonym %>%
      dplyr::filter(synonym_symbol %in% sp_symbol) %>% #note here is used synonym_symbol
      pull(symbol) #here is extracted the accepted names symbol
    sp_names <- plant_codes_synonym %>%
      dplyr::filter(symbol %in% sp_symbol2) %>%
      dplyr::filter(is.na(synonym_symbol)) %>%
      pull(scientific_name)
    
    pres_col2 <- pres_col[i]
    names(pres_col2) <- sp_names
    colnames(all_releves)[colnames(all_releves)==pres_col2] <- names(pres_col2)
    # all_releves <- dplyr::rename(all_releves, pres_col2)
    
  } else {
    message('This species code has no species names in USDA database')
  }
}

# One repeated column name
sum(duplicated(colnames(all_releves)))
colnames(all_releves)[duplicated(colnames(all_releves))] %>% table() %>% sort
repeatednames <- colnames(all_releves)[duplicated(colnames(all_releves))] %>% unique

tempdf <- data.frame(matrix(NA, nrow = nrow(all_releves), ncol=length(repeatednames)))
colnames(tempdf) <- repeatednames
for(i in 1:length(repeatednames)){
  message(i)
  cnames <- colnames(all_releves)
  fitl <- which(grepl(repeatednames[i], cnames))
  temp <- all_releves[fitl] %>% 
    apply(.,1, sum, na.rm = T) 
  temp <- ifelse(is.na(temp),0,temp)
  temp <- ifelse(temp>1,1,temp)
  tempdf[,i] <- temp
}
colSums(tempdf) %>% plot
all_releves <- all_releves[,-c(which(duplicated(names(all_releves))))]

for(i in 1:ncol(tempdf)){
  message(i)
  all_releves[names(tempdf[i])] <- tempdf[,i]
}

data.table::fwrite(all_releves, paste0(wd$data, 'plots/Thorne_releve_sp_names_raw.gz'))

#### Correcting Thorne geographic data to match environmental data ####
all_releves <- data.table::fread(paste0(wd$data, 'plots/Thorne_releve_sp_names_raw.gz')) %>% tibble()


# Thorne releve sources (all in different coordinate systems, lots of issues)
unique(all_releves$sarea)
all_releves %>% ggplot(aes(x=lon, y=lat, fill = sarea)) + 
  geom_hex(bins=200) + 
  coord_equal() + theme_minimal()

# HARR data are in longitude, latitude
harr <- all_releves %>% 
  filter(sarea == "HARR") %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4008, remove = FALSE) %>%
  st_transform(crs(env_stack))

# USFSN, PRGG, POTT, CCBLM, and SFHREL are in UTM zone 10
utm10_plots <- all_releves %>%
  filter(sarea == 'USFSN' | sarea == "PRGG" | sarea == "POTT" | sarea == "CCBLM" | sarea == "SFHREL") %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 32610, remove = FALSE) %>%
  st_transform(crs(env_stack))

# VTM, BMC, and LPBO are in California Albers NAD27
albers_plots <- all_releves %>% 
  filter(sarea == "VTM" | sarea == 'BMC' | sarea == 'LPBO') %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 3309, remove = FALSE) %>%
  st_transform(crs(env_stack))

# SEKINRI data are missing decimal places, so I multiply latitude and longitude by 1000 (Okayed by Jim)
sekinri <- all_releves %>%
  filter(sarea == 'SEKINRI') %>%
  mutate(lat = lat * 1000,
         lon = lon * 1000) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 32611, remove = FALSE) %>%
  st_transform(crs(env_stack))

# SAMO, SEKIAA, SEIRA, SEKIVM and JOTR are in UTM zone 11
utm_zone11 <- all_releves %>% 
  filter(sarea == 'SAMO'| sarea == 'SEKIAA'| sarea == 'SEKIRA' | sarea == 'SEKIVM' | sarea == 'JOTR') %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 32611, remove = FALSE) %>%
  st_transform(crs(env_stack))

# Thorne releve data with corrected and transformed coordinate systems
all_releves2 <- dplyr::bind_rows(harr, utm10_plots, albers_plots, sekinri, utm_zone11)

data.table::fwrite(all_releves2, paste0(wd$data, 'plots/Thorne_releve_sp_names_correct_coordinates.gz'))

# Removing VTM plots (collected in the 1930s), adding consistent coordinates (x and y albers), and retaining only plots in the CFP
final_releves <- all_releves2 %>% 
  filter(sarea != 'VTM') %>%
  dplyr::mutate(x_albers = sf::st_coordinates(.)[,1],
                y_albers = sf::st_coordinates(.)[,2]) %>%
  st_intersection(st_make_valid(cfp_trans)) %>%
  dplyr::mutate(survey = 'relevee',
                source = 'Jim Thorne') 

# Final Thorne Relevees
ggplot() +
  geom_sf(data = cfp_trans) +
  geom_sf(data = final_releves, aes(color = sarea)) +
  labs(title = paste('Thorne relevee plots:', nrow(final_releves)))

ggsave(paste0(wd$data, 'plots/Thorne_relevees_cfp_plot_map.jpeg'))

data.table::fwrite(final_releves, paste0(wd$data, 'plots/Thorne_releve_sp_names_novtm_clean.gz'))

##%######################################################%##
#                                                          #
####                 Thorne Rapids                      ####
#                                                          #
##%######################################################%##


all_rapids <-
  list.files(paste0(wd$data, 'Thorne_plots'),
             pattern = 'tbl_R',
             full.names = T)  %>%
  sapply(., vroom::vroom)


### Cleaning data
all_rapids <- dplyr::bind_rows(all_rapids) 

filt <- all_rapids %>% dplyr::select(-c(1:8)) %>% colnames %>% sort
filt #vector with species symbol

# Database with accepted names
plant_codes_accepted <-
  plant_codes %>% dplyr::filter(symbol %in% filt)

# Database with accepted and synonyms names for those columns names with synonyms
plant_codes_synonym <-
  plant_codes %>% dplyr::filter(synonym_symbol %in% filt)
plant_codes_synonym <-
  plant_codes %>% dplyr::filter(symbol %in% plant_codes_synonym$symbol)
plant_codes_synonym <-
  plant_codes_synonym %>% filter(synonym_symbol %in% filt |
                                   is.na(synonym_symbol))


### Clean columns names 
filt <- all_rapids %>% dplyr::select(c(1:8)) %>% names
newnames <- all_rapids[1, filt] %>% 
  janitor::clean_names() %>% names
colnames(all_rapids)[1:8] <- newnames

### Remove plots without coordinates
cat("Total number of rapids plots: ", nrow(all_rapids))
all_rapids <- all_rapids %>%
  dplyr::filter(complete.cases(all_rapids %>% dplyr::select(x_final, y_final)))
cat("Number of rapids plots with coordinates: ", nrow(all_rapids))

### plots in geographic space
all_rapids %>% ggplot(aes(x_final, y_final)) + 
  geom_hex(bins=200) + 
  coord_equal() + theme_minimal()


### Fill with zero in those columns with species names
filt <- all_rapids %>% dplyr::select(c(1:8)) %>% names
pres_col <- names(all_rapids)[-c(1:8)]

all_rapids <-
  all_rapids %>% 
  dplyr::mutate(dplyr::across(pres_col, ~ tidyr::replace_na(.x, 0)))

all_rapids <- all_rapids %>% arrange(new_id) #sort rows based on new_id

### Remove  identical rows coordinates and summarize information at new_id level 
all_rapids %>% dplyr::count(new_id) %>% arrange(desc(n)) #almost all new_id is repeated

df2 <- # unique data for each plot
  all_rapids %>% dplyr::select(-c('total', pres_col)) %>% unique

df3 <- all_rapids %>% # summarizing species records by plot 
  dplyr::select(c('new_id', c('total', pres_col))) %>%
  group_by(new_id) %>%
  summarise(across(c('total', pres_col), ~ sum(.x)))

all_rapids2 <- dplyr::left_join(df2, df3, by = "new_id")
#dim(all_rapids)
cat("Dimensions of clean data frame: ", dim(all_rapids2))
all_rapids <- all_rapids2
rm(all_rapids2)


### Exploring plot without any species
rowfilt <- all_rapids[pres_col] %>% colSums()

### Replacing species codes with species names
for( i in 1:length(pres_col)){
  message('processing species name ', i)
  
  sp_symbol <- pres_col[i]
  if(sp_symbol %in% plant_codes_accepted$symbol) {
    sp_names <- plant_codes_accepted %>%
      dplyr::filter(symbol %in% sp_symbol) %>%
      pull(scientific_name_with_author)
    sp_names <- flora::remove.authors(sp_names)
    pres_col2 <- pres_col[i]
    names(pres_col2) <- sp_names
    all_rapids <- dplyr::rename(all_rapids, pres_col2)
    
  } else if (sp_symbol %in% plant_codes_synonym$synonym_symbol) {
    message('Updating species names')
    sp_symbol2 <- plant_codes_synonym %>%
      dplyr::filter(synonym_symbol %in% sp_symbol) %>% #note here is used synonym_symbol
      pull(symbol) #here is extracted the accepted names symbol
    sp_names <- plant_codes_synonym %>%
      dplyr::filter(symbol %in% sp_symbol2) %>%
      dplyr::filter(is.na(synonym_symbol)) %>%
      pull(scientific_name_with_author)
    
    sp_names <- flora::remove.authors(sp_names)
    pres_col2 <- pres_col[i]
    names(pres_col2) <- sp_names
    all_rapids <- dplyr::rename(all_rapids, pres_col2)
    
  } else {
    message('This species code has no species names in USDA database')
  }
}

data.table::fwrite(all_rapids, paste0(wd$data, 'plots/Thorne_rapids_sp_names_raw.gz'))

#### Correcting Thorne geographic data to match environmental data ####
# all Thorne rapids are in NAD_1927_California_Teale_albers (EPSG = 3309)

all_rapids <- data.table::fread(paste0(wd$data, 'plots/Thorne_rapids_sp_names_raw.gz')) %>% tibble()

all_rapids2 <- all_rapids %>% 
  st_as_sf(coords = c('x_final', 'y_final'), crs = 3309, remove = FALSE) %>%
  st_transform(crs(env_stack)) %>%
  dplyr::mutate(x_albers = sf::st_coordinates(.)[,1],
                y_albers = sf::st_coordinates(.)[,2]) %>%
  st_intersection(st_make_valid(cfp_trans)) %>%
  dplyr::mutate(survey = 'rapid',
                source = 'Jim Thorne')

ggplot() +
  geom_sf(data = cfp_trans) +
  geom_sf(data = all_rapids2, aes(color = sarea)) +
  labs(title = paste('Thorne rapids:', nrow(all_rapids2)))

ggsave(paste0(wd$data, 'plots/Thorne_rapids_cfp_plot_map.jpeg'))

data.table::fwrite(all_rapids2, paste0(wd$data, 'plots/Thorne_rapids_sp_names_clean.gz'))

##%######################################################%##
#                                                          #
####               Cal Fish and Wildlife                ####
#                                                          #
##%######################################################%##

# Survey plot information
survey_plots <- st_read(paste0(wd$data, 'BIOS/Cal_survey_points/Cal_survey_plots.shp'))
# Waypoint as character (for joining)
survey_plots$WyptID <- as.character(survey_plots$WyptID)
survey_plots <- survey_plots %>% janitor::clean_names()

# Plant list of species in Cal Fish and Wildlife plots
plant_list <- vroom::vroom(paste0(wd$data, 'BIOS/Cal_survey_points/SurveyPlants.csv'))
plant_list <- plant_list %>% janitor::clean_names()

# Joining survey information with plant list
survey_df <- left_join(survey_plots, plant_list, by = c("wypt_id" = "waypoint_id")) %>%
  janitor::clean_names() %>%
  dplyr::mutate(x_coord = sf::st_coordinates(.)[, 1],
                y_coord = sf::st_coordinates(.)[, 2]) %>%
  st_set_geometry(NULL) 

# for joining to new presence absence data
plots_info <- survey_df %>% 
  dplyr::select(wypt_id, x_coord, y_coord, survey_type, survey_date, location_n, soil_text, micro_topo, macro_topo, project_cod) %>% unique

# converting survey_df to presence/absence dataset
pres_abs <- splist2presabs(
  survey_df,
  sites.col = "wypt_id",
  sp.col = "curr_plants_symbol",
  keep.n = FALSE
)

cfw_data <- left_join(plots_info, pres_abs, by = 'wypt_id')


#### Plant codes
plant_codes <- vroom::vroom(paste0(wd$data, 'plots/California USDA Plants codes cleaned.gz'))

# Rename some columns to match more data with plant_codes
cfw_data <- data.frame(cfw_data)
cnames <- colnames(cfw_data)
cnames2 <- gsub("X2", "", gsub("X2JM", "", cnames))
cnames2%in%plant_codes$symbol  %>% sum
cnames%in%plant_codes$symbol %>% sum
colnames(cfw_data) <- cnames2
any(duplicated(colnames(cfw_data))) # there area columns with the same names 
sum(duplicated(colnames(cfw_data)))

cfw_data %>% dplyr::select(c(1:8))
filt <- colnames(cfw_data)[-c(1:8)] %>%  sort
filt#vector with species symbol

# Database with accepted names
plant_codes_accepted <-
  plant_codes %>% 
  dplyr::filter(symbol %in% filt) %>% 
  dplyr::filter(is.na(synonym_symbol))

# Database with accepted and synonyms names for those columns names with synonyms
plant_codes_synonym <-
  plant_codes %>% dplyr::filter(synonym_symbol %in% filt)
plant_codes_synonym <-
  plant_codes %>% dplyr::filter(symbol %in% plant_codes_synonym$symbol)
plant_codes_synonym <-
  plant_codes_synonym %>% filter(synonym_symbol %in% filt |
                                   is.na(synonym_symbol))

### Select species 
filt <- names(cfw_data)[1:8]
pres_col <- names(cfw_data)[-c(1:8)]

for(i in 1:length(pres_col)){
  message('processing species name ', i)
  sp_symbol <- pres_col[i]
  if(sp_symbol %in% plant_codes_accepted$symbol) {
    sp_names <- plant_codes_accepted %>%
      dplyr::filter(symbol %in% sp_symbol) %>%
      pull(scientific_name)
    pres_col2 <- pres_col[i]
    names(pres_col2) <- sp_names
    colnames(cfw_data)[colnames(cfw_data)==pres_col2] <- names(pres_col2)
    # cfw_data <- dplyr::rename(cfw_data, pres_col2)
    
  } else if (sp_symbol %in% plant_codes_synonym$synonym_symbol) {
    message('Updating species names')
    sp_symbol2 <- plant_codes_synonym %>%
      dplyr::filter(synonym_symbol %in% sp_symbol) %>% #note here is used synonym_symbol
      pull(symbol) #here is extracted the accepted names symbol
    sp_names <- plant_codes_synonym %>%
      dplyr::filter(symbol %in% sp_symbol2) %>%
      dplyr::filter(is.na(synonym_symbol)) %>%
      pull(scientific_name)
    
    pres_col2 <- pres_col[i]
    names(pres_col2) <- sp_names
    colnames(cfw_data)[colnames(cfw_data)==pres_col2] <- 
      names(pres_col2)
    # cfw_data <- dplyr::rename(cfw_data, pres_col2)
    
  } else {
    message('This species code has no species names in USDA database')
  }
}

# Lot of columns with same names
sum(duplicated(colnames(cfw_data)))
colnames(cfw_data)[duplicated(colnames(cfw_data))] %>% table() %>% sort
repeatednames <- colnames(cfw_data)[duplicated(colnames(cfw_data))] %>% unique

tempdf <- data.frame(matrix(NA, nrow = nrow(cfw_data), ncol=length(repeatednames)))
colnames(tempdf) <- repeatednames
for(i in 1:length(repeatednames)){
  message(i)
  cnames <- colnames(cfw_data)
  fitl <- which(grepl(repeatednames[i], cnames))
  temp <- cfw_data[fitl] %>% 
    apply(.,1, sum, na.rm = T) 
  temp <- ifelse(is.na(temp),0,temp)
  temp <- ifelse(temp>1,1,temp)
  tempdf[,i] <- temp
}
colSums(tempdf) %>% plot
cfw_data <- cfw_data[,-c(which(duplicated(names(cfw_data))))]

for(i in 1:ncol(tempdf)){
  message(i)
  cfw_data[names(tempdf[i])] <- tempdf[,i]
}

# adding survey variable
cfw_data <- cfw_data %>%
  dplyr::mutate(
    survey = ifelse(
      survey_type == 'Transect' | survey_type == 'Releve' |
        survey_type == 'Multivisit Releve' |
        survey_type == 'releve' |
        survey_type == 'Multivisit Transect' |
        survey_type == 'Multi-visit releve' |
        survey_type == 'multi-visit transect' |
        survey_type == 'Relevee',
      'relevee',
      'rapid'
    ),
    source = 'Cal Fish and Wildlife'
  ) %>%
  st_as_sf(coords = c('x_coord', 'y_coord'), crs = 3310, remove = FALSE) 

ggplot() +
  geom_sf(data = cfp_trans) +
  geom_sf(data = cfw_data, aes(color = project_cod)) +
  labs(title = paste('CFW plots', nrow(cfw_data)))

ggsave(paste0(wd$data, 'plots/cfw_full_plot_map.jpeg'))

data.table::fwrite(cfw_data, paste0(wd$data, 'plots/cfw_sp_names_raw.gz'))

# Including only CFP data and transforming data to match environmental data
cfw_data <- data.table::fread(paste0(wd$data, 'plots/cfw_sp_names_raw.gz')) %>% tibble()

cfw_cfp <- cfw_data %>%
  st_as_sf(coords = c('x_coord', 'y_coord'), crs = 3310, remove = FALSE)  %>%
  st_transform(crs(env_stack)) %>%
  dplyr::mutate(x_albers = sf::st_coordinates(.)[,1],
                y_albers = sf::st_coordinates(.)[,2]) %>%
  st_intersection(st_make_valid(cfp_trans)) %>%
  relocate(x_albers:geometry, .before = Mirabilis.linearis)

ggplot() +
  geom_sf(data = cfp_trans) +
  geom_sf(data = cfw_cfp, aes(color = project_cod)) +
  labs(title = paste('CFW plots in CFP', nrow(cfw_cfp)))

ggsave(paste0(wd$data, 'plots/cfw_cfp_plot_map.jpeg'))

data.table::fwrite(cfw_cfp, paste0(wd$data, 'plots/cfw_sp_names_cfp_cleaned.gz'))

##%######################################################%##
#                                                          #
####                     CalFlora                       ####
#                                                          #
##%######################################################%##

all_calflora <-
  list.files(paste0(wd$data, 'CalFlora'),
             full.names = T)  %>%
  lapply(., vroom::vroom) 

all_calflora <- dplyr::bind_rows(all_calflora)

all_calflora <- all_calflora %>% janitor::clean_names() # 41052 plots

all_calflora2 <- dplyr::bind_rows(all_calflora) %>%
  dplyr::mutate(source = 'CalFlora',
         survey = 'CalFlora') %>%
  dplyr::filter(location_quality == 'high' | location_quality == 'medium') %>% # removing low quality locations
  subset(date > as.Date("1980-01-01"))

nrow(all_calflora) - nrow(all_calflora2) # How many plots with low or NA location quality & before 1980?

# removing observations in zoos and botanical gardens
non_wild <- c('botanic', 'Botanic', 'botanical', 'Botanical', 'Zoo', 'zoo')

non_wild_df <- all_calflora2[grepl(paste(non_wild,collapse="|"), 
                          all_calflora2$location_description),]

all_calflora3 <- all_calflora2[!grepl(paste(non_wild,collapse="|"), 
                                    all_calflora2$location_description),]

all_calflora <- all_calflora3

pres_abs <- splist2presabs(
  all_calflora,
  sites.col = "id",
  sp.col = "taxon",
  keep.n = FALSE
)

cf_pa <- left_join(all_calflora, pres_abs, by ='id')

cf_sf <- st_as_sf(cf_pa, coords = c("longitude", "latitude"), remove = FALSE) %>%
  st_set_crs(4008) %>%
  st_transform(crs(env_stack)) %>%
  st_intersection(st_make_valid(cfp_trans)) %>%
  dplyr::mutate(x_albers = sf::st_coordinates(.)[,1],
                y_albers = sf::st_coordinates(.)[,2]) %>%
  relocate(AREA_SQ_MI:y_albers, .before = Abies.bracteata)


ggplot() +
  geom_sf(data = cfp_trans) +
  geom_sf(data = cf_sf) +
  labs(title = paste('CalFlora plots in CFP', nrow(cf_sf)))

ggsave(paste0(wd$data, 'plots/calflora_cfp_plot_map.jpeg'))

data.table::fwrite(cf_sf, paste0(wd$data, 'plots/calflora_sp_names_cfp_cleaned.gz'))

##%######################################################%##
#                                                          #
####               Integrate all databases              ####
#                                                          #
##%######################################################%##
require(vroom)
require(ggplot2)
require(dplyr)

# Cleaned plot data
thorne_releves <- data.table::fread(paste0(wd$data, 'plots/Thorne_releve_sp_names_novtm_clean.gz')) %>% tibble 
names(thorne_releves) <- gsub(x = names(thorne_releves), pattern = "\\.", replacement = " ") 

thorne_rapids <- data.table::fread(paste0(wd$data, 'plots/Thorne_rapids_sp_names_clean.gz')) %>% tibble
names(thorne_rapids) <- gsub(x = names(thorne_rapids), pattern = "\\.", replacement = " ") 

cfw_plots <- data.table::fread(paste0(wd$data, 'plots/cfw_sp_names_cfp_cleaned.gz')) %>% tibble

cf_plots <- data.table::fread(paste0(wd$data, 'plots/calflora_sp_names_cfp_cleaned.gz')) %>% tibble
names(cf_plots) <- gsub(x = names(cf_plots), pattern = "\\.", replacement = " ") 


# all data
all_data <- dplyr::bind_rows(cf_plots, cfw_plots, thorne_rapids, thorne_releves) %>%
  st_as_sf(coords = c('x_albers', 'y_albers'), crs = crs(env_stack), remove = FALSE)


ggplot() +
  geom_sf(data = cfp_trans, fill = "antiquewhite") +
  geom_sf(data = all_data, aes(color = survey), alpha = .5) +
  labs(title = paste("All plots in the CFP:", paste(nrow(all_data)))) +
  theme(title = element_text(color = 'black', size = 12, family = "serif"),
        text = element_text(size = 12, family = "serif", color = 'black')) +
  facet_wrap(~source)

ggsave(paste0(wd$data, 'plots/all_plots_map.jpeg'))

data.table::fwrite(all_data, paste0(wd$data, 'plots/all_data_cfp_sp_names.gz'))


##%######################################################%##
#                                                          #
####                  Model data frame                  ####
#                                                          #
##%######################################################%##

# columns with study species names
study_cols <- unique(grep(paste(study_sp_names,collapse="|"), 
                          names(all_data), value=TRUE))

# data frame with study species and plot identifiers
study_df <- all_data %>%
  dplyr::select(
    id,
    wypt_id,
    new_id,
    x_albers,
    y_albers,
    date,
    source,
    JEP_REG,
    survey,
    geometry,
    all_of(study_cols)
  )

sp_extract <- raster::extract(x = env_stack, y = study_df, df = TRUE, sp = TRUE, cellnumbers = TRUE)

sp_data <- st_as_sf(sp_extract, remove = FALSE) %>%
  st_set_geometry(NULL)  %>%
  relocate(cells:landform, .before = Abies.bracteata)

sp_data[21:214][is.na(sp_data[21:214])] <- 0

# removing observations that have NA for these environmental variables (they are generated by BCM and do not have 100% full coverage)
sp_env <- sp_data[ , c("aet", "cwd", "tmin", "ppt_djf", "ppt_jja", "landform", "pH", "awc", "depth", "percent_clay")]  

sp_env_full <- sp_data[complete.cases(sp_env), ] 

nrow(sp_env_full)

for( x in 1:length(study_sp_names)){
  colnames(sp_env_full)[grepl(study_sp_names[x], colnames(sp_env_full), ignore.case = TRUE)] <- study_sp_names[x]
}

data.table::fwrite(sp_env_full, paste0(wd$data, 'plots/sp_pres_abs_clean.gz'))


