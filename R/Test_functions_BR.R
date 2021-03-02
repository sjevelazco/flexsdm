##%######################################################%##
#                                                          #
####                 Test some functions                ####
#                                                          #
##%######################################################%##

require(raster)
require(dplyr)
require(sf)

# working directory on Franklin Dell Desktop
wd <- list()

# commonly used paths in my working directory
wd$data   <- "/Users/Brooke Rose/Google Drive/Franklin_grant/project/data/"
wd$output <- "/Users/Brooke Rose/Google Drive/Franklin_grant/project/output/"
wd$scripts <- "/Users/Brooke Rose/Google Drive/Franklin_grant/project/scripts/"
wd$software <- "/Users/Brooke Rose/Google Drive/Franklin_grant/project/software/"

# working directory on Quercus Dell Desktop
wd <- list()

# commonly used paths in my working directory
wd$data   <- "/Users/Brooke/Google Drive/Franklin_grant/project/data/"
wd$output <- "/Users/Brooke/Google Drive/Franklin_grant/project/output/"
wd$scripts <- "/Users/Brooke/Google Drive/Franklin_grant/project/scripts/"
wd$software <- "/Users/Brooke/Google Drive/Franklin_grant/project/software/"


env_stack <- raster::brick(paste0(wd$output, 'predictors/env_stack/BCM1981_2010_CA_CFP.grd')) #environmental predictors

sp_data <- data.table::fread(paste0(wd$output, 'pres_abs_data.gz')) %>% tibble # species data

cfp <- st_read(paste0(wd$data,'shapefiles/CFP/CFP_GIS_California.shp')) %>%
  st_transform(crs = crs(env_stack)) # California portion of CFP

study_sp <- data.table::fread(paste0(wd$data, 'plots/study_species.csv')) %>% tibble()
study_sp <- study_sp %>% janitor::clean_names() %>% na.omit()
pres_abs_sp <- study_sp %>% dplyr::filter(model_type == 'presence-absence')
study_sp_names <- pres_abs_sp$species_name
study_sp_names <- gsub(x = study_sp_names, pattern = "\\.", replacement = " ")

##%######################################################%##
#                                                          #
####                 Map Plot Function                  ####
#                                                          #
##%######################################################%##

source('./R/pretty_map_fun.R')

test_pa_map <- pretty_map_fun(plot_area = cfp,
                       occ_data = sp_data %>% dplyr::filter(species == 'Abies magnifica'),
                       x = 'x_albers',
                       y = 'y_albers',
                       epsg_code = 3310,
                       fill_att = "pr_ab",
                       title = "Abies magnifica",
                       subtitle = paste(
                         nrow(occ_data %>% filter(pr_ab == 1)),
                         "presences;",
                         nrow(occ_data %>% filter(pr_ab == 0)),
                         "absences"
                       ))


test_source_map <- pretty_map_fun(plot_area = cfp,
                          occ_data = sp_data %>% dplyr::filter(species == 'Abies magnifica'),
                          x = 'x_albers',
                          y = 'y_albers',
                          epsg_code = 3310,
                          fill_att = "source",
                          title = "Abies magnifica data sources",
                          subtitle = paste(
                            nrow(occ_data %>% filter(source == "Jim Thorne")),
                            "Thorne plots;",
                            nrow(occ_data %>% filter(source == "Cal Fish and Wildlife")),
                            "CFW plots;",
                            nrow(occ_data %>% filter(source == "CalFlora")),
                            "CalFlora plots")
)

# testing with no titles, does not work
no_titles_map <- pretty_map_fun(plot_area = cfp,
                          occ_data = sp_data %>% dplyr::filter(species == 'Abies magnifica'),
                          x = 'x_albers',
                          y = 'y_albers',
                          epsg_code = 3310,
                          fill_att = "pr_ab")


##%######################################################%##
#                                                          #
####           Outlier Detection Function               ####
#                                                          #
##%######################################################%##

# Outlier Function
source("https://raw.githubusercontent.com/sjevelazco/spatial_sp_traits/main/R/env_outliers.R")

pa_data <- list()

# P/A datasets (by species)
for (i in 1:length(study_sp_names)) {
  pa_data[[i]] <- data.table::fread(paste0(wd$output,
                                           study_sp_names[i],
                                           '/pres_abs_for_modeling.gz')) %>% tibble()
  print(study_sp_names[i])
}

outliers <- list()

for(i in 1:length(pa_data)){
outliers[[i]] <- env_outliers(da = pa_data[[i]],
                        species = 'species',
                        x = 'x_albers',
                        y = 'y_albers',
                        pr_ab = 'pr_ab',
                        envr = env_stack,
                        id = 'ID')
}

# R Session timed-out March 2, 2021 (I will come back to this problem later )
