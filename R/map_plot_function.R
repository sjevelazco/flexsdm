## Written by Brooke Rose


#' Title
#'
#' @param plot_area 
#' @param occ_data 
#' @param x 
#' @param y 
#' @param epsg_code 
#' @param fill_att 
#' @param title 
#' @param subtitle 
#'
#' @return
#' @export
#'
#' @examples
pretty_map <- function(plot_area, # shapefile/polygon for main extent of plot
                       occ_data, # plot location data (with a presence absence column)
                       x, # x column 
                       y, # y column
                       epsg_code, # EPSG code for the coordinate system you want to use
                       fill_att, # attribute in shapefile  with character for plotting point color (presence/absence, data source, number of outlier methods, etc.)
                       title, # title of plot
                       subtitle){ # subtitle of plot 
  require(sf)
  require(ggplot2)
  require(tidyverse)
  require(ggsn)
  
  if(data.class(plot_area) != 'sf'){
    plot_area <- st_as_sf(plot_area)
  }
  if(data.class(occ_data) != 'sf'){
    occ_data <- st_as_sf(occ_data, coords = c(x, y), crs = epsg_code, remove = FALSE)
  }
  
  ggplot() +
    geom_sf(data = plot_area,
            color = "black",
            fill = "papayawhip") +
    geom_sf(data = occ_data, aes(color = as.character(!!sym(fill_att)))) +
    labs(
      title = paste(title),
      subtitle = paste(subtitle),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      legend.direction = "vertical",
      legend.position = "right",
      legend.key = element_blank(),
      legend.background = element_rect(fill = "white"),
      text = element_text(
        size = 15,
        family = "serif",
        color = 'black'
      ),
      axis.text.x = element_text(color = 'black', size = 15),
      axis.text.y = element_text(color = 'black', size = 15)
    ) +
    north(plot_area, location = 'topright', scale = .11) +
    ggsn::scalebar(
      plot_area,
      dist = 150,
      dist_unit = "km",
      transform = FALSE,
      st.size = 3,
      family = 'serif',
      location = 'bottomleft'
    )
}