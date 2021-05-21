#' Extract environmental data based on x and y coordinates
#'
#' @param data data.frame. Database with presences, presence-absence, or pseudo-absence, records with x and y coordinates
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param env_layer SpatRaster. Raster with environmental variables.
#' @param predictors character. Vector with the variable names of predictor variables
#' Usage predictors = c("aet", "cwd", "tmin")
#' @param filter_na logical. If filter_na = TRUE (default), the rows with NA values for any of the
#' environmental variables are removed from the returned tibble.
#'
#' @return
#'
#' A tibble that returns the original data base of presence, presence-absence, or pseudo-absence location records with additional columns
#' for the extracted environmental variables at each xy location from the SpatRast object 'env_layer'
#'
#' @export
#'
#' @examples
sdm_extract <- function(data, x, y, predictors, env_layer, filter_na = TRUE) {

  # spatial data frame
  sp_data <-
    terra::vect(data,
      geom = c(x, y),
      crs = crs(env_layer)
    )

  # extract environmental data at xy locations, if filter_na = FALSE, does not remove rows with NAs
  if (filter_na == FALSE) {
    extract_data <- data.frame(
      data,
      terra::extract(env_layer,
        sp_data,
        cells = TRUE
      ) %>%
        dplyr::select(
          cell,
          dplyr::all_of(predictors)
        )
    )

    extract_data <- as_tibble(extract_data)
  }
  # extract environmental data at xy locations, removes rows with NAs for any environmental variable
  else {
    extract_data <- data.frame(
      data,
      terra::extract(env_layer,
        sp_data,
        cells = TRUE
      ) %>%
        dplyr::select(
          cell,
          dplyr::all_of(predictors)
        )
    )

    complete_vec <- complete.cases(extract_data[, predictors])
    extract_data <- extract_data[complete_vec, ]
  }

  return(extract_data)
}
