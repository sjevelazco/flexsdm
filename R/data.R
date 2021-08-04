# data

#' A data set containing presences and absences of tree virtual species
#'
#' @format A tibble with 1150 rows and 3 variables:
#' \describe{
#'   \item{species}{virtual species names}
#'   \item{x}{longitude of species occurrences}
#'   \item{y}{latitude of species occurrences}
#'   \item{pr_ab}{presence and absences denoted by 1 and 0 respectively}
#'   ...
#' }
#' @examples
#' \dontrun{
#' require(dplyr)
#' data("spp")
#' spp
#' }
"spp"

#' A data set containing environmental conditions of background points
#'
#' @format A tibble object with 5000 rows and 10 variables:
#' \describe{
#'   \item{pr_ab}{background point denoted by 0}
#'   \item{x y}{columns with coordinates}
#'   \item{from column aet to landform}{columns with environmental variables}
#'   ...
#' }
#' @examples
#' \dontrun{
#' require(dplyr)
#' data("backg")
#' backg
#' }
"backg"

#' A data set containing environmental condition of an Abies species
#'
#' @format A tibble object with 5000 rows and 10 variables:
#' \describe{
#'   \item{ID}{presence and absences records ID}
#'   \item{pr_ab}{presence and absences denoted by 1 and 0 respectively}
#'   \item{x y}{columns with coordinates}
#'   \item{from column aet to landform}{columns with environmental variables}
#'   ...
#' }
#' @examples
#' \dontrun{
#' require(dplyr)
#' data("abies")
#' abies
#' }
"abies"
