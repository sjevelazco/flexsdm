# data

#' A data set containing presences and absences of three virtual species
#'
#' @format A tibble with 1150 rows and 3 variables:
#' \describe{
#'   \item{species}{virtual species names}
#'   \item{x}{longitude of species occurrences}
#'   \item{y}{latitude of species occurrences}
#'   \item{pr_ab}{presences and absences denoted by 1 and 0 respectively}
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
#'   \item{x y}{columns with geographical coordinates}
#'   \item{from column aet to landform}{columns with values of environmental variables at coordinate locations}
#'   ...
#' }
#' @examples
#' \dontrun{
#' require(dplyr)
#' data("backg")
#' backg
#' }
"backg"

#' A data set containing localities and environmental condition of an Abies (fir tree) species in California, USA
#'
#' @format A tibble object with 5000 rows and 10 variables:
#' \describe{
#'   \item{ID}{presences and absences records ID}
#'   \item{pr_ab}{presence and absences denoted by 1 and 0 respectively}
#'   \item{x y}{columns with coordinates in Albers Equal Area Conic coordinate system}
#'   \item{from column aet to landform}{columns with values for environmental variables at each locality}
#'   ...
#' }
#' @examples
#' \dontrun{
#' require(dplyr)
#' data("abies")
#' abies
#' }
"abies"

#' A data set containing localities of Hesperocyparis stephensonii species in California, USA
#'
#' @format A tibble object with 14 rows and 4 variables:
#' \describe{
#'   \item{ID}{presences records ID}
#'   \item{x y}{columns with coordinates in Albers Equal Area Conic coordinate system}
#'   \item{pr_ab}{presence denoted by 1}
#'   ...
#' }
#' @examples
#' \dontrun{
#' require(dplyr)
#' data("hespero")
#' hespero
#' }
"hespero"

#' A data set containing presences of palms species from Southern Brazil
#'
#' @description
#' This data set contains presences of 11 palm species from Southern Brazil sourced by Calambás-Trochez et al. (2021).
#'
#' @format A tibble with 327 rows and 3 variables:
#' \describe{
#'   \item{species}{species names}
#'   \item{x}{longitude of species occurrences}
#'   \item{y}{latitude of species occurrences}
#'   ...
#' }
#'
#' @references
#' Calambás-Trochez, L.F., Velazco, S.J.E., Hoffmann, P.M., Brum, F.T., Carlucci, M.B., 2021. Climate and land-use changes coupled with low coverage of protected areas threaten palm species in South Brazilian grasslands. Perspectives in Ecology and Conservation 9. https://doi.org/10.1016/j.pecon.2021.03.010
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' data("palms")
#'
#' }
"palms"
