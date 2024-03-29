% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{spp}
\alias{spp}
\title{A data set containing presences and absences of three virtual species}
\format{
A tibble with 1150 rows and 3 variables:
\describe{
  \item{species}{virtual species names}
  \item{x}{longitude of species occurrences}
  \item{y}{latitude of species occurrences}
  \item{pr_ab}{presences and absences denoted by 1 and 0 respectively}
  ...
}
}
\usage{
spp
}
\description{
A data set containing presences and absences of three virtual species
}
\examples{
\dontrun{
require(dplyr)
data("spp")
spp
}
}
\keyword{datasets}
