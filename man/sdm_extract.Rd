\name{sdm_extract}
\alias{sdm_extract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SDM Extract
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sdm_extract(data, x, y, predictors, env_layer, filter_na = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{predictors}{
%%     ~~Describe \code{predictors} here~~
}
  \item{env_layer}{
%%     ~~Describe \code{env_layer} here~~
}
  \item{filter_na}{
%%     ~~Describe \code{filter_na} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, x, y, predictors, env_layer, filter_na = TRUE)
{
    sp_data <- terra::vect(data, geom = c(x, y), crs = crs(env_layer))
    if (filter_na == FALSE) {
        extract_data <- data.frame(data, terra::extract(env_layer,
            sp_data, cells = TRUE) \%>\% dplyr::select(cell, dplyr::all_of(predictors)))
        extract_data <- as_tibble(extract_data)
    }
    else {
        extract_data <- data.frame(data, terra::extract(env_layer,
            sp_data, cells = TRUE) \%>\% dplyr::select(cell, dplyr::all_of(predictors)))
        complete_vec <- complete.cases(extract_data[, predictors])
        extract_data <- extract_data[complete_vec, ]
    }
    return(extract_data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
