% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_gau.R
\name{fit_gau}
\alias{fit_gau}
\title{Fit and validate Gaussian Process models}
\usage{
fit_gau(
  data,
  response,
  predictors,
  predictors_f = NULL,
  background = NULL,
  partition,
  thr = NULL
)
}
\arguments{
\item{data}{data.frame. Database with response (0,1) and predictors values.}

\item{response}{character. Column name with species absence-presence data (0,1).}

\item{predictors}{character. Vector with the column names of quantitative
predictor variables (i.e. continuous variables).
Usage predictors = c("aet", "cwd", "tmin")}

\item{predictors_f}{character. Vector with the column names of qualitative
predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")}

\item{background}{data.frame. Database with response column only with 0 and predictors variables. All
column names must be consistent with data}

\item{partition}{character. Column name with training and validation partition groups.}

\item{thr}{character. Threshold used to get binary suitability values (i.e. 0,1), useful for threshold-dependent performance metrics. It is possible to use more than one threshold type. A vector must be provided for this argument. The following threshold criteria are available:
\itemize{
  \item lpt: The highest threshold at which there is no omission.
  \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
  \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the highest (aka threshold that maximizes the TSS).
  \item max_jaccard: The threshold at which the Jaccard index is the highest.
  \item max_sorensen: The threshold at which the Sorensen index is the highest.
  \item max_fpb: The threshold at which FPB (F-measure on presence-background data) is highest.
  \item sensitivity: Threshold based on a specified sensitivity value.
  Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to sensitivity value. If a sensitivity value is not specified, the default used is 0.9.
  }
If more than one threshold type is used they must be concatenated, e.g., thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity', sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function will use all threshold criteria if none is specified.}
}
\value{
A list object with:
\itemize{
\item model: A "graf" class object. This object can be used for predicting.
\item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
\item performance: Performance metrics (see \code{\link{sdm_eval}}).
Threshold dependent metrics are calculated based on the threshold criteria specified in the argument.
\item data_ens: Predicted suitability for each test partition. This database is used in \code{\link{fit_ensemble}}
}
}
\description{
Fit and validate Gaussian Process models
}
\examples{
\dontrun{
data("abies")

# Using k-fold partition method
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 3)
)
abies2

bg <- abies2
bg$pr_ab <- 0


gaup_t1 <- fit_gau(
  data = abies2,
  response = "pr_ab",
  predictors = c("aet", "ppt_jja", "pH", "awc", "depth"),
  predictors_f = c("landform"),
  partition = ".part",
  background = bg,
  thr = c("max_sens_spec", "equal_sens_spec", "max_sorensen")
)

gaup_t1$model
gaup_t1$predictors
gaup_t1$performance
gaup_t1$data_ens

# Using bootstrap partition method only with presence-absence
abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "boot", replicates = 5, proportion = 0.7)
)
abies2

gaup_t2 <- fit_gau(
  data = abies2,
  response = "pr_ab",
  predictors = c("ppt_jja", "pH", "awc"),
  predictors_f = c("landform"),
  partition = ".part",
  thr = c(type = c("lpt", "max_sens_spec", "sensitivity"), sens = "0.8")
)
gaup_t2
}
}
\seealso{
\code{\link{fit_gam}}, \code{\link{fit_glm}}, \code{\link{fit_gbm}},
\code{\link{fit_max}}, \code{\link{fit_net}}, \code{\link{fit_raf}}, and \code{\link{fit_svm}}.
}
