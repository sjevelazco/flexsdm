% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tune_gbm.R
\name{tune_gbm}
\alias{tune_gbm}
\title{Fit and validate Generalized Boosted Regression models with exploration of
hyper-parameters that optimize performance}
\usage{
tune_gbm(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition,
  grid = NULL,
  thr = NULL,
  metric = "TSS",
  n_cores = 1
)
}
\arguments{
\item{data}{data.frame. Database with response (0,1) and predictors values.}

\item{response}{character. Column name with species absence-presence data (0,1).}

\item{predictors}{character. Vector with the column names of quantitative predictor variables
(i.e. continuous variables). Usage predictors = c("aet", "cwd", "tmin")}

\item{predictors_f}{character. Vector with the column names of qualitative predictor variables
(i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")}

\item{fit_formula}{formula. A formula object with response and predictor
variables (e.g. formula(pr_ab ~ aet + ppt_jja + pH + awc + depth + landform)).
Note that the variables used here must be consistent with those used in
response, predictors, and predictors_f arguments. Default is NULL.}

\item{partition}{character. Column name with training and validation partition groups.}

\item{grid}{data.frame. A data frame object with algorithm hyper-parameter values to be tested.
It Is recommended to generate this data.frame with the grid() function. Hyper-parameters needed
for tuning are 'n.trees', 'shrinkage', and 'n.minobsinnode'.}

\item{thr}{character. Threshold used to get binary suitability values (i.e. 0,1) needed for
threshold-dependent performance metrics. It is possible to use more than one threshold type.
Provide a vector for this argument. The following threshold types are available:
\itemize{
  \item lpt: The highest threshold at which there is no omission.
  \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
  \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the
  highest (aka threshold that maximizes the TSS).
  \item max_jaccard: The threshold at which the Jaccard index is the highest.
  \item max_sorensen: The threshold at which the Sorensen index is highest.
  \item max_fpb: The threshold at which FPB (F-measure on presence-background data) is highest.
  \item sensitivity: Threshold based on a specified sensitivity value.
  Usage thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens' refers to
  sensitivity value. If no sensitivity value is specified, the default used is 0.9
  }
If more than one threshold type is used they must be concatenate, e.g.,
thr=c('lpt', 'max_sens_spec', 'max_jaccard'), or thr=c('lpt', 'max_sens_spec', 'sensitivity',
sens='0.8'), or thr=c('lpt', 'max_sens_spec', 'sensitivity'). Function will use all threshold
types if no threshold is specified.}

\item{metric}{character. Performance metric used for selecting the best combination of
hyper-parameter values. The following metrics can be used: SORENSEN, JACCARD, FPB, TSS,
KAPPA, AUC, and BOYCE. TSS is used as the default.}

\item{n_cores}{numeric. Number of cores use for parallelization. Default 1}
}
\value{
A list object with:
\itemize{
\item model: A "gbm" class object from gbm package. This object can be used for predicting.
\item predictors: A tibble with quantitative (c column names) and qualitative (f column names)
variables use for modeling.
\item performance: Hyper-parameter values and performance metric (see \code{\link{sdm_eval}})
for the best hyper-parameter combination.
\item hyper_performance: Performance metric (see \code{\link{sdm_eval}}) for each combination of
the hyper-parameters.
\item data_ens: Predicted suitability for each test partition based on the best model. This
database is used in \code{\link{fit_ensemble}}
}
}
\description{
Fit and validate Generalized Boosted Regression models with exploration of
hyper-parameters that optimize performance
}
\examples{
\dontrun{
data(abies)
abies

# Partition the data with the k-fold method

abies2 <- part_random(
  data = abies,
  pr_ab = "pr_ab",
  method = c(method = "kfold", folds = 5)
)

# pr_ab is the name of the column with species presence and absences (i.e. the response variable)
# from aet to landform are the predictors variables (landform is a qualitative variable)

# Hyper-parameter values for tuning
tune_grid <-
  expand.grid(
    n.trees = c(20, 50, 100),
    shrinkage = c(0.1, 0.5, 1),
    n.minobsinnode = c(1, 3, 5, 7, 9)
  )

gbm_t <-
  tune_gbm(
    data = abies2,
    response = "pr_ab",
    predictors = c(
      "aet", "cwd", "tmin", "ppt_djf", "ppt_jja",
      "ppt_jja", "pH", "awc", "depth"
    ),
    predictors_f = c("landform"),
    partition = ".part",
    grid = tune_grid,
    thr = "max_sens_spec",
    metric = "TSS",
    n_cores = 1
  )

# Outputs
gbm_t$model
gbm_t$predictors
gbm_t$performance
gbm_t$data_ens
gbm_t$hyper_performance

# Graphical exploration of performance of each hyper-parameter setting
require(ggplot2)
pg <- position_dodge(width = 0.5)
ggplot(gbm_t$hyper_performance, aes(factor(n.minobsinnode),
  TSS_mean,
  col = factor(shrinkage)
)) +
  geom_errorbar(aes(ymin = TSS_mean - TSS_sd, ymax = TSS_mean + TSS_sd),
    width = 0.2, position = pg
  ) +
  geom_point(position = pg) +
  geom_line(
    data = gbm_t$tune_performance,
    aes(as.numeric(factor(n.minobsinnode)),
      TSS_mean,
      col = factor(shrinkage)
    ), position = pg
  ) +
  facet_wrap(. ~ n.trees) +
  theme(legend.position = "bottom")
}

}
\seealso{
\code{\link{tune_max}}, \code{\link{tune_net}}, \code{\link{tune_raf}},
and \code{\link{tune_svm}}.
}
