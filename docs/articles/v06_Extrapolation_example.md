# flexsdm: Tools to explore extrapolation in SDMs

## Introduction

Many SDM applications require model extrapolation, e.g., predictions
beyond the range of the data set used to fit the model. For example,
models often must extrapolate when predicting habitat suitability under
novel environmental conditions induced by climate change or predicting
the spread of an invasive species outside of its native range based on
the species-environment relationship observed in its native range.

In *flexsdm*, we offer a new approach (known as
[**Shape**](https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.06992))
for evaluating the extrapolation and truncating spatial predictions
based on the degree of extrapolation measured. Shape is a model-agnostic
approach for calculating the degree of extrapolation for a given
projection data point by its multivariate distance to the nearest
training data point – capturing the often complex shape of data within
environmental space. These distances are then relativized by a factor
that reflects the dispersion of the training data in environmental
space. As implemented in *flexsdm*, the Shape approach also incorporates
an adjustable threshold to allow for binary discrimination between
acceptable and unacceptable degrees of extrapolation, based on the
user’s needs and applications. For more information about Shape metric,
we recommend reading the article [Velazco et al.,
2023](https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.06992).

In this vignette, we will walk through how to evaluate model
extrapolation for *Hesperocyparis stephensonii* (Cuyamaca cypress), a
conifer tree species that is endemic to southern California. This
species is listed as Critically Endangered by the IUCN and has an
extremely restricted distribution, as it is only found in the headwaters
of King Creek in San Diego County.

Note: this tutorial follows generally the same workflow as the vignette
for modeling the distribution of a rare species using an ensemble of
small models (ESM). However, instead of constructing ESMs, we will
evaluate model extrapolation if we were to predict our models to the
extent of the California Floristic Province (CFP).

## Data

For our models, we will use four environmental variables that influence
plant distributions in California: available evapotranspiration (aet),
climatic water deficit (cwd), maximum temperature of the warmest month
(tmx), and minimum temperature of the coldest month (tmn). Our
occurrence data include 21 geo-referenced observations downloaded from
the online database Calflora.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`flexsdm`](https://sjevelazco.github.io/flexsdm/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`terra`](https://rspatial.org/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`patchwork`](https://patchwork.data-imaginist.com)`)`\
\
`# environmental data`\
`somevar`` ``<-`` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"external/somevar.tif"``, package ``=`` ``"flexsdm"``)`\
`somevar`` ``<-`` ``terra``::`[`rast`](https://rspatial.github.io/terra/reference/rast.html)`(``somevar``)`\
[`names`](https://rspatial.github.io/terra/reference/names.html)`(``somevar``)`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"cwd"``, ``"tmn"``, ``"aet"``, ``"ppt_jja"``)`\
\
`# species occurence data (presence-only)`\
[`data`](https://rdrr.io/r/utils/data.html)`(``hespero``)`\
`hespero`` ``<-`` ``hespero`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` ``dplyr``::`[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``-``id``)`\
\
`# California ecoregions`\
`regions`` ``<-`` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"external/regions.tif"``, package ``=`` ``"flexsdm"``)`\
`regions`` ``<-`` ``terra``::`[`rast`](https://rspatial.github.io/terra/reference/rast.html)`(``regions``)`\
`regions`` ``<-`` ``terra``::`[`as.polygons`](https://rspatial.github.io/terra/reference/as.polygons.html)`(``regions``)`\
`sp_region`` ``<-`` ``terra``::`[`subset`](https://rspatial.github.io/terra/reference/subset.html)`(``regions``, ``regions``$``category`` ``==`` ``"SCR"``)`` ``# ecoregion where *Hesperocyparis stephensonii* is found`\
\
`# visualize the species occurrences`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(`\
`  ``sp_region``,`\
`  col ``=`` ``"gray80"``,`\
`  legend ``=`` ``FALSE``,`\
`  axes ``=`` ``FALSE``,`\
`  main ``=`` ``"Hesperocyparis stephensonii occurrences"`\
`)`\
[`points`](https://rspatial.github.io/terra/reference/lines.html)`(``hespero``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"x"``, ``"y"``)``]``, col ``=`` ``"black"``, pch ``=`` ``16``)`\
`cols`` ``<-`` `[`rep`](https://rdrr.io/r/base/rep.html)`(``"gray80"``, ``8``)`\
`cols``[``regions``$``category`` ``==`` ``"SCR"``]`` ``<-`` ``"yellow"`\
`terra``::`[`inset`](https://rspatial.github.io/terra/reference/inset.html)`(`\
`  ``regions``,`\
`  loc ``=`` ``"bottomleft"``,`\
`  scale ``=`` ``.3``,`\
`  col ``=`` ``cols`\
`)`

![](v06_Extrapolation_example_files/figure-html/raw%20data-1.png)

## Delimit calibration area

First, we must define our model’s calibration area. The *flexsdm*
package offers several methods for defining the model calibration area.
Here, we will use 25-km buffer areas around the presence points to
select our pseudo-absence locations.

\
`ca`` ``<-`` `[`calib_area`](https://sjevelazco.github.io/flexsdm/reference/calib_area.md)`(`\
`  data ``=`` ``hespero``,`\
`  x ``=`` ``"x"``,`\
`  y ``=`` ``"y"``,`\
`  method ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"buffer"``, width ``=`` ``25000``)``,`\
`  crs ``=`` `[`crs`](https://rspatial.github.io/terra/reference/crs.html)`(``somevar``)`\
`)`\
\
`# visualize the species occurrences & calibration area`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(`\
`  ``sp_region``,`\
`  col ``=`` ``"gray80"``,`\
`  legend ``=`` ``FALSE``,`\
`  axes ``=`` ``FALSE``,`\
`  main ``=`` ``"Calibration area and occurrences"`\
`)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``ca``, add ``=`` ``TRUE``)`\
[`points`](https://rspatial.github.io/terra/reference/lines.html)`(``hespero``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"x"``, ``"y"``)``]``, col ``=`` ``"black"``, pch ``=`` ``16``)`

![](v06_Extrapolation_example_files/figure-html/calibration%20area-1.png)

## Create pseudo-absence data

As is often the case with rare species, we only have species presence
data. However, most SDM methods require either pseudo-absence or
background point data. Here, we use our calibration area to produce
pseudo-absence data that can be used in our SDMs.

\
`# Sample the same number of species presences`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``10``)`\
`psa`` ``<-`` `[`sample_pseudoabs`](https://sjevelazco.github.io/flexsdm/reference/sample_pseudoabs.md)`(`\
`  data ``=`` ``hespero``,`\
`  x ``=`` ``"x"``,`\
`  y ``=`` ``"y"``,`\
`  n ``=`` `[`sum`](https://rdrr.io/r/base/sum.html)`(``hespero``$``pr_ab``)``, ``# number of pseudo-absence points equal to number of presences`\
`  method ``=`` ``"random"``,`\
`  rlayer ``=`` ``somevar``,`\
`  calibarea ``=`` ``ca`\
`)`\
\
`# Visualize species presences and pseudo-absences`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(`\
`  ``sp_region``,`\
`  col ``=`` ``"gray80"``,`\
`  legend ``=`` ``FALSE``,`\
`  axes ``=`` ``FALSE``,`\
`  xlim ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``289347``, ``353284``)``,`\
`  ylim ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``-``598052``, ``-``520709``)``,`\
`  main ``=`` ``"Presence = yellow, Pseudo-absence = black"`\
`)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``ca``, add ``=`` ``TRUE``)`\
[`points`](https://rspatial.github.io/terra/reference/lines.html)`(``psa``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"x"``, ``"y"``)``]``, cex ``=`` ``0.8``, pch ``=`` ``16``, col ``=`` ``"black"``)`` ``# Pseudo-absences`\
[`points`](https://rspatial.github.io/terra/reference/lines.html)`(``hespero``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"x"``, ``"y"``)``]``, col ``=`` ``"yellow"``, pch ``=`` ``16``, cex ``=`` ``1.5``)`` ``# Presences`

![](v06_Extrapolation_example_files/figure-html/pseudo-absence%20data-1.png)

\
\
\
`# Bind a presences and pseudo-absences`\
`hespero_pa`` ``<-`` `[`bind_rows`](https://dplyr.tidyverse.org/reference/bind_rows.html)`(``hespero``, ``psa``)`\
`hespero_pa`` ``# Presence-Pseudo-absence database`\
`#> ``# A tibble: 42 × 3`\
`#>          x        y pr_ab`\
`#>      ``<dbl>``    ``<dbl>`` ``<dbl>`\
`#> `` 1`` ``316``923. -``557``843.``     1`\
`#> `` 2`` ``317``155. -``559``234.``     1`\
`#> `` 3`` ``316``960. -``558``186.``     1`\
`#> `` 4`` ``314``347. -``559``648.``     1`\
`#> `` 5`` ``317``348. -``557``349.``     1`\
`#> `` 6`` ``316``753. -``559``679.``     1`\
`#> `` 7`` ``316``777. -``558``644.``     1`\
`#> `` 8`` ``317``050. -``559``043.``     1`\
`#> `` 9`` ``316``655. -``559``928.``     1`\
`#> ``10`` ``316``418. -``567``439.``     1`\
`#> ``# ℹ 32 more rows`

## Partition data for evaluating models

To evaluate model performance, we need to specify data for testing and
training. *flexsdm* offers a range of random and spatial and random data
partition methods for evaluating SDMs. Here we will use repeated K-fold
cross-validation, which is a suitable partition approach for validating
SDM with few data.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``10``)`\
\
`# Repeated K-fold method`\
`hespero_pa2`` ``<-`` `[`part_random`](https://sjevelazco.github.io/flexsdm/reference/part_random.md)`(`\
`  data ``=`` ``hespero_pa``,`\
`  pr_ab ``=`` ``"pr_ab"``,`\
`  method ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``method ``=`` ``"rep_kfold"``, folds ``=`` ``5``, replicates ``=`` ``3``)`\
`)`

## Extracting environmental values

Next, we extract the values of our four environmental predictors at the
presence and pseudo-absence locations.

\
`hespero_pa3`` ``<-`\
`  `[`sdm_extract`](https://sjevelazco.github.io/flexsdm/reference/sdm_extract.md)`(`\
`    data ``=`` ``hespero_pa2``,`\
`    x ``=`` ``"x"``,`\
`    y ``=`` ``"y"``,`\
`    env_layer ``=`` ``somevar``,`\
`    variables ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"cwd"``, ``"tmn"``, ``"aet"``, ``"ppt_jja"``)`\
`  ``)`

## Modeling

Let’s use three standard algorithms to model the distribution of
*Hesperocyparis stephensonii*: GLM, GBM, and SVM. In this case, we will
use the extent of the CFP as our prediction area so that we can evaluate
model extrapolation across a broad geographic area.

\
`mglm`` ``<-`\
`  `[`fit_glm`](https://sjevelazco.github.io/flexsdm/reference/fit_glm.md)`(`\
`    data ``=`` ``hespero_pa3``,`\
`    response ``=`` ``"pr_ab"``,`\
`    predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"cwd"``, ``"tmn"``, ``"aet"``, ``"ppt_jja"``)``,`\
`    partition ``=`` ``".part"``,`\
`    thr ``=`` ``"max_sens_spec"`\
`  ``)`\
\
`mgbm`` ``<-`` `[`fit_gbm`](https://sjevelazco.github.io/flexsdm/reference/fit_gbm.md)`(`\
`  data ``=`` ``hespero_pa3``,`\
`  response ``=`` ``"pr_ab"``,`\
`  predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"cwd"``, ``"tmn"``, ``"aet"``, ``"ppt_jja"``)``,`\
`  partition ``=`` ``".part"``,`\
`  thr ``=`` ``"max_sens_spec"`\
`)`\
\
`msvm`` ``<-`` `[`fit_svm`](https://sjevelazco.github.io/flexsdm/reference/fit_svm.md)`(`\
`  data ``=`` ``hespero_pa3``,`\
`  response ``=`` ``"pr_ab"``,`\
`  predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"cwd"``, ``"tmn"``, ``"aet"``, ``"ppt_jja"``)``,`\
`  partition ``=`` ``".part"``,`\
`  thr ``=`` ``"max_sens_spec"`\
`)`\
\
\
`mpred`` ``<-`` `[`sdm_predict`](https://sjevelazco.github.io/flexsdm/reference/sdm_predict.md)`(`\
`  models ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``mglm``, ``mgbm``, ``msvm``)``,`\
`  pred ``=`` ``somevar``,`\
`  con_thr ``=`` ``TRUE``,`\
`  predict_area ``=`` ``NULL`\
`)`

## Comparing our models

First, let’s take a look at the spatial predictions for our models. GLM
and GBM predict a lot of suitable habitat very far from where the
species is found!

\
[`par`](https://rdrr.io/r/graphics/par.html)`(``mfrow ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``3``)``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``mpred``$``glm``, main ``=`` ``"GLM"``)`\
`# points(hespero$x, hespero$y, pch = 19)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``mpred``$``gbm``, main ``=`` ``"GBM"``)`\
`# points(hespero$x, hespero$y, pch = 19)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``mpred``$``svm``, main ``=`` ``"SVM"``)`

![](v06_Extrapolation_example_files/figure-html/comparison%20maps-1.png)

\
`# points(hespero$x, hespero$y, pch = 19)`

## Partial dependence plots to explore the impact of predictor conditions on suitability

Extrapolation reflects an issue with how a model handles novel data.
Here, we see that the three algorithms explored in this tutorial predict
pretty different geographic patterns of habitat suitability based on the
same occurrence/pseudo-absence data and environmental predictors. Let’s
take a look at some partial dependence plots to see that the marginal
effect of each of the environmental predictors on suitability looks like
for each of our test models. This function allows you to visualize the
how a model may extrapolate outside the environmental conditions used in
training, by visualizing the “projection” data in a different color. In
this case, that will be our environmental predictors that cover the
extent of the CFP. *flexsdm* allows users to plot univariate partial
dependence plots
([`p_pdp`](https://sjevelazco.github.io/flexsdm/reference/p_pdp.html))
and bivariate partial dependence plots
([`p_bpdp`](https://sjevelazco.github.io/flexsdm/reference/p_bpdp.html));
both are shown below for each model. Note: the p_bpdp function allows
users the option to show the boundaries for the training data using
either a rectangle or convex hull approach. Here we will use the convex
hull approach.

Uni and bivariate partial dependence plots for the GLM:

\
[`p_pdp`](https://sjevelazco.github.io/flexsdm/reference/p_pdp.md)`(`\
`  model ``=`` ``mglm``$``model``,`\
`  training_data ``=`` ``hespero_pa3``,`\
`  projection_data ``=`` ``somevar`\
`)`

![](v06_Extrapolation_example_files/figure-html/glm%20partial%20dependence%20plots-1.png)

\
[`p_bpdp`](https://sjevelazco.github.io/flexsdm/reference/p_bpdp.md)`(`\
`  model ``=`` ``mglm``$``model``,`\
`  training_data ``=`` ``hespero_pa3``,`\
`  training_boundaries ``=`` ``"convexh"`\
`)`

![](v06_Extrapolation_example_files/figure-html/glm%20partial%20dependence%20plots-2.png)

Uni and bivariate partial dependence plots for the GBM:

\
[`p_pdp`](https://sjevelazco.github.io/flexsdm/reference/p_pdp.md)`(`\
`  model ``=`` ``mgbm``$``model``,`\
`  training_data ``=`` ``hespero_pa3``,`\
`  projection_data ``=`` ``somevar`\
`)`

![](v06_Extrapolation_example_files/figure-html/gbm%20partial%20dependence%20plots-1.png)

\
[`p_bpdp`](https://sjevelazco.github.io/flexsdm/reference/p_bpdp.md)`(`\
`  model ``=`` ``mgbm``$``model``,`\
`  training_data ``=`` ``hespero_pa3``,`\
`  training_boundaries ``=`` ``"convexh"``,`\
`  resolution ``=`` ``100`\
`)`

![](v06_Extrapolation_example_files/figure-html/gbm%20partial%20dependence%20plots-2.png)

Uni and bivariate partial dependence plots for the SVM:

\
[`p_pdp`](https://sjevelazco.github.io/flexsdm/reference/p_pdp.md)`(`\
`  model ``=`` ``msvm``$``model``,`\
`  training_data ``=`` ``hespero_pa3``,`\
`  projection_data ``=`` ``somevar`\
`)`

![](v06_Extrapolation_example_files/figure-html/svm%20partial%20dependence%20plots-1.png)

\
[`p_bpdp`](https://sjevelazco.github.io/flexsdm/reference/p_bpdp.md)`(`\
`  model ``=`` ``msvm``$``model``,`\
`  training_data ``=`` ``hespero_pa3``,`\
`  training_boundaries ``=`` ``"convexh"`\
`)`

![](v06_Extrapolation_example_files/figure-html/svm%20partial%20dependence%20plots-2.png)

These plots show a really interesting story! Most notably, the GLM and
GBM show consistently high habitat suitability for those areas that have
much higher actual evapotranspiration than the very narrow range of
values that were used to train the model. However, the SVM seems to do
the best job of not estimating very high habitat suitability for
environmental values that were outside of the training data.
Importantly, these models can behave very differently depending on the
modeling situation and context.

## Extrapolation evaluation

Remember that our species is highly restricted to southern California!
However, two of our models (GLM and GBM) predict very high habitat
suitability throughout other parts of the CFP, while the SVM provides
more conservative predictions. We see that GLM and GBM tend to predict
high habitat suitability in those areas that are very environmentally
different from our training conditions. But where are our models
extrapolating in environmental space? Let’s find out using the
“extra_eval” function in SDM. This function requires you to input the
model training data, a column specifying presence vs. absence locations,
projection data (can be a SpatRaster or a tibble containing data used
for model projection – this can reflect a larger region, separate
region, or different time period than what was used for model training),
a metric for calculating the degree of extrapolation (the default is
Mahalanobis distance, though euclidean is also an option- we will
explore both), number of cores for parallel processing, and an
aggregation factor, in case you want to measure extrapolation for a very
large data set.

First we look at the degree of extrapolation in geographic space using
the Shape method based on Mahalanobis distance. Also we will distinguish
between univariate and combinatorial extrapolation.

Using Mahalanobis distance:

\
`xp_m`` ``<-`\
`  `[`extra_eval`](https://sjevelazco.github.io/flexsdm/reference/extra_eval.md)`(`\
`    training_data ``=`` ``hespero_pa3``,`\
`    pr_ab ``=`` ``"pr_ab"``,`\
`    projection_data ``=`` ``somevar``,`\
`    metric ``=`` ``"mahalanobis"``,`\
`    univar_comb ``=`` ``TRUE``,`\
`    aggreg_factor ``=`` ``1`\
`  ``)`\
`xp_m`\
`#> class       : SpatRaster`\
`#> size        : 558, 394, 2  (nrow, ncol, nlyr)`\
`#> resolution  : 1890, 1890  (x, y)`\
`#> extent      : -373685.8, 370974.2, -604813.3, 449806.7  (xmin, xmax, ymin, ymax)`\
`#> coord. ref. : +proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs`\
`#> source(s)   : memory`\
`#> varnames    : somevar`\
`#>               `\
`#> names       : extrapolation, uni_comb`\
`#> min values  :             0,        1`\
`#> max values  :    3730.67743,        2`

The output of the extra_eval function is a SpatRaster, showing the
degree of extrapolation across the projection area, as estimated by the
Shape method.

\
`cl`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`  ``"#FDE725"``,`\
`  ``"#B3DC2B"``,`\
`  ``"#6DCC57"``,`\
`  ``"#36B677"``,`\
`  ``"#1F9D87"``,`\
`  ``"#25818E"``,`\
`  ``"#30678D"``,`\
`  ``"#3D4988"``,`\
`  ``"#462777"``,`\
`  ``"#440154"`\
`)`\
\
[`par`](https://rdrr.io/r/graphics/par.html)`(``mfrow ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``2``)``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``xp_m``$``extrapolation``, main ``=`` ``"Shape metric"``, col ``=`` ``cl``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(`\
`  ``xp_m``$``uni_comb``,`\
`  main ``=`` ``"Univariate (1) and \n combinatorial (2) extrapolation"``,`\
`  col ``=`` ``cl`\
`)`

![](v06_Extrapolation_example_files/figure-html/comparison%20extrapolation%20outputs-1.png)

We can also explore extrapolation or suitability patterns in
environmental and geographic space, using just one function. To do that,
we will use the
[`p_extra`](https://sjevelazco.github.io/flexsdm/reference/p_extra.html)
function. This function plots a ggplot object.

Let’s start with our extrapolation evaluation. These plots show that
areas with high extrapolation (dark blue) are far from the training data
(shown in black) in both environmental and geographic space.

The higher extrapolation values extrapolation area in the northwestern
portion of the CFP.

\
[`p_extra`](https://sjevelazco.github.io/flexsdm/reference/p_extra.md)`(`\
`  training_data ``=`` ``hespero_pa3``,`\
`  x ``=`` ``"x"``,`\
`  y ``=`` ``"y"``,`\
`  pr_ab ``=`` ``"pr_ab"``,`\
`  color_p ``=`` ``"black"``,`\
`  extra_suit_data ``=`` ``xp_m``,`\
`  projection_data ``=`` ``somevar``,`\
`  geo_space ``=`` ``TRUE``,`\
`  prop_points ``=`` ``0.05`\
`)`\
`#> Number of cell used to plot 3642 (5%)`

![](v06_Extrapolation_example_files/figure-html/graphical%20explore%20-%20Mahalanobis-1.png)

Let’s explore univariate and combinatorial extrapolation. The former is
defined as the projecting data outside range of training conditions,
while the combinatorial extrapolation area those projecting data within
the range of training conditions.

\
[`p_extra`](https://sjevelazco.github.io/flexsdm/reference/p_extra.md)`(`\
`  training_data ``=`` ``hespero_pa3``,`\
`  x ``=`` ``"x"``,`\
`  y ``=`` ``"y"``,`\
`  pr_ab ``=`` ``"pr_ab"``,`\
`  color_p ``=`` ``"black"``,`\
`  extra_suit_data ``=`` ``xp_m``$``uni_comb``,`\
`  projection_data ``=`` ``somevar``,`\
`  geo_space ``=`` ``TRUE``,`\
`  prop_points ``=`` ``0.05``,`\
`  color_gradient ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#B3DC2B"``, ``"#30678D"``)``,`\
`  alpha_p ``=`` ``0.2`\
`)`\
`#> Number of cell used to plot 3642 (5%)`

![](v06_Extrapolation_example_files/figure-html/graphical%20explore%20-%20uni_comb%20extrapolation-1.png)

## Truncating SDMs predictions based on extrapolation thresholds

Depending on the user’s end goal, you may want to exclude suitability
values that are environmentally “too” far from modeling training data.
The Shape method allows you to select any extrapolation threshold to
exclude suitability values.

Before truncating our models we can use the
[`p_extra`](https://sjevelazco.github.io/flexsdm/reference/p_extra.html)
function to explore binary extrapolation patter in the environmental and
geographical space. Here we will test the values 50, 100, and 500, for
comparison.

\
[`p_extra`](https://sjevelazco.github.io/flexsdm/reference/p_extra.md)`(`\
`  training_data ``=`` ``hespero_pa3``,`\
`  x ``=`` ``"x"``,`\
`  y ``=`` ``"y"``,`\
`  pr_ab ``=`` ``"pr_ab"``,`\
`  color_p ``=`` ``"black"``,`\
`  extra_suit_data ``=`` `[`as.numeric`](https://rdrr.io/r/base/numeric.html)`(``xp_m``$``extrapolation`` ``<`` ``50``)``,`\
`  projection_data ``=`` ``somevar``,`\
`  geo_space ``=`` ``TRUE``,`\
`  prop_points ``=`` ``0.05``,`\
`  color_gradient ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"gray"``, ``"#FDE725"``)``,`\
`  alpha_p ``=`` ``0.5`\
`)`` ``+`\
`  `[`plot_annotation`](https://patchwork.data-imaginist.com/reference/plot_annotation.html)`(`\
`    subtitle ``=`` ``"Binary extrapolation pattern with using a threshold of 50"`\
`  ``)`\
`#> Number of cell used to plot 3642 (5%)`

![](v06_Extrapolation_example_files/figure-html/explore%20extrapolation%20thresholds-1.png)

\
\
[`p_extra`](https://sjevelazco.github.io/flexsdm/reference/p_extra.md)`(`\
`  training_data ``=`` ``hespero_pa3``,`\
`  x ``=`` ``"x"``,`\
`  y ``=`` ``"y"``,`\
`  pr_ab ``=`` ``"pr_ab"``,`\
`  color_p ``=`` ``"black"``,`\
`  extra_suit_data ``=`` `[`as.numeric`](https://rdrr.io/r/base/numeric.html)`(``xp_m``$``extrapolation`` ``<`` ``100``)``,`\
`  projection_data ``=`` ``somevar``,`\
`  geo_space ``=`` ``TRUE``,`\
`  prop_points ``=`` ``0.05``,`\
`  color_gradient ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"gray"``, ``"#FDE725"``)``,`\
`  alpha_p ``=`` ``0.5`\
`)`` ``+`\
`  `[`plot_annotation`](https://patchwork.data-imaginist.com/reference/plot_annotation.html)`(`\
`    subtitle ``=`` ``"Binary extrapolation pattern with using a threshold of 100"`\
`  ``)`\
`#> Number of cell used to plot 3642 (5%)`

![](v06_Extrapolation_example_files/figure-html/explore%20extrapolation%20thresholds-2.png)

\
\
[`p_extra`](https://sjevelazco.github.io/flexsdm/reference/p_extra.md)`(`\
`  training_data ``=`` ``hespero_pa3``,`\
`  x ``=`` ``"x"``,`\
`  y ``=`` ``"y"``,`\
`  pr_ab ``=`` ``"pr_ab"``,`\
`  color_p ``=`` ``"black"``,`\
`  extra_suit_data ``=`` `[`as.numeric`](https://rdrr.io/r/base/numeric.html)`(``xp_m``$``extrapolation`` ``<`` ``500``)``,`\
`  projection_data ``=`` ``somevar``,`\
`  geo_space ``=`` ``TRUE``,`\
`  prop_points ``=`` ``0.05``,`\
`  color_gradient ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"gray"``, ``"#FDE725"``)``,`\
`  alpha_p ``=`` ``0.5`\
`)`` ``+`\
`  `[`plot_annotation`](https://patchwork.data-imaginist.com/reference/plot_annotation.html)`(`\
`    subtitle ``=`` ``"Binary extrapolation pattern with using a threshold of 500"`\
`  ``)`\
`#> Number of cell used to plot 3642 (5%)`

![](v06_Extrapolation_example_files/figure-html/explore%20extrapolation%20thresholds-3.png)

Values of 1 (yellow one) depict the environmental and geographical
regions will constraint our models suitability (truncate). Note that the
lower the threshold, the more restrictive the environmental and
geographic regions used to constrain the model.

Now we will use the function
[`extra_truncate`](https://sjevelazco.github.io/flexsdm/reference/extra_truncate.html)
to truncate the suitability predictions made by GLM, GBM, and SVM based
on extrapolation thresholds explored previously. As a note, threshold
selection will be very user-dependent, but this function allows you to
select multiple thresholds at one time to compare outputs. Users can
also select a “trunc_value” within the extra_truncate function, that
specifies the value that should be assigned to those cells that exceed
the extrapolation threshold (also specified in the function). The
default is 0 but users could also choose another value for which to
reduce suitability.

\
`glm_trunc`` ``<-`` `[`extra_truncate`](https://sjevelazco.github.io/flexsdm/reference/extra_truncate.md)`(`\
`  suit ``=`` ``mpred``$``glm``,`\
`  extra ``=`` ``xp_m``,`\
`  threshold ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``50``, ``100``, ``500``)``,`\
`  trunc_value ``=`` ``0`\
`)`\
\
`gbm_trunc`` ``<-`` `[`extra_truncate`](https://sjevelazco.github.io/flexsdm/reference/extra_truncate.md)`(`\
`  suit ``=`` ``mpred``$``gbm``,`\
`  extra ``=`` ``xp_m``,`\
`  threshold ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``50``, ``100``, ``500``)``,`\
`  trunc_value ``=`` ``0`\
`)`\
\
`svm_trunc`` ``<-`` `[`extra_truncate`](https://sjevelazco.github.io/flexsdm/reference/extra_truncate.md)`(`\
`  suit ``=`` ``mpred``$``svm``,`\
`  extra ``=`` ``xp_m``,`\
`  threshold ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``50``, ``100``, ``500``)``,`\
`  trunc_value ``=`` ``0`\
`)`

\
[`par`](https://rdrr.io/r/graphics/par.html)`(``mfrow ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``3``, ``3``)``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``glm_trunc``$``` `50` ```, main ``=`` ``"GLM; extra threshold = 50"``, col ``=`` ``cl``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``glm_trunc``$``` `100` ```, main ``=`` ``"GLM; extra threshold = 100"``, col ``=`` ``cl``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``glm_trunc``$``` `500` ```, main ``=`` ``"GLM; extra threshold = 500"``, col ``=`` ``cl``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``gbm_trunc``$``` `50` ```, main ``=`` ``"GBM; extra threshold = 50"``, col ``=`` ``cl``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``gbm_trunc``$``` `100` ```, main ``=`` ``"GBM; extra threshold = 100"``, col ``=`` ``cl``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``gbm_trunc``$``` `500` ```, main ``=`` ``"GBM; extra threshold = 500"``, col ``=`` ``cl``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``svm_trunc``$``` `50` ```, main ``=`` ``"SVM; extra threshold = 50"``, col ``=`` ``cl``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``svm_trunc``$``` `100` ```, main ``=`` ``"SVM; extra threshold = 100"``, col ``=`` ``cl``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``svm_trunc``$``` `500` ```, main ``=`` ``"SVM; extra threshold = 500"``, col ``=`` ``cl``)`

![](v06_Extrapolation_example_files/figure-html/comparison%20truncated%20outputs-1.png)

Based on these maps, you can see that the lower the extrapolation
threshold, the more restricted the habitat suitability patterns, while
higher values retain a greater amount of suitable habitat. Selecting the
best threshold will depend on modeling goals and objectives, and .

Want to learn more about Shape and other extrapolation metrics? Read the
article “Velazco, S. J. E., Brooke, M. R., De Marco Jr., P., Regan, H.
M., & Franklin, J. (2023). How far can I extrapolate my species
distribution model? Exploring Shape, a novel method. *Ecography*, *11*,
e06992. <https://doi.org/10.1111/ecog.06992>”
