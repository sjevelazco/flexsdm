# flexsdm: Overview of Modeling functions

## Introduction

Species distribution modeling (SDM) has become a standard tool in
multiple research areas, including ecology, conservation biology,
biogeography, paleobiogeography, and epidemiology. SDM is an area of
active theoretical and methodological research. The *flexsdm* package
provides users the ability to manipulate and parameterize models in a
variety of ways that meet their unique research needs.

This flexibility enables users to define their own complete or partial
modeling procedure specific for their modeling situations (e.g., number
of variables, number of records, different algorithms and ensemble
methods, algorithms tuning, etc.).

In this vignette, users will learn about the second set of functions in
the *flexsdm* package that fall under the “modeling” umbrella. These
functions were designed to construct and validate different types of
models and can be grouped into fit\_\* , tune\_\* , and esm\_\* family
functions. In addition there is a function to perform ensemble modeling.

The fit\_\* functions construct and validate models with default
hyper-parameter values. The tune\_\* functions construct and validate
models by searching for the best combination of hyper-parameter values,
and esm\_ functions can be used for constructing and validating Ensemble
of Small Models. Finally, the fit_ensemble() function is for fitting and
validating ensemble models.

These are the functions for model construction and validation:

**fit\_\* functions family**

- fit_gam() Fit and validate Generalized Additive Models

- fit_gau() Fit and validate Gaussian Process models

- fit_gbm() Fit and validate Generalized Boosted Regression models

- fit_glm() Fit and validate Generalized Linear Models

- fit_max() Fit and validate Maximum Entropy models

- fit_net() Fit and validate Neural Networks models

- fit_raf() Fit and validate Random Forest models

- fit_svm() Fit and validate Support Vector Machine models

**tune\_\* functions family**

- tune_gbm() Fit and validate Generalized Boosted Regression models with
  exploration of hyper-parameters

- tune_max() Fit and validate Maximum Entropy models with exploration of
  hyper-parameters

- tune_net() Fit and validate Neural Networks models with exploration of
  hyper-parameters

- tune_raf() Fit and validate Random Forest models with exploration of
  hyper-parameters

- tune_svm() Fit and validate Support Vector Machine models with
  exploration of hyper-parameters

**model ensemble**

- fit_ensemble() Fit and validate ensemble models with different
  ensemble methods

**esm\_\* functions family**

- esm_gam() Fit and validate Generalized Additive Models with Ensemble
  of Small Model approach

- esm_gau() Fit and validate Gaussian Process models Models with
  Ensemble of Small Model approach

- esm_gbm() Fit and validate Generalized Boosted Regression models with
  Ensemble of Small Model approach

- esm_glm() Fit and validate Generalized Linear Models with Ensemble of
  Small Model approach

- esm_max() Fit and validate Maximum Entropy models with Ensemble of
  Small Model approach

- esm_net() Fit and validate Neural Networks models with Ensemble of
  Small Model approach

- esm_svm() Fit and validate Support Vector Machine models with Ensemble
  of Small Model approach

## Installation

First, install the flexsdm package. You can install the released version
of *flexsdm* from [github](https://github.com/sjevelazco/flexsdm) with:

\
`# devtools::install_github('sjevelazco/flexsdm')`\
[`require`](https://rdrr.io/r/base/library.html)`(`[`flexsdm`](https://sjevelazco.github.io/flexsdm/)`)`\
`#> Loading required package: flexsdm`\
[`require`](https://rdrr.io/r/base/library.html)`(`[`terra`](https://rspatial.org/)`)`\
`#> Loading required package: terra`\
`#> terra 1.9.46`\
`#> `\
`#> Attaching package: 'terra'`\
`#> The following object is masked from 'package:knitr':`\
`#> `\
`#>     spin`\
[`require`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`\
`#> Loading required package: dplyr`\
`#> `\
`#> Attaching package: 'dplyr'`\
`#> The following objects are masked from 'package:terra':`\
`#> `\
`#>     intersect, union`\
`#> The following objects are masked from 'package:stats':`\
`#> `\
`#>     filter, lag`\
`#> The following objects are masked from 'package:base':`\
`#> `\
`#>     intersect, setdiff, setequal, union`

## Project directory setup

Decide where on your computer you would like to store the inputs and
outputs of your project (this will be your main directory). Use an
existing one or use dir.create() to create your main directory. Then
specify whether or not to include folders for projections, calibration
areas, algorithms, ensembles, and thresholds. For more details see
[Vignette
01_pre_modeling](https://sjevelazco.github.io/flexsdm/articles/v01_pre_modeling.html)

## Data, species occurrence and background data

In this tutorial, we will be using species occurrences and environmental
data that are available through the *flexsdm* package. The “abies”
example dataset includes a pr_ab column (presence = 1, and absence = 0),
location columns (x, y) and other environmental data. You can load the
“abies” data into your local R environment by using the code below:

(THIS EXAMPLE LOOKS A LITTLE STRANGE BECAUSE WE ARE ALSO USING
BACKGROUND DATA, WHILE THE ABIES DATASET CLEARLY HAS ABSENCES…)

\
[`data`](https://rdrr.io/r/utils/data.html)`(``"abies"``)`\
[`data`](https://rdrr.io/r/utils/data.html)`(``"backg"``)`\
`# Use a subset of the data to keep this tutorial fast to run;`\
`# in a real analysis you may want to use the full dataset`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``10``)`\
`abies`` ``<-`` ``abies`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  ``dplyr``::`[`group_by`](https://dplyr.tidyverse.org/reference/group_by.html)`(``pr_ab``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  ``dplyr``::`[`slice_sample`](https://dplyr.tidyverse.org/reference/slice.html)`(``prop ``=`` ``0.3``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  ``dplyr``::`[`ungroup`](https://dplyr.tidyverse.org/reference/group_by.html)`(``)`\
`backg`` ``<-`` ``backg`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` ``dplyr``::`[`slice_sample`](https://dplyr.tidyverse.org/reference/slice.html)`(``n ``=`` ``1000``)`\
\
`dplyr``::`[`glimpse`](https://pillar.r-lib.org/reference/glimpse.html)`(``abies``)`\
`#> Rows: 420`\
`#> Columns: 13`\
`#> $ id       ``<int>`` 12040``, ``10361``, ``9402``, ``9815``, ``10524``, ``8860``, ``6431``, ``11730``, ``808``, ``1105…`\
`#> $ pr_ab    ``<dbl>`` 0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0…`\
`#> $ x        ``<dbl>`` -308908.77``, ``-254286.44``, ``-286978.67``, ``-291848.83``, ``-256658.45``, ``1…`\
`#> $ y        ``<dbl>`` 384247.811``, ``417157.885``, ``386206.009``, ``445594.587``, ``184437.725``, ``-…`\
`#> $ aet      ``<dbl>`` 572.9367``, ``259.6567``, ``587.2900``, ``443.1700``, ``355.3867``, ``354.0600``, ``4…`\
`#> $ cwd      ``<dbl>`` 332.0133``, ``469.4567``, ``375.9467``, ``454.9833``, ``567.6433``, ``733.3933``, ``5…`\
`#> $ tmin     ``<dbl>`` 4.8400``, ``2.9333``, ``6.4533``, ``4.3933``, ``5.8667``, ``3.9733``, ``4.8733``, ``6.726…`\
`#> $ ppt_djf  ``<dbl>`` 521.4311``, ``151.2758``, ``332.6133``, ``331.5974``, ``303.1179``, ``181.9209``, ``1…`\
`#> $ ppt_jja  ``<dbl>`` 48.7567``, ``15.0839``, ``15.6589``, ``19.0647``, ``10.5549``, ``9.8277``, ``7.6569``, ``…`\
`#> $ pH       ``<dbl>`` 5.631732``, ``6.202818``, ``5.500000``, ``6.000000``, ``5.200000``, ``0.000000``, ``5…`\
`#> $ awc      ``<dbl>`` 0.10841342``, ``0.09496477``, ``0.16000000``, ``0.07000000``, ``0.08000000``, ``0…`\
`#> $ depth    ``<dbl>`` 63.69832``, ``68.67968``, ``178.00000``, ``107.00000``, ``61.00000``, ``0.00000``, ``…`\
`#> $ landform ``<fct>`` 6``, ``4``, ``7``, ``11``, ``6``, ``10``, ``14``, ``12``, ``7``, ``7``, ``14``, ``10``, ``7``, ``11``, ``14``, ``10``, ``10``, ``…`\
`dplyr``::`[`glimpse`](https://pillar.r-lib.org/reference/glimpse.html)`(``backg``)`\
`#> Rows: 1,000`\
`#> Columns: 13`\
`#> $ pr_ab        ``<dbl>`` 0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``…`\
`#> $ x            ``<dbl>`` 322239.164``, ``18219.164``, ``138909.164``, ``-52520.836``, ``-265550.83…`\
`#> $ y            ``<dbl>`` -554188.334``, ``-319558.334``, ``-407038.334``, ``71941.666``, ``234481.…`\
`#> $ aet          ``<dbl>`` 594.6933``, ``196.3567``, ``349.2067``, ``504.1167``, ``511.7433``, ``450.810…`\
`#> $ cwd          ``<dbl>`` 564.1600``, ``1172.8033``, ``1030.2833``, ``698.0667``, ``330.4567``, ``479.9…`\
`#> $ tmin         ``<dbl>`` 4.7500``, ``7.1067``, ``11.1400``, ``6.9300``, ``3.5200``, ``3.9767``, ``9.7767``, ``…`\
`#> $ ppt_djf      ``<dbl>`` 129.3690``, ``36.1234``, ``104.8891``, ``187.1644``, ``312.0807``, ``86.2167``,``…`\
`#> $ ppt_jja      ``<dbl>`` 12.7164``, ``1.0201``, ``2.1754``, ``6.2898``, ``11.0941``, ``22.0211``, ``1.1523…`\
`#> $ pH           ``<dbl>`` 5.8000002``, ``7.8000002``, ``6.1700625``, ``5.6999998``, ``3.0788457``, ``5.…`\
`#> $ awc          ``<dbl>`` 0.13000000``, ``0.18000001``, ``0.17149687``, ``0.12000000``, ``0.0479597…`\
`#> $ depth        ``<dbl>`` 77.00000``, ``201.00000``, ``50.92476``, ``152.00000``, ``119.32758``, ``182.…`\
`#> $ percent_clay ``<dbl>`` 17.100000``, ``25.400000``, ``22.295300``, ``21.000000``, ``13.046448``, ``41…`\
`#> $ landform     ``<fct>`` 6``, ``13``, ``7``, ``14``, ``7``, ``7``, ``10``, ``13``, ``6``, ``11``, ``10``, ``10``, ``6``, ``10``, ``13``, ``11``,``…`

If you want to replace the abies dataset with your own data, make sure
that your dataset contains the environmental conditions related to
presence-absence data. We use a pre-modeling family function for a
k-fold partition method (to be used for cross-validation). In the
partition method the number of folds or replications must be the same
for presence-absence and for background points datasets.

\
`abies2`` ``<-`` `[`part_random`](https://sjevelazco.github.io/flexsdm/reference/part_random.md)`(`\
`  data ``=`` ``abies``,`\
`  pr_ab ``=`` ``"pr_ab"``,`\
`  method ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``method ``=`` ``"kfold"``, folds ``=`` ``3``)`\
`)`\
\
`dplyr``::`[`glimpse`](https://pillar.r-lib.org/reference/glimpse.html)`(``abies2``)`\
`#> Rows: 420`\
`#> Columns: 14`\
`#> $ id       ``<int>`` 12040``, ``10361``, ``9402``, ``9815``, ``10524``, ``8860``, ``6431``, ``11730``, ``808``, ``1105…`\
`#> $ pr_ab    ``<dbl>`` 0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0…`\
`#> $ x        ``<dbl>`` -308908.77``, ``-254286.44``, ``-286978.67``, ``-291848.83``, ``-256658.45``, ``1…`\
`#> $ y        ``<dbl>`` 384247.811``, ``417157.885``, ``386206.009``, ``445594.587``, ``184437.725``, ``-…`\
`#> $ aet      ``<dbl>`` 572.9367``, ``259.6567``, ``587.2900``, ``443.1700``, ``355.3867``, ``354.0600``, ``4…`\
`#> $ cwd      ``<dbl>`` 332.0133``, ``469.4567``, ``375.9467``, ``454.9833``, ``567.6433``, ``733.3933``, ``5…`\
`#> $ tmin     ``<dbl>`` 4.8400``, ``2.9333``, ``6.4533``, ``4.3933``, ``5.8667``, ``3.9733``, ``4.8733``, ``6.726…`\
`#> $ ppt_djf  ``<dbl>`` 521.4311``, ``151.2758``, ``332.6133``, ``331.5974``, ``303.1179``, ``181.9209``, ``1…`\
`#> $ ppt_jja  ``<dbl>`` 48.7567``, ``15.0839``, ``15.6589``, ``19.0647``, ``10.5549``, ``9.8277``, ``7.6569``, ``…`\
`#> $ pH       ``<dbl>`` 5.631732``, ``6.202818``, ``5.500000``, ``6.000000``, ``5.200000``, ``0.000000``, ``5…`\
`#> $ awc      ``<dbl>`` 0.10841342``, ``0.09496477``, ``0.16000000``, ``0.07000000``, ``0.08000000``, ``0…`\
`#> $ depth    ``<dbl>`` 63.69832``, ``68.67968``, ``178.00000``, ``107.00000``, ``61.00000``, ``0.00000``, ``…`\
`#> $ landform ``<fct>`` 6``, ``4``, ``7``, ``11``, ``6``, ``10``, ``14``, ``12``, ``7``, ``7``, ``14``, ``10``, ``7``, ``11``, ``14``, ``10``, ``10``, ``…`\
`#> $ .part    ``<int>`` 3``, ``1``, ``3``, ``2``, ``2``, ``2``, ``1``, ``1``, ``2``, ``2``, ``1``, ``3``, ``3``, ``3``, ``3``, ``2``, ``2``, ``3``, ``2``, ``2``, ``2…`

Now, in the abies2 object we have a new column called “.part” with the 3
k-folds (1, 2, 3), indicating which partition each record (row) is a
member of. Next, we have to apply the same partition method and number
of folds to the environmental conditions of the background points.

\
`backg2`` ``<-`` `[`part_random`](https://sjevelazco.github.io/flexsdm/reference/part_random.md)`(`\
`  data ``=`` ``backg``,`\
`  pr_ab ``=`` ``"pr_ab"``,`\
`  method ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``method ``=`` ``"kfold"``, folds ``=`` ``3``)`\
`)`\
\
`dplyr``::`[`glimpse`](https://pillar.r-lib.org/reference/glimpse.html)`(``backg2``)`\
`#> Rows: 1,000`\
`#> Columns: 14`\
`#> $ pr_ab        ``<dbl>`` 0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``0``, ``…`\
`#> $ x            ``<dbl>`` 322239.164``, ``18219.164``, ``138909.164``, ``-52520.836``, ``-265550.83…`\
`#> $ y            ``<dbl>`` -554188.334``, ``-319558.334``, ``-407038.334``, ``71941.666``, ``234481.…`\
`#> $ aet          ``<dbl>`` 594.6933``, ``196.3567``, ``349.2067``, ``504.1167``, ``511.7433``, ``450.810…`\
`#> $ cwd          ``<dbl>`` 564.1600``, ``1172.8033``, ``1030.2833``, ``698.0667``, ``330.4567``, ``479.9…`\
`#> $ tmin         ``<dbl>`` 4.7500``, ``7.1067``, ``11.1400``, ``6.9300``, ``3.5200``, ``3.9767``, ``9.7767``, ``…`\
`#> $ ppt_djf      ``<dbl>`` 129.3690``, ``36.1234``, ``104.8891``, ``187.1644``, ``312.0807``, ``86.2167``,``…`\
`#> $ ppt_jja      ``<dbl>`` 12.7164``, ``1.0201``, ``2.1754``, ``6.2898``, ``11.0941``, ``22.0211``, ``1.1523…`\
`#> $ pH           ``<dbl>`` 5.8000002``, ``7.8000002``, ``6.1700625``, ``5.6999998``, ``3.0788457``, ``5.…`\
`#> $ awc          ``<dbl>`` 0.13000000``, ``0.18000001``, ``0.17149687``, ``0.12000000``, ``0.0479597…`\
`#> $ depth        ``<dbl>`` 77.00000``, ``201.00000``, ``50.92476``, ``152.00000``, ``119.32758``, ``182.…`\
`#> $ percent_clay ``<dbl>`` 17.100000``, ``25.400000``, ``22.295300``, ``21.000000``, ``13.046448``, ``41…`\
`#> $ landform     ``<fct>`` 6``, ``13``, ``7``, ``14``, ``7``, ``7``, ``10``, ``13``, ``6``, ``11``, ``10``, ``10``, ``6``, ``10``, ``13``, ``11``,``…`\
`#> $ .part        ``<int>`` 2``, ``2``, ``3``, ``1``, ``2``, ``3``, ``2``, ``2``, ``1``, ``3``, ``1``, ``1``, ``3``, ``1``, ``1``, ``3``, ``2``, ``2``, ``2``, ``…`

In backg2 object we have a new column called “.part” with the 3 k-folds
(1, 2, 3).

### 1. Fit and validate models

We fit and validate models: I. a maximum entropy model with default
hyper-parameter values (flexsdm::fit_max) and II. a random forest model
with exploration of hyper-parameters (flexsdm::tune_raf).

I. Maximum Entropy models with default hyper-parameter values.

\
`max_t1`` ``<-`` `[`fit_max`](https://sjevelazco.github.io/flexsdm/reference/fit_max.md)`(`\
`  data ``=`` ``abies2``,`\
`  response ``=`` ``"pr_ab"``,`\
`  predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"aet"``, ``"ppt_jja"``, ``"pH"``, ``"awc"``, ``"depth"``)``,`\
`  predictors_f ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"landform"``)``,`\
`  partition ``=`` ``".part"``,`\
`  background ``=`` ``backg2``,`\
`  thr ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"max_sens_spec"``, ``"equal_sens_spec"``, ``"max_sorensen"``)``,`\
`  clamp ``=`` ``TRUE``,`\
`  classes ``=`` ``"default"``,`\
`  pred_type ``=`` ``"cloglog"``,`\
`  regmult ``=`` ``1`\
`)`\
`#> Formula used for model fitting:`\
`#> ~aet + ppt_jja + pH + awc + depth + I(aet^2) + I(ppt_jja^2) + I(pH^2) + I(awc^2) + I(depth^2) + hinge(aet) + hinge(ppt_jja) + hinge(pH) + hinge(awc) + hinge(depth) + ppt_jja:aet + pH:aet + awc:aet + depth:aet + pH:ppt_jja + awc:ppt_jja + depth:ppt_jja + awc:pH + depth:pH + depth:awc + categorical(landform) - 1`\
`#> Replica number: 1/1`\
`#> Partition number: 1/3`\
`#> Partition number: 2/3`\
`#> Partition number: 3/3`

This function returns a list object with the following elements:

\
[`names`](https://rspatial.github.io/terra/reference/names.html)`(``max_t1``)`\
`#> [1] "model"            "predictors"       "performance"      "performance_part"`\
`#> [5] "data_ens"`

model: A “MaxEnt” class object. This object can be used for predicting.

\
[`options`](https://rdrr.io/r/base/options.html)`(``max.print ``=`` ``20``)`\
`max_t1``$``model`\
`#> `\
`#> Call:  glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",      weights = weights, lambda = 10^(seq(4, 0, length.out = 200)) *          sum(reg)/length(reg) * sum(p)/sum(weights), standardize = F,      penalty.factor = reg) `\
`#> `\
`#>     Df  %Dev  Lambda`\
`#> 1    0  0.00 22.3600`\
`#> 2    0  0.00 21.3500`\
`#> 3    0  0.00 20.3800`\
`#> 4    0  0.00 19.4600`\
`#> 5    0  0.00 18.5800`\
`#> 6    0  0.00 17.7400`\
`#>  [ reached 'max' / getOption("max.print") -- omitted 194 rows ]`

predictors: A tibble with quantitative (c column names) and qualitative
(f column names) variables use for modeling.

\
`max_t1``$``predictors`\
`#> ``# A tibble: 1 × 6`\
`#>   c1    c2      c3    c4    c5    f       `\
`#>   ``<chr>`` ``<chr>``   ``<chr>`` ``<chr>`` ``<chr>`` ``<chr>``   `\
`#> ``1`` aet   ppt_jja pH    awc   depth landform`

performance: The performance metric (see sdm_eval). Those metrics that
are threshold dependent are calculated based on the threshold specified
in the argument. We can see all the selected threshold values.

\
`max_t1``$``performance`\
`#> ``# A tibble: 3 × 33`\
`#>   model threshold     thr_value n_presences n_absences TPR_mean  TPR_sd TNR_mean`\
`#>   ``<chr>`` ``<chr>``             ``<dbl>``       ``<int>``      ``<int>``    ``<dbl>``   ``<dbl>``    ``<dbl>`\
`#> ``1`` max   equal_sens_s…     0.581         210        210    0.652 0.008``25``    0.652`\
`#> ``2`` max   max_sens_spec     0.491         210        210    0.767 0.083``7``     0.586`\
`#> ``3`` max   max_sorensen      0.468         210        210    0.943 0.037``8``     0.319`\
`#> ``# ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,`\
`#> ``#   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,`\
`#> ``#   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,`\
`#> ``#   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,`\
`#> ``#   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,`\
`#> ``#   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,`\
`#> ``#   IMAE_mean <dbl>, IMAE_sd <dbl>`

Predicted suitability for each test partition (row) based on the best
model. This database is used in fit_ensemble.

\
`max_t1``$``data_ens`\
`#> ``# A tibble: 420 × 5`\
`#>    rnames replicates part  pr_ab   pred`\
`#>    ``<chr>``  ``<chr>``      ``<chr>`` ``<dbl>``  ``<dbl>`\
`#> `` 1`` 2      .part      1         0 0.498 `\
`#> `` 2`` 7      .part      1         0 0.297 `\
`#> `` 3`` 8      .part      1         0 0.413 `\
`#> `` 4`` 11     .part      1         0 0.552 `\
`#> `` 5`` 22     .part      1         0 0.269 `\
`#> `` 6`` 24     .part      1         0 0.149 `\
`#> `` 7`` 29     .part      1         0 0.033``9`\
`#> `` 8`` 35     .part      1         0 0.143 `\
`#> `` 9`` 36     .part      1         0 0.371 `\
`#> ``10`` 37     .part      1         0 0.506 `\
`#> ``# ℹ 410 more rows`

II- Random forest models with exploration of hyper-parameters.

First, we create a data.frame that provides hyper-parameters values to
be tested. It Is recommended to generate this data.frame.
Hyper-parameter needed for tuning is ‘mtry’. The maximum mtry must be
equal to total number of predictors.

\
`tune_grid`` ``<-`\
`  `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(`\
`    mtry ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``2``, ``4``, ``6``)``,`\
`    ntree ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``500``, ``800``)`\
`  ``)`

We use the same data object abies2, with the same k-fold partition
method:

\
`rf_t`` ``<-`\
`  `[`tune_raf`](https://sjevelazco.github.io/flexsdm/reference/tune_raf.md)`(`\
`    data ``=`` ``abies2``,`\
`    response ``=`` ``"pr_ab"``,`\
`    predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`      ``"aet"``,`\
`      ``"cwd"``,`\
`      ``"tmin"``,`\
`      ``"ppt_djf"``,`\
`      ``"ppt_jja"``,`\
`      ``"pH"``,`\
`      ``"awc"``,`\
`      ``"depth"`\
`    ``)``,`\
`    predictors_f ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"landform"``)``,`\
`    partition ``=`` ``".part"``,`\
`    grid ``=`` ``tune_grid``,`\
`    thr ``=`` ``"max_sens_spec"``,`\
`    metric ``=`` ``"TSS"``,`\
`  ``)`\
`#> Formula used for model fitting:`\
`#> pr_ab ~ aet + cwd + tmin + ppt_djf + ppt_jja + pH + awc + depth + landform`\
`#> Tuning model...`\
`#> Replica number: 1/1`\
`#> Formula used for model fitting:`\
`#> pr_ab ~ aet + cwd + tmin + ppt_djf + ppt_jja + pH + awc + depth + landform`\
`#> Replica number: 1/1`\
`#> Partition number: 1/3`\
`#> Partition number: 2/3`\
`#> Partition number: 3/3`

Let’s see what the output object contains. This function returns a list
object with the following elements:

\
[`names`](https://rspatial.github.io/terra/reference/names.html)`(``rf_t``)`\
`#> [1] "model"             "predictors"        "performance"      `\
`#> [4] "performance_part"  "hyper_performance" "data_ens"`

model: A “randomForest” class object. This object can be used to see the
formula details, a basic summary o fthe model, and for predicting.

\
`rf_t``$``model`\
`#> `\
`#> Call:`\
`#>  randomForest(formula = formula1, data = data, mtry = mtry, ntree = ntree,      importance = TRUE, ) `\
`#>                Type of random forest: classification`\
`#>                      Number of trees: 500`\
`#> No. of variables tried at each split: 2`\
`#> `\
`#>         OOB estimate of  error rate: 14.05%`\
`#> Confusion matrix:`\
`#>     0   1 class.error`\
`#> 0 176  34   0.1619048`\
`#> 1  25 185   0.1190476`

predictors: A tibble with quantitative (c column names) and qualitative
(f column names) variables use for modeling.

\
`rf_t``$``predictors`\
`#> ``# A tibble: 1 × 9`\
`#>   c1    c2    c3    c4      c5      c6    c7    c8    f       `\
`#>   ``<chr>`` ``<chr>`` ``<chr>`` ``<chr>``   ``<chr>``   ``<chr>`` ``<chr>`` ``<chr>`` ``<chr>``   `\
`#> ``1`` aet   cwd   tmin  ppt_djf ppt_jja pH    awc   depth landform`

performance: The performance metric (see sdm_eval). Those metrics that
are threshold dependent are calculated based on the threshold specified
in the argument. We can see all the selected threshold values.

\
`rf_t``$``performance`\
`#> ``# A tibble: 1 × 35`\
`#>    mtry ntree model threshold   thr_value n_presences n_absences TPR_mean TPR_sd`\
`#>   ``<dbl>`` ``<dbl>`` ``<chr>`` ``<chr>``           ``<dbl>``       ``<int>``      ``<int>``    ``<dbl>``  ``<dbl>`\
`#> ``1``     2   500 raf   max_sens_s…     0.592         210        210    0.876 0.059``5`\
`#> ``# ℹ 26 more variables: TNR_mean <dbl>, TNR_sd <dbl>, W_TPR_TNR_mean <dbl>,`\
`#> ``#   W_TPR_TNR_sd <dbl>, SORENSEN_mean <dbl>, SORENSEN_sd <dbl>,`\
`#> ``#   JACCARD_mean <dbl>, JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>,`\
`#> ``#   OR_mean <dbl>, OR_sd <dbl>, TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>,`\
`#> ``#   KAPPA_sd <dbl>, MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,`\
`#> ``#   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,`\
`#> ``#   IMAE_mean <dbl>, IMAE_sd <dbl>`

Predicted suitability for each test partition (row) based on the best
model. This database is used in fit_ensemble.

\
`rf_t``$``data_ens`\
`#> ``# A tibble: 420 × 5`\
`#>    rnames replicates part  pr_ab  pred`\
`#>    ``<chr>``  ``<chr>``      ``<chr>`` ``<fct>`` ``<dbl>`\
`#> `` 1`` 2      .part      1     0     0.734`\
`#> `` 2`` 7      .part      1     0     0.172`\
`#> `` 3`` 8      .part      1     0     0.1  `\
`#> `` 4`` 11     .part      1     0     0.642`\
`#> `` 5`` 22     .part      1     0     0.028`\
`#> `` 6`` 24     .part      1     0     0.516`\
`#> `` 7`` 29     .part      1     0     0.03 `\
`#> `` 8`` 35     .part      1     0     0.516`\
`#> `` 9`` 36     .part      1     0     0.1  `\
`#> ``10`` 37     .part      1     0     0.352`\
`#> ``# ℹ 410 more rows`

These model objects can be used in flexsdm::fit_ensemble().

### 2. Model Ensemble

In this example we fit and validate and ensemble model using the two
model objects that were just created.

\
`# Fit and validate ensemble model`\
`an_ensemble`` ``<-`` `[`fit_ensemble`](https://sjevelazco.github.io/flexsdm/reference/fit_ensemble.md)`(`\
`  models ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``max_t1``, ``rf_t``)``,`\
`  ens_method ``=`` ``"meansup"``,`\
`  thr ``=`` ``NULL``,`\
`  thr_model ``=`` ``"max_sens_spec"``,`\
`  metric ``=`` ``"TSS"`\
`)`\
`#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%`

\
`# Outputs`\
[`names`](https://rspatial.github.io/terra/reference/names.html)`(``an_ensemble``)`\
`#> [1] "models"           "thr_metric"       "predictors"       "performance"     `\
`#> [5] "performance_part"`\
\
`an_ensemble``$``thr_metric`\
`#> [1] "max_sens_spec" "TSS_mean"`\
`an_ensemble``$``predictors`\
`#> ``# A tibble: 2 × 9`\
`#>   c1    c2      c3    c4      c5      f        c6    c7    c8   `\
`#>   ``<chr>`` ``<chr>``   ``<chr>`` ``<chr>``   ``<chr>``   ``<chr>``    ``<chr>`` ``<chr>`` ``<chr>`\
`#> ``1`` aet   ppt_jja pH    awc     depth   landform ``NA``    ``NA``    ``NA``   `\
`#> ``2`` aet   cwd     tmin  ppt_djf ppt_jja landform pH    awc   depth`\
`an_ensemble``$``performance`\
`#> ``# A tibble: 7 × 33`\
`#>   model   threshold   thr_value n_presences n_absences TPR_mean  TPR_sd TNR_mean`\
`#>   ``<chr>``   ``<chr>``           ``<dbl>``       ``<int>``      ``<int>``    ``<dbl>``   ``<dbl>``    ``<dbl>`\
`#> ``1`` meansup equal_sens…     0.538         210        210    0.852 0.008``25``    0.852`\
`#> ``2`` meansup lpt             0.124         210        210    1     0          0.519`\
`#> ``3`` meansup max_fpb         0.548         210        210    0.890 0.064``4``     0.838`\
`#> ``4`` meansup max_jaccard     0.548         210        210    0.890 0.064``4``     0.838`\
`#> ``5`` meansup max_sens_s…     0.56          210        210    0.876 0.059``5``     0.852`\
`#> ``6`` meansup max_sorens…     0.548         210        210    0.890 0.064``4``     0.838`\
`#> ``7`` meansup sensitivity     0.466         210        210    0.905 0.008``25``    0.771`\
`#> ``# ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,`\
`#> ``#   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,`\
`#> ``#   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,`\
`#> ``#   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,`\
`#> ``#   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,`\
`#> ``#   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,`\
`#> ``#   IMAE_mean <dbl>, IMAE_sd <dbl>`

### 3. Fit and validate models with Ensemble of Small Model approach

This method consists of creating bivariate models with all pair-wise
combinations of predictors and perform an ensemble based on the average
of suitability weighted by Somers’ D metric (D = 2 x (AUC -0.5)). ESM is
recommended for modeling species with very few occurrences. This
function does not allow categorical variables because the use of these
types of variables could be problematic when applied to species with few
occurrences. For more detail see Breiner et al. (2015, 2018)

\
[`data`](https://rdrr.io/r/utils/data.html)`(``"abies"``)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`\
\
`# Create a smaller subset of occurrences`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``10``)`\
`abies2`` ``<-`` ``abies`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`na.omit`](https://rspatial.github.io/terra/reference/na.omit.html)`(``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`group_by`](https://dplyr.tidyverse.org/reference/group_by.html)`(``pr_ab``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  ``dplyr``::`[`slice_sample`](https://dplyr.tidyverse.org/reference/slice.html)`(``n ``=`` ``10``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`group_by`](https://dplyr.tidyverse.org/reference/group_by.html)`(``)`

We can use different methods in the flexsdm::part_random function
according to our data. See
[part_random](https://sjevelazco.github.io/flexsdm/reference/part_random.html)
for more details.

\
`# Using k-fold partition method for model cross validation`\
`abies2`` ``<-`` `[`part_random`](https://sjevelazco.github.io/flexsdm/reference/part_random.md)`(`\
`  data ``=`` ``abies2``,`\
`  pr_ab ``=`` ``"pr_ab"``,`\
`  method ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``method ``=`` ``"kfold"``, folds ``=`` ``3``)`\
`)`\
`abies2`\
`#> ``# A tibble: 20 × 14`\
`#>       id pr_ab        x        y   aet   cwd   tmin ppt_djf ppt_jja    pH    awc`\
`#>    ``<int>`` ``<dbl>``    ``<dbl>``    ``<dbl>`` ``<dbl>`` ``<dbl>``  ``<dbl>``   ``<dbl>``   ``<dbl>`` ``<dbl>``  ``<dbl>`\
`#> `` 1`` ``12``040     0 -``308``909.``  ``384``248.  573.  332.  4.84     521.   48.8   5.63 0.108 `\
`#> `` 2`` ``10``361     0 -``254``286.``  ``417``158.  260.  469.  2.93     151.   15.1   6.20 0.095``0`\
`#> `` 3``  ``9``402     0 -``286``979.``  ``386``206.  587.  376.  6.45     333.   15.7   5.5  0.160 `\
`#> `` 4``  ``9``815     0 -``291``849.``  ``445``595.  443.  455.  4.39     332.   19.1   6    0.070``0`\
`#> `` 5`` ``10``524     0 -``256``658.``  ``184``438.  355.  568.  5.87     303.   10.6   5.20 0.080``0`\
`#> `` 6``  ``8``860     0  ``121``343. -``164``170.``  354.  733.  3.97     182.    9.83  0    0     `\
`#> `` 7``  ``6``431     0  ``107``903. -``122``968.``  461.  578.  4.87     161.    7.66  5.90 0.090``0`\
`#> `` 8`` ``11``730     0 -``333``903.``  ``431``238.  561.  364.  6.73     387.   25.2   5.80 0.130 `\
`#> `` 9``   808     0 -``150``163.``  ``357``180.  339.  564.  2.64     220.   15.3   6.40 0.100 `\
`#> ``10`` ``11``054     0 -``293``663.``  ``340``981.  477.  396.  3.89     332.   26.4   4.60 0.063``4`\
`#> ``11``  ``2``960     1  -``49``273.``  ``181``752.  512.  275.  0.920    319.   17.3   5.92 0.090``0`\
`#> ``12``  ``3``065     1  ``126``907. -``198``892.``  322.  544.  0.700    203.   10.6   5.60 0.110 `\
`#> ``13``  ``5``527     1  ``116``751. -``181``089.``  261.  537.  0.363    178.    7.43  0    0     `\
`#> ``14``  ``4``035     1  -``31``777.``  ``115``940.  394.  440.  2.07     298.   11.2   6.01 0.076``9`\
`#> ``15``  ``4``081     1   -``5``158.``   ``90``159.  301.  502.  0.703    203.   14.6   6.11 0.063``3`\
`#> ``16``  ``3``087     1  ``102``151. -``143``976.``  299.  425. -``2.08``     205.   13.4   3.88 0.110 `\
`#> ``17``  ``3``495     1  -``19``586.``   ``89``803.  438.  419.  2.13     189.   15.2   6.19 0.095``9`\
`#> ``18``  ``4``441     1   ``49``405.  -``60``502.``  362.  582.  2.42     218.    7.84  5.64 0.078``6`\
`#> ``19``   301     1 -``132``516.``  ``270``845.  367.  196. -``2.56``     422.   26.3   6.70 0.030``0`\
`#> ``20``  ``3``162     1   ``59``905.  -``53``634.``  319.  626.  1.99     212.    4.50  4.51 0.039``6`\
`#> ``# ℹ 3 more variables: depth <dbl>, landform <fct>, .part <int>`

This function constructs Generalized Additive Models using the Ensembles
of Small Models (ESM) approach (Breiner et al., 2015, 2018).

\
`# We set the model without threshold specification and with the kfold created above`\
`esm_gam_t1`` ``<-`` `[`esm_gam`](https://sjevelazco.github.io/flexsdm/reference/esm_gam.md)`(`\
`  data ``=`` ``abies2``,`\
`  response ``=`` ``"pr_ab"``,`\
`  predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"aet"``,`\
`    ``"cwd"``,`\
`    ``"tmin"``,`\
`    ``"ppt_djf"``,`\
`    ``"ppt_jja"``,`\
`    ``"pH"``,`\
`    ``"awc"``,`\
`    ``"depth"`\
`  ``)``,`\
`  partition ``=`` ``".part"``,`\
`  thr ``=`` ``NULL``,`\
`  k ``=`` ``2`\
`)`\
`#>   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%`

This function returns a list object with the following elements:

\
[`names`](https://rspatial.github.io/terra/reference/names.html)`(``esm_gam_t1``)`\
`#> [1] "esm_model"        "predictors"       "performance"      "performance_part"`

esm_model: A list with “GAM” class object for each bivariate model. This
object can be used for predicting using the ESM approachwith sdm_predict
function.

\
[`options`](https://rdrr.io/r/base/options.html)`(``max.print ``=`` ``10``)`` ``# If you don't want to see printed all the output`\
`esm_gam_t1``$``esm_model`\
`` #> $`0.398148148148148` ``\
`#> `\
`#> Family: binomial `\
`#> Link function: logit `\
`#> `\
`#> Formula:`\
`#> pr_ab ~ s(aet, k = 2) + s(cwd, k = 2)`\
`#> `\
`#> Estimated degrees of freedom:`\
`#> 1.81 1.73  total = 4.54 `\
`#> `\
`#> UBRE score: 0.2827405     `\
`#> `\
`` #> $`1` ``\
`#> `\
`#> Family: binomial `\
`#> Link function: logit `\
`#> `\
`#> Formula:`\
`#> pr_ab ~ s(aet, k = 2) + s(tmin, k = 2)`\
`#> `\
`#> Estimated degrees of freedom:`\
`#> 1 1  total = 3 `\
`#> `\
`#> UBRE score: -0.7     `\
`#> `\
`` #> $`0.0925925925925926` ``\
`#> `\
`#> Family: binomial `\
`#> Link function: logit `\
`#> `\
`#> Formula:`\
`#> pr_ab ~ s(aet, k = 2) + s(pH, k = 2)`\
`#> `\
`#> Estimated degrees of freedom:`\
`#> 1 1  total = 3 `\
`#> `\
`#> UBRE score: 0.4818906     `\
`#> `\
`` #> $`0.564814814814815` ``\
`#> `\
`#> Family: binomial `\
`#> Link function: logit `\
`#> `\
`#> Formula:`\
`#> pr_ab ~ s(aet, k = 2) + s(depth, k = 2)`\
`#> `\
`#> Estimated degrees of freedom:`\
`#> 1.00 1.12  total = 3.12 `\
`#> `\
`#> UBRE score: 0.08148345     `\
`#> `\
`` #> $`1` ``\
`#> `\
`#> Family: binomial `\
`#> Link function: logit `\
`#> `\
`#> Formula:`\
`#> pr_ab ~ s(cwd, k = 2) + s(tmin, k = 2)`\
`#> `\
`#> Estimated degrees of freedom:`\
`#> 1 1  total = 3 `\
`#> `\
`#> UBRE score: -0.7     `\
`#> `\
`` #> $`0.314814814814815` ``\
`#> `\
`#> Family: binomial `\
`#> Link function: logit `\
`#> `\
`#> Formula:`\
`#> pr_ab ~ s(cwd, k = 2) + s(ppt_djf, k = 2)`\
`#> `\
`#> Estimated degrees of freedom:`\
`#> 1 1  total = 3 `\
`#> `\
`#> UBRE score: 0.4360771     `\
`#> `\
`` #> $`0.875` ``\
`#> `\
`#> Family: binomial `\
`#> Link function: logit `\
`#> `\
`#> Formula:`\
`#> pr_ab ~ s(cwd, k = 2) + s(ppt_jja, k = 2)`\
`#> `\
`#> Estimated degrees of freedom:`\
`#> 1 1  total = 3 `\
`#> `\
`#> UBRE score: 0.107549     `\
`#> `\
`` #> $`0.553240740740741` ``\
`#> `\
`#> Family: binomial `\
`#> Link function: logit `\
`#> `\
`#> Formula:`\
`#> pr_ab ~ s(cwd, k = 2) + s(depth, k = 2)`\
`#> `\
`#> Estimated degrees of freedom:`\
`#> 1.00 1.38  total = 3.38 `\
`#> `\
`#> UBRE score: 0.3416266     `\
`#> `\
`` #> $`1` ``\
`#> `\
`#> Family: binomial `\
`#> Link function: logit `\
`#> `\
`#> Formula:`\
`#> pr_ab ~ s(tmin, k = 2) + s(ppt_djf, k = 2)`\
`#> `\
`#> Estimated degrees of freedom:`\
`#> 1 1  total = 3 `\
`#> `\
`#> UBRE score: -0.7     `\
`#> `\
`` #> $`1` ``\
`#> `\
`#> Family: binomial `\
`#> Link function: logit `\
`#> `\
`#> Formula:`\
`#> pr_ab ~ s(tmin, k = 2) + s(ppt_jja, k = 2)`\
`#> `\
`#> Estimated degrees of freedom:`\
`#> 1 1  total = 3 `\
`#> `\
`#> UBRE score: -0.7     `\
`#> `\
`#>  [ reached 'max' / getOption("max.print") -- omitted 9 entries ]`

predictors: A tibble with variables use for modeling.

\
`esm_gam_t1``$``predictors`\
`#> ``# A tibble: 1 × 8`\
`#>   c1    c2    c3    c4      c5      c6    c7    c8   `\
`#>   ``<chr>`` ``<chr>`` ``<chr>`` ``<chr>``   ``<chr>``   ``<chr>`` ``<chr>`` ``<chr>`\
`#> ``1`` aet   cwd   tmin  ppt_djf ppt_jja pH    awc   depth`

performance: Performance metric (see sdm_eval). Those threshold
dependent metrics are calculated based on the threshold specified in the
argument.

\
`esm_gam_t1``$``performance`\
`#> ``# A tibble: 7 × 33`\
`#>   model   threshold    thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean`\
`#>   ``<chr>``   ``<chr>``            ``<dbl>``       ``<int>``      ``<int>``    ``<dbl>``  ``<dbl>``    ``<dbl>`\
`#> ``1`` esm_gam equal_sens_…     0.607          10         10        1      0        1`\
`#> ``2`` esm_gam lpt              0.607          10         10        1      0        1`\
`#> ``3`` esm_gam max_fpb          0.607          10         10        1      0        1`\
`#> ``4`` esm_gam max_jaccard      0.607          10         10        1      0        1`\
`#> ``5`` esm_gam max_sens_sp…     0.607          10         10        1      0        1`\
`#> ``6`` esm_gam max_sorensen     0.607          10         10        1      0        1`\
`#> ``7`` esm_gam sensitivity      0.614          10         10        1      0        1`\
`#> ``# ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,`\
`#> ``#   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,`\
`#> ``#   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,`\
`#> ``#   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,`\
`#> ``#   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,`\
`#> ``#   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,`\
`#> ``#   IMAE_mean <dbl>, IMAE_sd <dbl>`

Now, we test the rep_kfold partition method. In this method ‘folds’
refers to the number of partitions for data partitioning and ‘replicate’
refers to the number of replicates. Both assume values \>=1.

\
`# Remove the previous k-fold partition`\
`abies2`` ``<-`` ``abies2`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``-`[`starts_with`](https://tidyselect.r-lib.org/reference/starts_with.html)`(``"."``)``)`\
\
`# Test with rep_kfold partition using 3 folds and 2 replicates`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``10``)`\
`abies2`` ``<-`` `[`part_random`](https://sjevelazco.github.io/flexsdm/reference/part_random.md)`(`\
`  data ``=`` ``abies2``,`\
`  pr_ab ``=`` ``"pr_ab"``,`\
`  method ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``method ``=`` ``"rep_kfold"``, folds ``=`` ``3``, replicates ``=`` ``2``)`\
`)`\
`abies2`\
`#> ``# A tibble: 20 × 15`\
`#>       id pr_ab        x        y   aet   cwd   tmin ppt_djf ppt_jja    pH    awc`\
`#>    ``<int>`` ``<dbl>``    ``<dbl>``    ``<dbl>`` ``<dbl>`` ``<dbl>``  ``<dbl>``   ``<dbl>``   ``<dbl>`` ``<dbl>``  ``<dbl>`\
`#> `` 1`` ``12``040     0 -``308``909.``  ``384``248.  573.  332.  4.84     521.   48.8   5.63 0.108 `\
`#> `` 2`` ``10``361     0 -``254``286.``  ``417``158.  260.  469.  2.93     151.   15.1   6.20 0.095``0`\
`#> `` 3``  ``9``402     0 -``286``979.``  ``386``206.  587.  376.  6.45     333.   15.7   5.5  0.160 `\
`#> `` 4``  ``9``815     0 -``291``849.``  ``445``595.  443.  455.  4.39     332.   19.1   6    0.070``0`\
`#> `` 5`` ``10``524     0 -``256``658.``  ``184``438.  355.  568.  5.87     303.   10.6   5.20 0.080``0`\
`#> `` 6``  ``8``860     0  ``121``343. -``164``170.``  354.  733.  3.97     182.    9.83  0    0     `\
`#> `` 7``  ``6``431     0  ``107``903. -``122``968.``  461.  578.  4.87     161.    7.66  5.90 0.090``0`\
`#> `` 8`` ``11``730     0 -``333``903.``  ``431``238.  561.  364.  6.73     387.   25.2   5.80 0.130 `\
`#> `` 9``   808     0 -``150``163.``  ``357``180.  339.  564.  2.64     220.   15.3   6.40 0.100 `\
`#> ``10`` ``11``054     0 -``293``663.``  ``340``981.  477.  396.  3.89     332.   26.4   4.60 0.063``4`\
`#> ``11``  ``2``960     1  -``49``273.``  ``181``752.  512.  275.  0.920    319.   17.3   5.92 0.090``0`\
`#> ``12``  ``3``065     1  ``126``907. -``198``892.``  322.  544.  0.700    203.   10.6   5.60 0.110 `\
`#> ``13``  ``5``527     1  ``116``751. -``181``089.``  261.  537.  0.363    178.    7.43  0    0     `\
`#> ``14``  ``4``035     1  -``31``777.``  ``115``940.  394.  440.  2.07     298.   11.2   6.01 0.076``9`\
`#> ``15``  ``4``081     1   -``5``158.``   ``90``159.  301.  502.  0.703    203.   14.6   6.11 0.063``3`\
`#> ``16``  ``3``087     1  ``102``151. -``143``976.``  299.  425. -``2.08``     205.   13.4   3.88 0.110 `\
`#> ``17``  ``3``495     1  -``19``586.``   ``89``803.  438.  419.  2.13     189.   15.2   6.19 0.095``9`\
`#> ``18``  ``4``441     1   ``49``405.  -``60``502.``  362.  582.  2.42     218.    7.84  5.64 0.078``6`\
`#> ``19``   301     1 -``132``516.``  ``270``845.  367.  196. -``2.56``     422.   26.3   6.70 0.030``0`\
`#> ``20``  ``3``162     1   ``59``905.  -``53``634.``  319.  626.  1.99     212.    4.50  4.51 0.039``6`\
`#> ``# ℹ 4 more variables: depth <dbl>, landform <fct>, .part1 <int>, .part2 <int>`

We use the new rep_kfold partition in the gam model

\
`esm_gam_t2`` ``<-`` `[`esm_gam`](https://sjevelazco.github.io/flexsdm/reference/esm_gam.md)`(`\
`  data ``=`` ``abies2``,`\
`  response ``=`` ``"pr_ab"``,`\
`  predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"aet"``,`\
`    ``"cwd"``,`\
`    ``"tmin"``,`\
`    ``"ppt_djf"`\
`  ``)``,`\
`  partition ``=`` ``".part"``,`\
`  thr ``=`` ``NULL``,`\
`  k ``=`` ``2`\
`)`\
`#>   |                                                                              |                                                                      |   0%  |                                                                              |============                                                          |  17%  |                                                                              |=======================                                               |  33%  |                                                                              |===================================                                   |  50%  |                                                                              |===============================================                       |  67%  |                                                                              |==========================================================            |  83%  |                                                                              |======================================================================| 100%`

Test with random bootstrap partitioning. In method ‘replicate’ refers to
the number of replicates (assumes a value \>=1), ‘proportion’ refers to
the proportion of occurrences used for model fitting (assumes a value
\>0 and \<=1). With this method we can configure the proportion of
training and testing data according to the species occurrences. In this
example, proportion=‘0.7’ indicates that 70% of data will be used for
model training, while 30% will be used for model testing. For this
method, the function will return .partX columns with “train” or “test”
words as the entries.

\
`# Remove the previous k-fold partition`\
`abies2`` ``<-`` ``abies2`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``-`[`starts_with`](https://tidyselect.r-lib.org/reference/starts_with.html)`(``"."``)``)`\
\
`# Test with bootstrap partition using 3 replicates`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``10``)`\
`abies2`` ``<-`` `[`part_random`](https://sjevelazco.github.io/flexsdm/reference/part_random.md)`(`\
`  data ``=`` ``abies2``,`\
`  pr_ab ``=`` ``"pr_ab"``,`\
`  method ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``method ``=`` ``"boot"``, replicates ``=`` ``3``, proportion ``=`` ``0.7``)`\
`)`\
`abies2`\
`#> ``# A tibble: 20 × 16`\
`#>       id pr_ab        x        y   aet   cwd   tmin ppt_djf ppt_jja    pH    awc`\
`#>    ``<int>`` ``<dbl>``    ``<dbl>``    ``<dbl>`` ``<dbl>`` ``<dbl>``  ``<dbl>``   ``<dbl>``   ``<dbl>`` ``<dbl>``  ``<dbl>`\
`#> `` 1`` ``12``040     0 -``308``909.``  ``384``248.  573.  332.  4.84     521.   48.8   5.63 0.108 `\
`#> `` 2`` ``10``361     0 -``254``286.``  ``417``158.  260.  469.  2.93     151.   15.1   6.20 0.095``0`\
`#> `` 3``  ``9``402     0 -``286``979.``  ``386``206.  587.  376.  6.45     333.   15.7   5.5  0.160 `\
`#> `` 4``  ``9``815     0 -``291``849.``  ``445``595.  443.  455.  4.39     332.   19.1   6    0.070``0`\
`#> `` 5`` ``10``524     0 -``256``658.``  ``184``438.  355.  568.  5.87     303.   10.6   5.20 0.080``0`\
`#> `` 6``  ``8``860     0  ``121``343. -``164``170.``  354.  733.  3.97     182.    9.83  0    0     `\
`#> `` 7``  ``6``431     0  ``107``903. -``122``968.``  461.  578.  4.87     161.    7.66  5.90 0.090``0`\
`#> `` 8`` ``11``730     0 -``333``903.``  ``431``238.  561.  364.  6.73     387.   25.2   5.80 0.130 `\
`#> `` 9``   808     0 -``150``163.``  ``357``180.  339.  564.  2.64     220.   15.3   6.40 0.100 `\
`#> ``10`` ``11``054     0 -``293``663.``  ``340``981.  477.  396.  3.89     332.   26.4   4.60 0.063``4`\
`#> ``11``  ``2``960     1  -``49``273.``  ``181``752.  512.  275.  0.920    319.   17.3   5.92 0.090``0`\
`#> ``12``  ``3``065     1  ``126``907. -``198``892.``  322.  544.  0.700    203.   10.6   5.60 0.110 `\
`#> ``13``  ``5``527     1  ``116``751. -``181``089.``  261.  537.  0.363    178.    7.43  0    0     `\
`#> ``14``  ``4``035     1  -``31``777.``  ``115``940.  394.  440.  2.07     298.   11.2   6.01 0.076``9`\
`#> ``15``  ``4``081     1   -``5``158.``   ``90``159.  301.  502.  0.703    203.   14.6   6.11 0.063``3`\
`#> ``16``  ``3``087     1  ``102``151. -``143``976.``  299.  425. -``2.08``     205.   13.4   3.88 0.110 `\
`#> ``17``  ``3``495     1  -``19``586.``   ``89``803.  438.  419.  2.13     189.   15.2   6.19 0.095``9`\
`#> ``18``  ``4``441     1   ``49``405.  -``60``502.``  362.  582.  2.42     218.    7.84  5.64 0.078``6`\
`#> ``19``   301     1 -``132``516.``  ``270``845.  367.  196. -``2.56``     422.   26.3   6.70 0.030``0`\
`#> ``20``  ``3``162     1   ``59``905.  -``53``634.``  319.  626.  1.99     212.    4.50  4.51 0.039``6`\
`#> ``# ℹ 5 more variables: depth <dbl>, landform <fct>, .part1 <chr>, .part2 <chr>,`\
`#> ``#   .part3 <chr>`

Use the new rep_kfold partition in the gam model

\
`esm_gam_t3`` ``<-`` `[`esm_gam`](https://sjevelazco.github.io/flexsdm/reference/esm_gam.md)`(`\
`  data ``=`` ``abies2``,`\
`  response ``=`` ``"pr_ab"``,`\
`  predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"aet"``,`\
`    ``"cwd"``,`\
`    ``"tmin"``,`\
`    ``"ppt_djf"`\
`  ``)``,`\
`  partition ``=`` ``".part"``,`\
`  thr ``=`` ``NULL`\
`)`\
`#>   |                                                                              |                                                                      |   0%  |                                                                              |============                                                          |  17%  |                                                                              |=======================                                               |  33%  |                                                                              |===================================                                   |  50%  |                                                                              |===============================================                       |  67%  |                                                                              |==========================================================            |  83%  |                                                                              |======================================================================| 100%`
