# 

title: ‘flexsdm: Red Fir example’ output: rmarkdown::html_vignette
vignette: \> % % % —

## Example of full modeling process

### Study species & overview of methods

Here, we used the *flexsdm* package to model the current distribution of
California red fir (*Abies magnifica*). Red fir is a high-elevation
conifer species that’s geographic range extends through the Sierra
Nevada in California, USA, into the southern portion of the Cascade
Range of Oregon. For this species, we used presence data compiled from
several public datasets curated by natural resources agencies. We built
the distribution models using four hydro-climatic variables: actual
evapotranspiration, climatic water deficit, maximum temperature of the
warmest month, and minimum temperature of the coldest month. All
variables were resampled (aggregated) to a 1890 m spatial resolution to
improve processing time.

### Delimit of a calibration area

Delimiting the calibration area (aka accessible area) is an essential
step in SDMs both in methodological and theoretical terms. The
calibration area will affect several characteristics of a SDM like the
range of environmental variables, the number of absences, the
distribution of background points and pseudo-absences, and
unfortunately, some performance metrics like AUC and TSS. There are
several ways to delimit a calibration area. In
[calib_area()](https://sjevelazco.github.io/flexsdm/reference/calib_area.html).
We used a method that the calibration area is delimited by a 100-km
buffer around presences (shown in the figure below).

\
`# devtools::install_github('sjevelazco/flexsdm')`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`flexsdm`](https://sjevelazco.github.io/flexsdm/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`terra`](https://rspatial.org/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`\
\
`somevar`` ``<-`` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"external/somevar.tif"``, package ``=`` ``"flexsdm"``)`\
`somevar`` ``<-`` ``terra``::`[`rast`](https://rspatial.github.io/terra/reference/rast.html)`(``somevar``)`` ``# environmental data`\
[`names`](https://rspatial.github.io/terra/reference/names.html)`(``somevar``)`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"aet"``, ``"cwd"``, ``"tmx"``, ``"tmn"``)`\
[`data`](https://rdrr.io/r/utils/data.html)`(``abies``)`\
`abies_p`` ``<-`` ``abies`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``x``, ``y``, ``pr_ab``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``pr_ab`` ``==`` ``1``)`` ``# filter only for presence locations`\
\
`ca`` ``<-`\
`  `[`calib_area`](https://sjevelazco.github.io/flexsdm/reference/calib_area.md)`(`\
`    data ``=`` ``abies_p``,`\
`    x ``=`` ``"x"``,`\
`    y ``=`` ``"y"``,`\
`    method ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"buffer"``, width ``=`` ``100000``)``,`\
`    crs ``=`` `[`crs`](https://rspatial.github.io/terra/reference/crs.html)`(``somevar``)`\
`  ``)`` ``# create a calibration area with 100 km buffer around occurrence points`\
\
`# visualize the species occurrences`\
`layer1`` ``<-`` ``somevar``[[``1``]``]`\
`layer1``[``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``layer1``)``]`` ``<-`` ``1`\
\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``layer1``, col ``=`` ``"gray80"``, legend ``=`` ``FALSE``, axes ``=`` ``FALSE``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(`[`crop`](https://rspatial.github.io/terra/reference/crop.html)`(``ca``, ``layer1``)``, add ``=`` ``TRUE``)`\
[`points`](https://rspatial.github.io/terra/reference/lines.html)`(``abies_p``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"x"``, ``"y"``)``]``, col ``=`` ``"#00000480"``)`

![](v04_Red_fir_example_files/figure-html/raw%20data-1.png)

### Occurrence filtering

Sample bias in species occurrence data has long been a recognized issue
in SDM. However, environmental filtering of observation data can improve
model predictions by reducing redundancy in environmental
(e.g. climatic) hyper-space (Varela et al. 2014). Here we will use the
function occfilt_env() to thin the red fir occurrences based on
environmental space. This function is unique to *flexsdm*, and in
contrast with other packages is able to use any number of environmental
dimensions and does not perform a PCA before filtering.

Next we apply environmental occurrence filtering using 5 bins and
display the resulting filtered occurrence data

\
`abies_p``$``id`` ``<-`` ``1``:`[`nrow`](https://rspatial.github.io/terra/reference/dimensions.html)`(``abies_p``)`` ``# adding unique id to each row`\
`abies_pf`` ``<-`` ``abies_p`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`occfilt_env`](https://sjevelazco.github.io/flexsdm/reference/occfilt_env.md)`(`\
`    data ``=`` ``.``,`\
`    x ``=`` ``"x"``,`\
`    y ``=`` ``"y"``,`\
`    id ``=`` ``"id"``,`\
`    nbins ``=`` ``5``,`\
`    env_layer ``=`` ``somevar`\
`  ``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`left_join`](https://dplyr.tidyverse.org/reference/mutate-joins.html)`(``abies_p``, by ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"id"``, ``"x"``, ``"y"``)``)`\
`#> Extracting values from raster ...`\
`#> 27 records were removed because they have NAs for some variables`\
`#> Number of unfiltered records: 673`\
`#> Number of filtered records: 94`\
\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``layer1``, col ``=`` ``"gray80"``, legend ``=`` ``FALSE``, axes ``=`` ``FALSE``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(`[`crop`](https://rspatial.github.io/terra/reference/crop.html)`(``ca``, ``layer1``)``, add ``=`` ``TRUE``)`\
[`points`](https://rspatial.github.io/terra/reference/lines.html)`(``abies_p``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"x"``, ``"y"``)``]``, col ``=`` ``"#00000480"``)`\
[`points`](https://rspatial.github.io/terra/reference/lines.html)`(``abies_pf``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"x"``, ``"y"``)``]``, col ``=`` ``"#5DC86180"``)`

![](v04_Red_fir_example_files/figure-html/occurrence%20filtering-1.png)

### Block partition with 4 folds

Data partitioning, or splitting data into testing and training groups,
is a key step in building SDMs. *flexsdm* offers multiple options for
data partitioning and here we use a spatial block method. Geographically
structured data partitioning methods are especially useful if users want
to evaluate model transferability to different regions or time periods.
The part_sblock() function explores spatial blocks with different raster
cells sizes and returns the one that is best suited for the input datset
based on spatial autocorrelation, environmental similarity, and the
number of presence/absence records in each block partition. The
function’s output provides users with 1) a tibble with presence/absence
locations and the assigned partition number, 2) a tibble with
information about the best partition, and 3) a SpatRaster showing the
selected grid. Here we want to divide the data into 4 different
partitions using the spatial block method.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``10``)`\
`occ_part`` ``<-`` ``abies_pf`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`part_sblock`](https://sjevelazco.github.io/flexsdm/reference/part_sblock.md)`(`\
`    data ``=`` ``.``,`\
`    env_layer ``=`` ``somevar``,`\
`    pr_ab ``=`` ``"pr_ab"``,`\
`    x ``=`` ``"x"``,`\
`    y ``=`` ``"y"``,`\
`    n_part ``=`` ``4``,`\
`    min_res_mult ``=`` ``3``,`\
`    max_res_mult ``=`` ``200``,`\
`    num_grids ``=`` ``10``,`\
`    prop ``=`` ``1`\
`  ``)`\
`#> The following grid cell sizes will be tested:`\
`#> 5670 | 47040 | 88410 | 129780 | 171150 | 212520 | 253890 | 295260 | 336630 | 378000`\
`#> Creating basic raster mask...`\
`#> Searching for the optimal grid size...`\
`abies_pf`` ``<-`` ``occ_part``$``part`\
\
`# Transform best block partition to a raster layer with same resolution and extent than`\
`# predictor variables`\
`block_layer`` ``<-`` `[`get_block`](https://sjevelazco.github.io/flexsdm/reference/get_block.md)`(``env_layer ``=`` ``somevar``, best_grid ``=`` ``occ_part``$``grid``)`\
\
`cl`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#64146D"``, ``"#9E2962"``, ``"#F47C15"``, ``"#FCFFA4"``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``block_layer``, col ``=`` ``cl``, legend ``=`` ``FALSE``, axes ``=`` ``FALSE``)`\
[`points`](https://rspatial.github.io/terra/reference/lines.html)`(``abies_pf``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"x"``, ``"y"``)``]``)`

![](v04_Red_fir_example_files/figure-html/block%20partition-1.png)

\
\
\
`# Number of presences per block`\
`abies_pf`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  ``dplyr``::`[`group_by`](https://dplyr.tidyverse.org/reference/group_by.html)`(``.part``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  ``dplyr``::`[`count`](https://dplyr.tidyverse.org/reference/count.html)`(``)`\
`#> ``# A tibble: 4 × 2`\
`#> ``# Groups:   .part [4]`\
`#>   .part     n`\
`#>   ``<int>`` ``<int>`\
`#> ``1``     1    32`\
`#> ``2``     2    13`\
`#> ``3``     3    35`\
`#> ``4``     4    14`\
`# Additional information of the best block`\
`occ_part``$``best_part_info`\
`#> ``# A tibble: 1 × 5`\
`#>   n_grid cell_size spa_auto env_sim  sd_p`\
`#>    ``<int>``     ``<dbl>``    ``<dbl>``   ``<dbl>`` ``<dbl>`\
`#> ``1``      7    ``295``260    0.409    188.  11.6`

### Pseudo-absence/background points (using partition previously created as a mask)

In this example, we only have species presence data. However, most SDM
methods require either pseudo-absence or background data. Here, we will
use the spatial block partition we just created to generate
pseudo-absence and background points.

\
`# Spatial blocks where species occurs`\
`# Sample background points throughout study area with random method, allocating 10X the number of presences a background`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``10``)`\
`bg`` ``<-`` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``1``:``4``, ``function``(``x``)`` ``{`\
`  `[`sample_background`](https://sjevelazco.github.io/flexsdm/reference/sample_background.md)`(`\
`    data ``=`` ``abies_pf``,`\
`    x ``=`` ``"x"``,`\
`    y ``=`` ``"y"``,`\
`    n ``=`` `[`sum`](https://rdrr.io/r/base/sum.html)`(``abies_pf``$``.part`` ``==`` ``x``)`` ``*`` ``10``,`\
`    method ``=`` ``"random"``,`\
`    rlayer ``=`` ``block_layer``,`\
`    maskval ``=`` ``x``,`\
`    calibarea ``=`` ``ca`\
`  ``)`\
`}``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`bind_rows`](https://dplyr.tidyverse.org/reference/bind_rows.html)`(``)`\
`bg`` ``<-`` `[`sdm_extract`](https://sjevelazco.github.io/flexsdm/reference/sdm_extract.md)`(``data ``=`` ``bg``, x ``=`` ``"x"``, y ``=`` ``"y"``, env_layer ``=`` ``block_layer``)`\
\
`# Sample a number of pseudo-absences equal to the presence in each partition`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``10``)`\
`psa`` ``<-`` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``1``:``4``, ``function``(``x``)`` ``{`\
`  `[`sample_pseudoabs`](https://sjevelazco.github.io/flexsdm/reference/sample_pseudoabs.md)`(`\
`    data ``=`` ``abies_pf``,`\
`    x ``=`` ``"x"``,`\
`    y ``=`` ``"y"``,`\
`    n ``=`` `[`sum`](https://rdrr.io/r/base/sum.html)`(``abies_pf``$``.part`` ``==`` ``x``)``,`\
`    method ``=`` ``"random"``,`\
`    rlayer ``=`` ``block_layer``,`\
`    maskval ``=`` ``x``,`\
`    calibarea ``=`` ``ca`\
`  ``)`\
`}``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`bind_rows`](https://dplyr.tidyverse.org/reference/bind_rows.html)`(``)`\
`psa`` ``<-`` `[`sdm_extract`](https://sjevelazco.github.io/flexsdm/reference/sdm_extract.md)`(``data ``=`` ``psa``, x ``=`` ``"x"``, y ``=`` ``"y"``, env_layer ``=`` ``block_layer``)`\
\
`cl`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#280B50"``, ``"#9E2962"``, ``"#F47C15"``, ``"#FCFFA4"``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``block_layer``, col ``=`` ``"gray80"``, legend ``=`` ``FALSE``, axes ``=`` ``FALSE``)`\
[`points`](https://rspatial.github.io/terra/reference/lines.html)`(``bg``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"x"``, ``"y"``)``]``, col ``=`` ``cl``[``bg``$``.part``]``, cex ``=`` ``0.8``)`` ``# Background points`\
[`points`](https://rspatial.github.io/terra/reference/lines.html)`(``psa``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"x"``, ``"y"``)``]``, bg ``=`` ``cl``[``psa``$``.part``]``, cex ``=`` ``0.8``, pch ``=`` ``21``)`` ``# Pseudo-absences`

![](v04_Red_fir_example_files/figure-html/pseudo/absence%20and%20background%20data-1.png)

\
\
`# Bind a presences and pseudo-absences`\
`abies_pa`` ``<-`` `[`bind_rows`](https://dplyr.tidyverse.org/reference/bind_rows.html)`(``abies_pf``, ``psa``)`\
`abies_pa`` ``# Presence-Pseudo-absence database`\
`#> ``# A tibble: 188 × 4`\
`#>           x        y pr_ab .part`\
`#>       ``<dbl>``    ``<dbl>`` ``<dbl>`` ``<dbl>`\
`#> `` 1``  -``12``558.``   ``68``530.     1     3`\
`#> `` 2``  ``115``217. -``145``937.``     1     1`\
`#> `` 3``    ``3``634.   ``22``501.     1     3`\
`#> `` 4``   ``44``972.  -``60``781.``     1     1`\
`#> `` 5``  -``34``463.``  ``160``313.     1     3`\
`#> `` 6``   ``83``108.  -``27``300.``     1     1`\
`#> `` 7``  ``118``707. -``179``991.``     1     1`\
`#> `` 8``  -``49``722.``  ``141``124.     1     3`\
`#> `` 9``   ``46``612.  -``59``242.``     1     1`\
`#> ``10`` -``119``068.``  ``279``241.     1     2`\
`#> ``# ℹ 178 more rows`\
`bg`` ``# Background points`\
`#> ``# A tibble: 940 × 4`\
`#>          x        y pr_ab .part`\
`#>      ``<dbl>``    ``<dbl>`` ``<dbl>`` ``<dbl>`\
`#> `` 1``  ``46``839.  -``23``638.``     0     1`\
`#> `` 2`` -``34``431.``  -``89``788.``     0     1`\
`#> `` 3``  ``92``199. -``256``108.``     0     1`\
`#> `` 4`` ``118``659. -``282``568.``     0     1`\
`#> `` 5`` ``122``439. -``123``808.``     0     1`\
`#> `` 6`` -``60``891.``  -``50``098.``     0     1`\
`#> `` 7``  ``20``379.  -``46``318.``     0     1`\
`#> `` 8`` ``158``349. -``165``388.``     0     1`\
`#> `` 9`` ``133``779. -``259``888.``     0     1`\
`#> ``10``  ``50``619. -``218``308.``     0     1`\
`#> ``# ℹ 930 more rows`

Extract environmental data for the presence-absence and background data
. View the distributions of present points, pseudo-absence points, and
background points using the blocks as a reference map.

\
`abies_pa`` ``<-`` ``abies_pa`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`sdm_extract`](https://sjevelazco.github.io/flexsdm/reference/sdm_extract.md)`(`\
`    data ``=`` ``.``,`\
`    x ``=`` ``"x"``,`\
`    y ``=`` ``"y"``,`\
`    env_layer ``=`` ``somevar``,`\
`    filter_na ``=`` ``TRUE`\
`  ``)`\
`bg`` ``<-`` ``bg`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`sdm_extract`](https://sjevelazco.github.io/flexsdm/reference/sdm_extract.md)`(`\
`    data ``=`` ``.``,`\
`    x ``=`` ``"x"``,`\
`    y ``=`` ``"y"``,`\
`    env_layer ``=`` ``somevar``,`\
`    filter_na ``=`` ``TRUE`\
`  ``)`

### Fit models with tune_max, fit_gau, and fit_glm

Now, fit our models. The *flexsdm* package offers a wide range of
modeling options, from traditional statistical methods like GLMs and
GAMs, to machine learning methods like random forests and support vector
machines. For each modeling method, *flexsdm* provides both fit\_ and
tune\_ functions, which allow users to use default settings or adjust
hyperparameters depending on their research goals. Here, we will test
out tune_max() (tuned Maximum Entropy model), fit_gau() (fit Guassian
Process model), and fit_glm (fit Generalized Linear Model). For each
model, we selected three threshold values to generate binary suitability
predictions: the threshold that maximizes TSS (max_sens_spec), the
threshold at which sensitivity and specificity are equal
(equal_sens_spec), and the threshold at which the Sorenson index is
highest (max_sorenson). In this example, we selected TSS as the
performance metric used for selecting the best combination of
hyper-parameter values in the tuned Maximum Entropy model.

\
`t_max`` ``<-`` `[`tune_max`](https://sjevelazco.github.io/flexsdm/reference/tune_max.md)`(`\
`  data ``=`` ``abies_pa``,`\
`  response ``=`` ``"pr_ab"``,`\
`  predictors ``=`` `[`names`](https://rspatial.github.io/terra/reference/names.html)`(``somevar``)``,`\
`  background ``=`` ``bg``,`\
`  partition ``=`` ``".part"``,`\
`  grid ``=`` `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(`\
`    regmult ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0.5``, ``1.5``, ``2.5``)``,`\
`    classes ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"l"``, ``"lq"``)`\
`  ``)``,`\
`  thr ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"max_sens_spec"``, ``"equal_sens_spec"``, ``"max_sorensen"``)``,`\
`  metric ``=`` ``"TSS"``,`\
`  clamp ``=`` ``TRUE``,`\
`  pred_type ``=`` ``"cloglog"`\
`)`\
`#> Tuning model...`\
`#> Replica number: 1/1`\
`#> Partition number: 1/4`\
`#> Partition number: 2/4`\
`#> Partition number: 3/4`\
`#> Partition number: 4/4`\
`#> Fitting best model`\
`#> Formula used for model fitting:`\
`#> ~aet + cwd + tmx + tmn + I(aet^2) + I(cwd^2) + I(tmx^2) + I(tmn^2) - 1`\
`#> Replica number: 1/1`\
`#> Partition number: 1/4`\
`#> Partition number: 2/4`\
`#> Partition number: 3/4`\
`#> Partition number: 4/4`\
`f_gau`` ``<-`` `[`fit_gau`](https://sjevelazco.github.io/flexsdm/reference/fit_gau.md)`(`\
`  data ``=`` ``abies_pa``,`\
`  response ``=`` ``"pr_ab"``,`\
`  predictors ``=`` `[`names`](https://rspatial.github.io/terra/reference/names.html)`(``somevar``)``,`\
`  partition ``=`` ``".part"``,`\
`  thr ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"max_sens_spec"``, ``"equal_sens_spec"``, ``"max_sorensen"``)`\
`)`\
`#> Replica number: 1/1`\
`#> Partition number: 1/4`\
`#> Partition number: 2/4`\
`#> Partition number: 3/4`\
`#> Partition number: 4/4`\
`f_glm`` ``<-`` `[`fit_glm`](https://sjevelazco.github.io/flexsdm/reference/fit_glm.md)`(`\
`  data ``=`` ``abies_pa``,`\
`  response ``=`` ``"pr_ab"``,`\
`  predictors ``=`` `[`names`](https://rspatial.github.io/terra/reference/names.html)`(``somevar``)``,`\
`  partition ``=`` ``".part"``,`\
`  thr ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"max_sens_spec"``, ``"equal_sens_spec"``, ``"max_sorensen"``)``,`\
`  poly ``=`` ``2`\
`)`\
`#> Formula used for model fitting:`\
`#> pr_ab ~ aet + cwd + tmx + tmn + I(aet^2) + I(cwd^2) + I(tmx^2) + I(tmn^2)`\
`#> Replica number: 1/1`\
`#> Partition number: 1/4`\
`#> Partition number: 2/4`\
`#> Partition number: 3/4`\
`#> Partition number: 4/4`

### Fit an ensemble model

Spatial predictions from different SDM algorithms can vary
substantially, and ensemble modeling has become increasingly popular.
With the fit_ensemble() function, users can easily produce an ensemble
SDM based on any of the individual fit\_ and tune\_ models included the
package. In this example, we fit an ensemble model for red fir based on
the weighted average of the three individual models. We used the same
threshold values and performance metric that were implemented in the
individual models.

\
`ens_m`` ``<-`` `[`fit_ensemble`](https://sjevelazco.github.io/flexsdm/reference/fit_ensemble.md)`(`\
`  models ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``t_max``, ``f_gau``, ``f_glm``)``,`\
`  ens_method ``=`` ``"meanw"``,`\
`  thr ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"max_sens_spec"``, ``"equal_sens_spec"``, ``"max_sorensen"``)``,`\
`  thr_model ``=`` ``"max_sens_spec"``,`\
`  metric ``=`` ``"TSS"`\
`)`\
`#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%`\
`ens_m``$``performance`\
`#> ``# A tibble: 3 × 33`\
`#>   model threshold      thr_value n_presences n_absences TPR_mean TPR_sd TNR_mean`\
`#>   ``<chr>`` ``<chr>``              ``<dbl>``       ``<int>``      ``<int>``    ``<dbl>``  ``<dbl>``    ``<dbl>`\
`#> ``1`` meanw equal_sens_sp…     0.483          94         94    0.777 0.113     0.777`\
`#> ``2`` meanw max_sens_spec      0.204          94         94    0.848 0.069``9``    0.764`\
`#> ``3`` meanw max_sorensen       0.204          94         94    0.977 0.029``7``    0.594`\
`#> ``# ℹ 25 more variables: TNR_sd <dbl>, W_TPR_TNR_mean <dbl>, W_TPR_TNR_sd <dbl>,`\
`#> ``#   SORENSEN_mean <dbl>, SORENSEN_sd <dbl>, JACCARD_mean <dbl>,`\
`#> ``#   JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>, OR_mean <dbl>, OR_sd <dbl>,`\
`#> ``#   TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>, KAPPA_sd <dbl>,`\
`#> ``#   MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,`\
`#> ``#   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,`\
`#> ``#   IMAE_mean <dbl>, IMAE_sd <dbl>`

The output of *flexsdm* model objects allows you to easily compare
metrics across models, such as AUC or TSS. For example, we can use the
sdm_summarize() function to merge model performance tables.

\
`model_perf`` ``<-`` `[`sdm_summarize`](https://sjevelazco.github.io/flexsdm/reference/sdm_summarize.md)`(`[`list`](https://rdrr.io/r/base/list.html)`(``t_max``, ``f_gau``, ``f_glm``, ``ens_m``)``)`\
`model_perf`\
`#> ``# A tibble: 10 × 36`\
`#>    model_ID model threshold     thr_value n_presences n_absences TPR_mean TPR_sd`\
`#>       ``<int>`` ``<chr>`` ``<chr>``             ``<dbl>``       ``<int>``      ``<int>``    ``<dbl>``  ``<dbl>`\
`#> `` 1``        1 max   max_sens_spec     0.408          94         94    0.884 0.104 `\
`#> `` 2``        2 gau   equal_sens_s…     0.589          94         94    0.754 0.092``8`\
`#> `` 3``        2 gau   max_sens_spec     0.553          94         94    0.848 0.069``9`\
`#> `` 4``        2 gau   max_sorensen      0.553          94         94    0.906 0.077``9`\
`#> `` 5``        3 glm   equal_sens_s…     0.595          94         94    0.724 0.103 `\
`#> `` 6``        3 glm   max_sens_spec     0.578          94         94    0.818 0.159 `\
`#> `` 7``        3 glm   max_sorensen      0.334          94         94    0.975 0.033``8`\
`#> `` 8``        4 meanw equal_sens_s…     0.483          94         94    0.777 0.113 `\
`#> `` 9``        4 meanw max_sens_spec     0.204          94         94    0.848 0.069``9`\
`#> ``10``        4 meanw max_sorensen      0.204          94         94    0.977 0.029``7`\
`#> ``# ℹ 28 more variables: TNR_mean <dbl>, TNR_sd <dbl>, W_TPR_TNR_mean <dbl>,`\
`#> ``#   W_TPR_TNR_sd <dbl>, SORENSEN_mean <dbl>, SORENSEN_sd <dbl>,`\
`#> ``#   JACCARD_mean <dbl>, JACCARD_sd <dbl>, FPB_mean <dbl>, FPB_sd <dbl>,`\
`#> ``#   OR_mean <dbl>, OR_sd <dbl>, TSS_mean <dbl>, TSS_sd <dbl>, KAPPA_mean <dbl>,`\
`#> ``#   KAPPA_sd <dbl>, MCC_mean <dbl>, MCC_sd <dbl>, AUC_mean <dbl>, AUC_sd <dbl>,`\
`#> ``#   BOYCE_mean <dbl>, BOYCE_sd <dbl>, CRPS_mean <dbl>, CRPS_sd <dbl>,`\
`#> ``#   IMAE_mean <dbl>, IMAE_sd <dbl>, regmult <dbl>, classes <fct>`

### Project the ensemble model

Next we project the ensemble model in space across the entire extent of
our environmental layer, the California Floristic Province, using the
sdm_predict() function. This function can be use to predict species
suitability across any area for species’ current or future suitability.
In this example, we only project the ensemble model with one threshold,
though users have the option to project multiple models with multiple
threshold values. Here, we also specify that we want the function to
return a SpatRast with continuous suitability values above the threshold
(con_thr = TRUE).

\
`pr_1`` ``<-`` `[`sdm_predict`](https://sjevelazco.github.io/flexsdm/reference/sdm_predict.md)`(`\
`  models ``=`` ``ens_m``,`\
`  pred ``=`` ``somevar``,`\
`  thr ``=`` ``"max_sens_spec"``,`\
`  con_thr ``=`` ``TRUE``,`\
`  predict_area ``=`` ``NULL`\
`)`\
`#> Predicting ensembles`\
\
`unconstrained`` ``<-`` ``pr_1``$``meanw``[[``1``]``]`\
[`names`](https://rspatial.github.io/terra/reference/names.html)`(``unconstrained``)`` ``<-`` ``"unconstrained"`\
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
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``unconstrained``, col ``=`` ``cl``, legend ``=`` ``FALSE``, axes ``=`` ``FALSE``)`

![](v04_Red_fir_example_files/figure-html/predict%20models-1.png)

### Constrain the model with msdm_posterior

Finally, *flexsdm* offers users function that help correct
overprediction of SDM based on occurrence records and suitability
patterns. In this example we constrained the ensemble model using the
method “occurrence based restriction”, which assumes that suitable
patches that intercept species occurrences are more likely a part of
species distributions than suitable patches that do not intercept any
occurrences. Because all methods of the msdm_posteriori() function work
with presences it is important to always use the original database
(i.e., presences that have not been spatially or environmentally
filtered). All of the methods available in the msdm_posteriori()
function are based on Mendes et al. (2020).

\
`thr_val`` ``<-`` ``ens_m``$``performance`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  ``dplyr``::`[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``threshold`` ``==`` ``"max_sens_spec"``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`pull`](https://dplyr.tidyverse.org/reference/pull.html)`(``thr_value``)`\
`m_pres`` ``<-`` `[`msdm_posteriori`](https://sjevelazco.github.io/flexsdm/reference/msdm_posteriori.md)`(`\
`  records ``=`` ``abies_p``,`\
`  x ``=`` ``"x"``,`\
`  y ``=`` ``"y"``,`\
`  pr_ab ``=`` ``"pr_ab"``,`\
`  cont_suit ``=`` ``pr_1``$``meanw``[[``1``]``]``,`\
`  method ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"obr"``)``,`\
`  thr ``=`` ``thr_val``,`\
`  buffer ``=`` ``NULL`\
`)`\
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
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``m_pres``[[``1``]``]``, col ``=`` ``cl``, legend ``=`` ``FALSE``, axes ``=`` ``FALSE``)`

![](v04_Red_fir_example_files/figure-html/constrain%20with%20msdm-1.png)
