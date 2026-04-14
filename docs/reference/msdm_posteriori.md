# Methods to correct overprediction of species distribution models based on occurrences and suitability patterns.

These methods reduce overprediction of species distribution models based
on a posteriori methods (see Mendes et al 2020), i.e., the combination
of the patterns of species occurrences and predicted suitability

## Usage

``` r
msdm_posteriori(
  records,
  x,
  y,
  pr_ab,
  cont_suit,
  method = c("obr", "pres", "lq", "mcp", "bmcp"),
  pres_as_patch = FALSE,
  thr = "equal_sens_spec",
  con_thr = FALSE,
  buffer = NULL,
  crs = NULL
)
```

## Arguments

- records:

  tibble or data.frame. A database with spatial coordinates of species
  presences and absences (or pseudo-absence) used to create species
  distribution models.

- x:

  character. Column name with spatial x coordinates.

- y:

  character. Column name with spatial y coordinates.

- pr_ab:

  character. Column name with presence and absence data (i.e. 1 and 0)

- cont_suit:

  SpatRaster. Raster with continuous suitability predictions
  "species_specific" type calculates the minimum pairwise-distances
  between all occurrences and then selects the maximum distance, i.e.,
  the value of the buffer will be the maximum distance from the minimum
  distance. This procedure depends on the spatial pattern of each
  species' occurrences; thus, for each species, a value of buffer width
  will be calculated (usage buffer="species_specific").

- method:

  character. A character string indicating which constraint method will
  be used (see in details).

- pres_as_patch:

  character. For the 'lq' method, assume that cells with presences but
  below the threshold are patches. Default FALSE.

- thr:

  character or numeric. Threshold used to get binary suitability values
  (i.e. 0,1), needed for threshold-dependent performance metrics. Only
  one threshold type can be specified. It is necessary to provide a
  vector for this argument. The following threshold criteria are
  available:

  - lpt: The highest threshold at which there is no omission.

  - equal_sens_spec: Threshold at which the sensitivity and specificity
    are equal.

  - max_sens_spec: Threshold at which the sum of the sensitivity and
    specificity is the highest (aka threshold that maximizes the TSS).

  - max_jaccard: The threshold at which the Jaccard index is the
    highest.

  - max_sorensen: The threshold at which the Sorensen index is highest.

  - max_fpb: The threshold at which FPB is highest.

  - sensitivity: Threshold based on a specified sensitivity value. Usage
    thr = c('sensitivity', sens='0.6') or thr = c('sensitivity'). 'sens'
    refers to sensitivity value. If it is not specified a sensitivity
    values, function will use by default 0.9

  - Also, it is possible to specify the threshold value using a numeric
    value (thr = 0.623)

    Default "equal_sens_spec".

- con_thr:

  logical. If TRUE, returns continuous suitability values for cells
  above the threshold, with other cells as 0. If False, returns a binary
  map (1 for above threshold, 0 for below). Default FALSE.

- buffer:

  numeric. Buffer width use in 'bmcp' approach. The buffer width will be
  interpreted in m if Coordinate reference system used in "crs" argument
  has a longitude/latitude, or map units in other cases. Usage
  buffer=50000. Default NULL

- crs:

  character. Coordinate reference system used for calculating buffer in
  "bmcp" method.

## Value

This function return a SpatRaster with continuous and binary prediction.

## Details

These function help reduce overprediction of species distribution models
based on the combination of the patterns of species occurrences and
predicted suitability. It is recommended to use these approaches only
for current distribution not for models projected for different time
periods (past or future).

Five methods are implemented:

Abbreviation list

- SDM: species distribution model

- l: suitability patches that intercept species occurrences

- k: suitability patches that do not intercept species occurrences

- T: threshold distances used to select suitability patches

These methods reduce overprediction of species distribution models
already fitted based on the occurrences and suitability patterns of
species (see 'thr' arguments)

Method 'obr' (Occurrences Based Restriction). This method assumes that
suitable patches intercepting species occurrences (l) are more likely to
be part of species distributions than suitable patches that do not
intercept any occurrence (k). Distance from all k patches to the closest
l patch is calculated, then k patches are removed that exceed a
species-specific distance threshold from SDMs models. This threshold (T)
is calculated as the maximum distance in a vector of minimum pairwise
distances between occurrences. Whenever a suitable pixel is within a k
patch that is more than distance T from the closest l patch, the
suitability of the pixel is reduced to zero. It is assumed that this
simple threshold is a surrogate for the species-specific dispersal
ability. If T is low, either the species has been sampled throughout its
distribution, or the species is geographically restricted, justifying a
narrow inclusion of k patches (Mendes et al., 2020).

Method 'lq' (Lower Quantile). This method is similar to the 'obr'
method, except by the procedure to define a distance threshold to
withdrawn k patches, which is the lower quartile distance between k
patches to the closest l patch (i.e., threshold distance is based on
patches with and without presences and not in presences points as in
'obr' method) . Whenever a suitable pixel is within a k patch, i.e., not
within this lower quartile, the suitability of the pixel is reduced to
zero. This means that 75% of k patches were withdrawn from the model
(Mendes et al., 2020).

Method 'pres' (only occurrences based restriction). This is a more
restrictive variant of the 'obr' method. It only retains those pixels in
suitability patches intercepting occurrences (k) (Mendes et al., 2020).

Method 'mcp' (Minimum Convex Polygon). Compiled and adapted from Kremen
et al. (2008), this method excludes from SDM predictions suitable pixels
that do not intercept a minimum convex polygon, with interior angles
smaller than 180, enclosing all occurrences of a species.

Method 'bmcp' (Buffered Minimum Convex Polygon). Compiled and adapted
from Kremen et al. (2008), it is similar to the 'mcp' method except for
the inclusion of a buffer zone surrounding minimum convex polygons. For
this method a buffer width value must be provided in "buffer" argument
and CRS in "crs" argument.

For further methodological and performance information of these methods
see Mendes et al. (2020).

If using one these constraining methods, cite Mendes et al (2020).

## References

- Mendes, P.; Velazco S.J.E.; Andrade, A.F.A.; De Marco, P. (2020)
  Dealing with overprediction in species distribution models: how adding
  distance constraints can improve model accuracy, Ecological Modelling,
  in press. https://doi.org/10.1016/j.ecolmodel.2020.109180

- Kremen, C., Cameron, A., Moilanen, A., Phillips, S. J., Thomas, C. D.,
  Beentje, H., . Zjhra, M. L. (2008). Aligning Conservation Priorities
  Across Taxa in Madagascar with High-Resolution Planning Tools.
  Science, 320(5873), 222-226. doi:10.1126/science.1155193

## See also

[`msdm_priori`](https://sjevelazco.github.io/flexsdm/reference/msdm_priori.md)

## Examples

``` r
if (FALSE) { # \dontrun{
require(dplyr)
require(terra)

data("spp")
somevar <- system.file("external/somevar.tif", package = "flexsdm")
somevar <- terra::rast(somevar)


# Preparing data for modeling a species
set.seed(10)
occ <- spp %>%
  dplyr::filter(species == "sp2") %>% # filter a species
  sdm_extract(
    data = ., x = "x", y = "y",
    env_layer = somevar, filter_na = TRUE
  ) %>% # extrac variables values
  part_random(.,
    pr_ab = "pr_ab",
    method = c(method = "kfold", folds = 10)
  ) # add columns with partition

# Fit a model
m_glm <- fit_glm(
  data = occ,
  response = "pr_ab",
  predictors = names(somevar),
  partition = ".part",
  thr = "equal_sens_spec",
)


# Lets predict this model
m_pred <- sdm_predict(models = m_glm, pred = somevar, thr = NULL, con_thr = FALSE)
plot(m_pred[[1]])
m_pred[[1]] %>% plot()

# Lets extract the raster from this list
m_pred <- m_pred[[1]]

### bmcp method
m_bmcp <- msdm_posteriori(
  records = occ,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  method = "bmcp",
  cont_suit = m_pred,
  thr = "equal_sens_spec",
  buffer = 30000,
  crs = crs(m_pred)
)

plot(m_bmcp)


### mcp method
m_mcp <- msdm_posteriori(
  records = occ,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  method = "mcp",
  cont_suit = m_pred,
  thr = "equal_sens_spec",
  buffer = NULL
)

plot(m_mcp)


### pres method
m_pres <- msdm_posteriori(
  records = occ,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  method = "pres",
  cont_suit = m_pred,
  thr = "equal_sens_spec",
  buffer = NULL
)

plot(m_pres)


### lq method
m_lq <- msdm_posteriori(
  records = occ,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  method = "lq",
  cont_suit = m_pred,
  thr = "equal_sens_spec",
  buffer = NULL
)

plot(m_lq)


### obr method
m_obr <- msdm_posteriori(
  records = occ,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  method = "obr",
  cont_suit = m_pred,
  thr = "equal_sens_spec",
  buffer = NULL
)

plot(m_obr)
} # }
```
