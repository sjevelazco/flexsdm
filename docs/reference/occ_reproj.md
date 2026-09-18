# Reproject location coordinates to a new coordinate reference system

Reprojects the x/y (longitude/latitude) coordinates of a set of species
occurrence records from one coordinate reference system (CRS) to
another, and returns the results either as new columns appended to the
input data or as replacements for the original coordinate columns.

## Usage

``` r
occ_reproj(
  data,
  x = "x",
  y = "y",
  from = "EPSG:4326",
  to = NULL,
  replace_cols = FALSE
)
```

## Arguments

- data:

  data.frame or tibble. A table containing occurrence records with at
  least two columns for the x (longitude) and y (latitude) coordinates.

- x:

  character. Name of the column in `data` containing the x (longitude)
  coordinate. Default `"x"`.

- y:

  character. Name of the column in `data` containing the y (latitude)
  coordinate. Default `"y"`.

- from:

  character. CRS of the input coordinates, provided as an EPSG code
  (e.g., `"EPSG:4326"`) or any other CRS string accepted by
  [`terra::vect()`](https://rspatial.github.io/terra/reference/vect.html).
  Default `"EPSG:4326"` (WGS84 unprojected).

- to:

  character. CRS to project the coordinates to, provided as an EPSG code
  (e.g., `"EPSG:32617"`) or any other CRS string accepted by
  [`terra::project()`](https://rspatial.github.io/terra/reference/project.html).
  This argument is mandatory; the function will stop with an error if
  `to = NULL`.

- replace_cols:

  logical. If `TRUE`, the original `x` and `y` columns in `data` are
  overwritten with the reprojected coordinates. If `FALSE` (default),
  the reprojected coordinates are appended as new columns `x_new` and
  `y_new`, leaving the original columns unchanged.

## Value

A tibble (or data.frame, if `replace_cols = TRUE`) with the reprojected
coordinates. If `replace_cols = FALSE`, two new columns, `x_new` and
`y_new`, are added to `data`. If `replace_cols = TRUE`, the original `x`
and `y` columns are replaced in place and the object class of `data` is
preserved.

## Details

This function is a convenience wrapper around
[`terra::vect()`](https://rspatial.github.io/terra/reference/vect.html)
and
[`terra::project()`](https://rspatial.github.io/terra/reference/project.html)
for reprojecting tabular occurrence coordinates without needing to first
convert `data` into a spatial object. It is useful, for example, when
occurrence records are stored in geographic coordinates (e.g.,
`"EPSG:4326"`) but need to be converted to a projected CRS (e.g., a UTM
zone) for distance-based calculations or to match the CRS of
environmental raster layers used elsewhere in a species distribution
modeling workflow.

## See also

[`project`](https://rspatial.github.io/terra/reference/project.html),
[`vect`](https://rspatial.github.io/terra/reference/vect.html)

## Examples

``` r
if (FALSE) { # \dontrun{
data <- data.frame(
  species = c("sp1", "sp2", "sp3"),
  x = c(-74.1, -73.9, -74.0),
  y = c(4.6, 4.7, 4.65)
)

# Append reprojected coordinates as new columns
occ_reproj(data, x = "x", y = "y", from = "EPSG:4326", to = "EPSG:32618")

# Replace original coordinate columns with reprojected values
occ_reproj(
  data,
  x = "x", y = "y",
  from = "EPSG:4326", to = "EPSG:32618",
  replace_cols = TRUE
)
} # }
```
