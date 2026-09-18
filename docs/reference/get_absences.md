# Build a presence-absence database from multi-species occurrence data

Converts a database containing occurrence records for multiple species
into a presence-absence database, either for every species present in
the data or for a specified subset of target species. For each target
species, its own records are coded as presences and the records of all
other species are recoded as absences for that species.

## Usage

``` r
get_absences(
  data,
  x = "x",
  y = "y",
  species = "species",
  target_species = NULL,
  pr_ab_name = "pr_ab"
)
```

## Arguments

- data:

  data.frame or tibble. A table with occurrence records for one or more
  species, including columns for species identity and x/y coordinates.

- x:

  character. Name of the column in `data` with the x (longitude)
  coordinate. Default `"x"`.

- y:

  character. Name of the column in `data` with the y (latitude)
  coordinate. Default `"y"`.

- species:

  character. Name of the column in `data` with the species
  name/identifier. Default `"species"`.

- target_species:

  character vector. One or more species names to build presence-absence
  data for. If `NULL` (default), a presence-absence database is built
  for every unique species found in `data[[species]]`.

- pr_ab_name:

  character. Name to give the presence-absence column in the output,
  coded `1` for presences and `0` for absences. Default `"pr_ab"`.

## Value

A tibble with columns `species`, `x`, `y` (using the names supplied in
the `species`, `x`, and `y` arguments) and the presence-absence column
named according to `pr_ab_name`. For each species in `target_species`
(or each unique species in `data`, if `target_species = NULL`), the
output contains one presence-absence block: records belonging to that
species are coded `1`, and records of all other species are recoded as
belonging to that species with a value of `0`. Blocks for all target
species are stacked row-wise, so the returned tibble has
`nrow(data) * length(sp)` rows, where `sp` is the set of target species.

## Details

Only the `species`, `x`, and `y` columns are retained in the output; any
other columns present in `data` (e.g., environmental covariates or
metadata) are dropped. If you need those columns downstream, re-join
them after calling this function.

Because absences for a given target species are drawn from every other
species' occurrence records, duplicate x/y coordinates across species in
the input will produce duplicate absence coordinates in the output.

## Examples

``` r
if (FALSE) { # \dontrun{
data <- data.frame(
  species = c("sp1", "sp1", "sp2", "sp2", "sp3"),
  x = c(-74.1, -74.2, -73.9, -73.8, -74.0),
  y = c(4.6, 4.65, 4.7, 4.72, 4.68)
)

# Presence-absence database for every species in the data
get_absences(data, x = "x", y = "y", species = "species", target_species = NULL, pr_ab_name = "pr_ab")

# Presence-absence database for a single target species
get_absences(data, x = "x", y = "y", species = "species", target_species = "sp1", pr_ab_name = "pr_ab")

# Presence-absence database for a subset of species, with a custom
# presence-absence column name
get_absences(data, x = "x", y = "y", species = "species", target_species = c("sp1", "sp2"), pr_ab_name = "occ")
} # }
```
